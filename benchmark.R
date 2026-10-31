#!/usr/bin/env Rscript
# Benchmark: future (multisession) vs mirai vs mizu across four
# representative tidymodels workloads, against a sequential lapply()
# baseline.
#
#   Rscript benchmark.R                 full run (~25 min)
#   BENCH_QUICK=1 Rscript benchmark.R   smoke test (tiny sizes)
#   BENCH_SWEEP=1 Rscript benchmark.R   worker-count sweep (2/4/8/12) on
#                                       workloads B and C; writes
#                                       results/sweep.rds and exits
#   bash run_all.sh                     resumable driver: one process per
#                                       workload, skips completed phases
#                                       (robust to killed/interrupted runs)
#
# Writes results/timings.rds, results/monitors.rds, results/meta.rds.

quick <- Sys.getenv("BENCH_QUICK") == "1"
sweep <- Sys.getenv("BENCH_SWEEP") == "1"
# BENCH_ONLY=<workload id> runs a single workload and merges into the
# existing results files (resume after an interrupted full run).
only <- Sys.getenv("BENCH_ONLY")

note <- function(...) {
  cat(...)
  flush.console()
}

suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(rsample))

dir.create("results", showWarnings = FALSE)

workers <- min(parallelly::availableCores(), 12)
n_timed <- if (quick) 1L else 3L

# Not required with current future: since 1.40.0 the multisession
# backend imposes no limit on globals size by default (pre-1.40 the
# default was 500 MiB per future; 1.75.0 made +Inf the default on
# every code path). This setting is a harmless remnant of that era.
options(future.globals.maxSize = 8 * 1024^3)

# ------------------------------------------------------------------------------
# Data

# Workload A: regression, as in topepo/bootstrap-parallel-benchmark
set.seed(1)
dat_reg <- tibble(y = rnorm(10000)) |>
  mutate(.pred = y + rnorm(10000, sd = 0.1))

# Workload B: two-class classification, 20 predictors
set.seed(2)
n <- 10000
x <- matrix(rnorm(n * 20), ncol = 20)
colnames(x) <- paste0("x", 1:20)
beta <- c(1.2, -0.9, 0.8, 0.7, -0.6, rep(0, 15))
eta <- as.vector(x %*% beta) + 0.9 * x[, 1] * x[, 2] + rnorm(n, sd = 1.5)
dat_cls <- as_tibble(x) |>
  mutate(
    class = factor(
      ifelse(plogis(eta) > 0.5, "yes", "no"),
      levels = c("yes", "no")
    )
  )

# Workload C: held-out validation set for bootstrap permutation
# importance, 12k rows x 30 predictors; each task carries one bootstrap
# replicate (~2.8 MB) and permutes every predictor in turn
set.seed(11)
nv <- if (quick) 3000L else 12000L
np <- 30L
xm <- matrix(rnorm(nv * np), ncol = np)
colnames(xm) <- paste0("v", 1:np)
b <- c(1.1, -0.9, 0.8, -0.7, 0.6, 0.5, -0.4, 0.3, rep(0, np - 8))
eta_v <- as.vector(xm %*% b) + rnorm(nv, sd = 2)
dat_valid <- as_tibble(xm) |>
  mutate(
    y = factor(
      ifelse(plogis(eta_v) > 0.5, "yes", "no"),
      levels = c("no", "yes")
    )
  )
# Strip the fit to prediction essentials (butcher-style): a glm fit
# carries its training frame, QR matrix, and per-observation vectors
# that predict() never uses — shipping them would multiply every
# transfer several-fold.
fit_perm <- glm(
  y ~ .,
  data = dat_valid,
  family = binomial(),
  model = FALSE,
  x = FALSE,
  y = FALSE
)
fit_perm[c(
  "effects",
  "fitted.values",
  "residuals",
  "linear.predictors",
  "weights",
  "prior.weights",
  "data"
)] <- NULL
fit_perm$qr$qr <- matrix(0, 0L, 0L) # predict() needs the component, not the matrix

# Workload D: insurance-style severity model. Training data (20k x 20)
# is bootstrap-refit; each refit scores a fixed 250k-row portfolio
# (~38 MB shared read-only input) and returns the 250k predicted losses
# (~2 MB per-task result) for per-policy bootstrap intervals.
set.seed(12)
nt <- if (quick) 5000L else 20000L
xt <- matrix(rnorm(nt * 20), ncol = 20)
colnames(xt) <- paste0("x", 1:20)
bt <- c(0.8, -0.6, 0.5, 0.4, -0.3, rep(0, 15))
mu <- exp(as.vector(xt %*% bt) + rnorm(nt, sd = 0.5))
dat_train <- as_tibble(xt) |>
  mutate(loss = rgamma(nt, shape = 2, scale = mu / 2))

npf <- if (quick) 20000L else 250000L
xp <- matrix(rnorm(npf * 20), ncol = 20)
colnames(xp) <- paste0("x", 1:20)
portfolio <- as_tibble(xp)

# ------------------------------------------------------------------------------
# Task functions. Self-contained on purpose: packages are require()d inside
# so every backend pays the same worker-side setup, and no backend-specific
# globals/packages export mechanism is needed.

# A: yardstick metrics on one bootstrap split (~3 ms)
task_metrics <- function(split) {
  require(rsample)
  require(yardstick)
  mtr <- metric_set(rmse, rsq, mae, ccc)
  mtr(analysis(split), y, .pred)
}

# B: parsnip workflow (recipe + logistic regression) on one bootstrap
# split, scored on the out-of-bag assessment set (~45 ms)
task_fit_glm <- function(split) {
  require(rsample)
  require(recipes)
  require(parsnip)
  require(workflows)
  require(yardstick)
  require(dplyr)
  wf <- workflow() |>
    add_recipe(
      recipe(class ~ ., data = analysis(split)) |>
        step_normalize(all_numeric_predictors())
    ) |>
    add_model(logistic_reg() |> set_engine("glm"))
  fitted <- fit(wf, data = analysis(split))
  assess <- assessment(split)
  pred <- predict(fitted, assess, type = "prob") |>
    bind_cols(predict(fitted, assess, type = "class")) |>
    bind_cols(assess |> select(class))
  bind_rows(
    roc_auc(pred, class, .pred_yes),
    accuracy(pred, class, .pred_class)
  )
}

# C: one bootstrap replicate of the validation set; permute each
# predictor in turn and record the AUC drop against the pre-fit model
# (~85 ms). The kind is pinned because mirai daemons default to
# L'Ecuyer-CMRG while R sessions default to Mersenne-Twister.
task_permimp <- function(task, model) {
  require(rsample)
  require(yardstick)
  require(dplyr)
  set.seed(task$seed, kind = "Mersenne-Twister")
  dat <- analysis(task$split)
  truth <- dat$y
  pred <- as.numeric(predict(model, newdata = dat, type = "response"))
  base_auc <- roc_auc_vec(truth, pred, event_level = "second")
  cols <- setdiff(colnames(dat), "y")
  importance <- numeric(length(cols))
  for (j in seq_along(cols)) {
    orig <- dat[[cols[j]]]
    dat[[cols[j]]] <- orig[sample.int(length(orig))]
    p <- as.numeric(predict(model, newdata = dat, type = "response"))
    importance[j] <- base_auc - roc_auc_vec(truth, p, event_level = "second")
    dat[[cols[j]]] <- orig
  }
  tibble(col = cols, base_auc = base_auc, importance = importance)
}

# D: refit the severity model on one bootstrap split and predict
# expected loss for every policy in the shared portfolio (~90 ms);
# returns a 250k-length numeric vector
task_bootpred <- function(split, portfolio) {
  require(rsample)
  fit <- glm(
    loss ~ .,
    data = analysis(split),
    family = Gamma(link = "log"),
    control = list(maxit = 25)
  )
  as.numeric(predict(fit, newdata = portfolio, type = "response"))
}

# ------------------------------------------------------------------------------
# Workloads

set.seed(3)
rs_a <- bootstraps(dat_reg, times = if (quick) 200 else 5000)
set.seed(4)
rs_b <- bootstraps(dat_cls, times = if (quick) 20 else 250)
set.seed(13)
rs_d <- bootstraps(dat_train, times = if (quick) 8 else 400)

set.seed(14)
rs_c <- bootstraps(dat_valid, times = if (quick) 12 else 240)
tasks_c <- lapply(seq_len(nrow(rs_c)), function(i) {
  list(split = rs_c$splits[[i]], seed = 3000L + i)
})

workloads <- list(
  list(
    id = "boot_metrics",
    tasks = rs_a$splits,
    fn = task_metrics,
    args = list(),
    n_seq = if (quick) 1L else 3L
  ),
  list(
    id = "resample_fit",
    tasks = rs_b$splits,
    fn = task_fit_glm,
    args = list(),
    n_seq = if (quick) 1L else 3L
  ),
  list(
    id = "perm_importance",
    tasks = tasks_c,
    fn = task_permimp,
    args = list(model = fit_perm),
    n_seq = if (quick) 1L else 3L
  ),
  list(
    id = "boot_pred_interval",
    tasks = rs_d$splits,
    fn = task_bootpred,
    args = list(portfolio = portfolio),
    n_seq = if (quick) 1L else 3L
  )
)

if (nzchar(only)) {
  workloads <- Filter(function(w) w$id == only, workloads)
  if (length(workloads) == 0L) {
    stop("unknown BENCH_ONLY workload: ", only)
  }
}

# ------------------------------------------------------------------------------
# Backend harness. Cluster/pool/daemons are set up outside the timed region
# (as in the original benchmark); each backend gets one warmup run (also the
# correctness check), n_timed timed runs, and one syrup-monitored run.
# Constant arguments (workload D's portfolio; C's only constant is the
# fitted model) go through each backend's natural map call: future exports
# them as globals once per worker chunk, mirai serializes them with every
# task, mizu stages them once per call. mizu's map is handle-first:
# mizu_map() stages and submits, returning a mizu_map handle immediately;
# mizu_collect() blocks for the results. Both halves sit inside the timed
# region, so the measured work matches the other backends' blocking maps.

setup_backend <- function(backend, workers) {
  switch(
    backend,
    future = {
      future::plan(future::multisession, workers = workers)
      NULL
    },
    mirai = {
      mirai::daemons(workers)
      NULL
    },
    mizu = ,
    mizu_stream = mizu::mizu_pool(workers)
  )
}

run_backend <- function(backend, handle, tasks, fn, args) {
  switch(
    backend,
    future = do.call(future.apply::future_lapply, c(list(tasks, fn), args)),
    mirai = do.call(
      mirai::mirai_map,
      c(list(tasks, fn), list(.args = args))
    )[],
    mizu = mizu::mizu_collect(
      do.call(mizu::mizu_map, c(list(handle, tasks, fn), args))
    ),
    # D variant: x streams to workers in per-chunk slices under a
    # sliding window instead of being staged wholesale
    mizu_stream = mizu::mizu_collect(
      do.call(
        mizu::mizu_map,
        c(list(handle, tasks, fn), args, list(.stream = TRUE))
      )
    )
  )
}

teardown_backend <- function(backend, handle) {
  switch(
    backend,
    future = future::plan(future::sequential),
    mirai = mirai::daemons(0),
    mizu = ,
    mizu_stream = mizu::mizu_pool_stop(handle)
  )
}

# ------------------------------------------------------------------------------

timings <- list()
monitors <- list()
correct <- list()
refs <- list()
prev_meta <- list()

backends_all <- c("future", "mirai", "mizu")
backends_run <- backends_all
backends_warmup <- character(0)
seq_run <- TRUE

# An extra mizu configuration measured on workload D only (see
# run_backend): a streaming map whose per-chunk transport mirrors
# mirai's model.
wl_backends <- function(wl_id) {
  if (identical(wl_id, "boot_pred_interval")) {
    c(backends_all, "mizu_stream")
  } else {
    backends_all
  }
}

# ------------------------------------------------------------------------------
# Worker-count sweep (BENCH_SWEEP=1): workloads B and C only, timed at
# 2/4/8/12 workers per backend plus a sequential baseline. Validates the
# report's Amdahl fits: each curve is pinned by the 12-worker measurement,
# so the lower worker counts test whether the implied serial fraction is
# really constant. No syrup monitoring; writes results/sweep.rds and exits
# without touching the main results files.

if (sweep) {
  sweep_workers <- sort(unique(c(2L, 4L, 8L, workers)))
  sweep_wls <- Filter(
    function(w) w$id %in% c("resample_fit", "perm_importance"),
    workloads
  )
  sweep_timings <- list()
  sweep_correct <- list()
  record_sweep <- function(workload, backend, workers, run, tt) {
    sweep_timings[[length(sweep_timings) + 1L]] <<- tibble(
      workload = workload,
      backend = backend,
      workers = workers,
      run = run,
      user = unname(tt[["user.self"]]),
      system = unname(tt[["sys.self"]]),
      elapsed = unname(tt[["elapsed"]])
    )
  }

  for (wl in sweep_wls) {
    note(sprintf("\n=== sweep: %s (%d tasks) ===\n", wl$id, length(wl$tasks)))

    tt <- system.time(ref <- do.call(lapply, c(list(wl$tasks, wl$fn), wl$args)))
    record_sweep(wl$id, "sequential", 1L, 1L, tt)
    note(sprintf("sequential 1: %.2fs\n", tt[["elapsed"]]))
    if (wl$n_seq > 1L) {
      for (i in 2:wl$n_seq) {
        tt <- system.time(do.call(lapply, c(list(wl$tasks, wl$fn), wl$args)))
        record_sweep(wl$id, "sequential", 1L, i, tt)
        note(sprintf("sequential %d: %.2fs\n", i, tt[["elapsed"]]))
      }
    }

    for (w in sweep_workers) {
      for (be in backends_all) {
        h <- setup_backend(be, w)
        out <- run_backend(be, h, wl$tasks, wl$fn, wl$args)
        ok <- isTRUE(all.equal(unname(ref), unname(out), tolerance = 1e-8))
        sweep_correct[[paste(wl$id, be, w, sep = " / ")]] <- ok
        note(sprintf("%s w=%d warmup (correct: %s)\n", be, w, ok))
        if (!ok) {
          warning("result mismatch: ", wl$id, " / ", be, " / workers=", w)
        }
        for (i in seq_len(n_timed)) {
          tt <- system.time(run_backend(be, h, wl$tasks, wl$fn, wl$args))
          record_sweep(wl$id, be, w, i, tt)
          note(sprintf("%s w=%d %d: %.2fs\n", be, w, i, tt[["elapsed"]]))
        }
        teardown_backend(be, h)
      }
    }
  }

  saveRDS(
    list(
      timings = bind_rows(sweep_timings),
      sweep_workers = sweep_workers,
      correct = sweep_correct,
      timestamp = Sys.time()
    ),
    "results/sweep.rds"
  )
  note("\nsweep done\n")
  quit("no")
}

if (nzchar(only)) {
  # resume mode: keep rows for other workloads and for this workload's
  # already-complete backends; re-run only what is missing (a killed
  # process loses at most one backend phase). A backend with complete
  # timings but a missing correctness entry (killed before the meta
  # write) gets a warmup-only re-verification.
  wl_id <- workloads[[1]]$id
  if (file.exists("results/meta.rds")) {
    prev_meta <- readRDS("results/meta.rds")
    correct <- prev_meta$correct
  }
  if (file.exists("results/timings.rds")) {
    prev_t <- readRDS("results/timings.rds")
    prev_m <- readRDS("results/monitors.rds")
    mine <- prev_t[prev_t$workload == wl_id, ]
    seq_run <- sum(mine$backend == "sequential") < workloads[[1]]$n_seq
    backends_run <- Filter(
      function(b) sum(mine$backend == b) < n_timed + 1L,
      wl_backends(wl_id)
    )
    backends_warmup <- Filter(
      function(b) {
        sum(mine$backend == b) >= n_timed + 1L &&
          !paste(wl_id, b, sep = " / ") %in% names(correct)
      },
      wl_backends(wl_id)
    )
    complete <- setdiff(
      c("sequential", wl_backends(wl_id)),
      c(if (seq_run) "sequential", backends_run)
    )
    timings <- list(prev_t[
      prev_t$workload != wl_id | prev_t$backend %in% complete,
    ])
    monitors <- list(prev_m[
      prev_m$workload != wl_id | prev_m$backend %in% complete,
    ])
  }
}

extras_todo <- Filter(
  function(w) length(w$args) > 0L && !w$id %in% names(prev_meta$extras),
  workloads
)

if (
  nzchar(only) &&
    !seq_run &&
    length(backends_run) == 0L &&
    length(backends_warmup) == 0L &&
    length(extras_todo) == 0L
) {
  note(sprintf("=== %s already complete, skipping ===\n", workloads[[1]]$id))
  quit("no")
}

record <- function(workload, backend, run, monitored, tt) {
  timings[[length(timings) + 1L]] <<- tibble(
    workload = workload,
    backend = backend,
    run = run,
    monitored = monitored,
    user = unname(tt[["user.self"]]),
    system = unname(tt[["sys.self"]]),
    elapsed = unname(tt[["elapsed"]])
  )
}

checkpoint <- function() {
  saveRDS(bind_rows(timings), "results/timings.rds")
  saveRDS(bind_rows(monitors), "results/monitors.rds")
}

main_pid <- Sys.getpid()

for (wl in workloads) {
  note(sprintf("\n=== %s (%d tasks) ===\n", wl$id, length(wl$tasks)))
  if (nzchar(only)) {
    run_these <- intersect(wl_backends(wl$id), backends_run)
    warm_these <- intersect(wl_backends(wl$id), backends_warmup)
  } else {
    run_these <- wl_backends(wl$id)
    warm_these <- backends_warmup
  }

  ref <- do.call(lapply, c(list(wl$tasks, wl$fn), wl$args))
  refs[[wl$id]] <- ref

  for (be in warm_these) {
    h <- setup_backend(be, workers)
    out <- run_backend(be, h, wl$tasks, wl$fn, wl$args)
    ok <- isTRUE(all.equal(unname(ref), unname(out), tolerance = 1e-8))
    correct[[paste(wl$id, be, sep = " / ")]] <- ok
    note(sprintf("%s warmup re-verify (correct: %s)\n", be, ok))
    teardown_backend(be, h)
  }

  if (seq_run) {
    for (i in seq_len(wl$n_seq)) {
      tt <- system.time(do.call(lapply, c(list(wl$tasks, wl$fn), wl$args)))
      record(wl$id, "sequential", i, FALSE, tt)
      note(sprintf("sequential %d: %.2fs\n", i, tt[["elapsed"]]))
    }
  }

  for (be in run_these) {
    h <- setup_backend(be, workers)

    out <- run_backend(be, h, wl$tasks, wl$fn, wl$args)
    ok <- isTRUE(all.equal(unname(ref), unname(out), tolerance = 1e-8))
    correct[[paste(wl$id, be, sep = " / ")]] <- ok
    note(sprintf("%s warmup (correct: %s)\n", be, ok))
    if (!ok) {
      warning("result mismatch: ", wl$id, " / ", be)
    }
    # release the warmup result (a view pins its whole map region)
    rm(out)
    invisible(gc())

    for (i in seq_len(n_timed)) {
      tt <- system.time(run_backend(be, h, wl$tasks, wl$fn, wl$args))
      record(wl$id, be, i, FALSE, tt)
      note(sprintf("%s %d: %.2fs\n", be, i, tt[["elapsed"]]))
      # release discarded results (views pin their regions) before the
      # next run; untimed
      invisible(gc())
    }

    tt <- system.time(
      mon <- syrup::syrup(
        run_backend(be, h, wl$tasks, wl$fn, wl$args),
        interval = 1 / 20
      )
    )
    record(wl$id, be, n_timed + 1L, TRUE, tt)
    monitors[[length(monitors) + 1L]] <- mon |>
      filter(!is.na(pct_cpu)) |>
      mutate(across(any_of(c("rss", "vms")), as.numeric)) |>
      mutate(workload = wl$id, backend = be, main_pid = main_pid)
    note(sprintf("%s monitored: %.2fs\n", be, tt[["elapsed"]]))
    rm(mon)
    invisible(gc())

    teardown_backend(be, h)
    checkpoint()
  }
}

# ------------------------------------------------------------------------------
# Footnote measurement (workloads C and D): mirai's documented pattern for
# large shared objects is everywhere() — one out-of-band export per daemon
# instead of per-task serialization. Timed separately: the export, then one
# map run with the shared object resolved from the daemons' global env.

extras <- if (length(prev_meta) > 0L) prev_meta$extras else list()
if (length(extras_todo) > 0L) {
  mirai::daemons(workers)
  for (wl in extras_todo) {
    f <- wl$fn
    te <- system.time(
      mirai::everywhere({}, shared = wl$args)
    )[["elapsed"]]
    tm <- system.time(
      out <- mirai::mirai_map(
        wl$tasks,
        function(t) do.call(f, c(list(t), shared)),
        f = f
      )[]
    )[["elapsed"]]
    ok <- isTRUE(all.equal(
      unname(refs[[wl$id]]),
      unname(out),
      tolerance = 1e-8
    ))
    extras[[wl$id]] <- list(export = te, map = tm, correct = ok)
    note(sprintf(
      "%s mirai+everywhere: export %.2fs, map %.2fs (correct: %s)\n",
      wl$id,
      te,
      tm,
      ok
    ))
  }
  mirai::daemons(0)
}

# ------------------------------------------------------------------------------

n_tasks <- vapply(workloads, function(w) length(w$tasks), integer(1))
names(n_tasks) <- vapply(workloads, function(w) w$id, character(1))
if (length(prev_meta) > 0L) {
  n_tasks <- c(prev_meta$n_tasks, n_tasks)
  n_tasks <- n_tasks[!duplicated(names(n_tasks), fromLast = TRUE)]
}

meta <- list(
  workers = workers,
  quick = quick,
  n_tasks = n_tasks,
  correct = correct,
  extras = extras,
  r_version = R.version.string,
  sysinfo = as.list(Sys.info()[c("sysname", "machine", "nodename")]),
  cores = parallelly::availableCores(),
  packages = vapply(
    c(
      "future",
      "future.apply",
      "mirai",
      "mizu",
      "rsample",
      "yardstick",
      "parsnip",
      "workflows",
      "recipes",
      "syrup",
      "parallelly"
    ),
    function(p) as.character(packageVersion(p)),
    character(1)
  ),
  timestamp = Sys.time()
)
saveRDS(meta, "results/meta.rds")

note("\ndone\n")

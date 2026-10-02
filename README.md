# tidymodels Parallel Backend Benchmark

Benchmarks three parallel backends — future (multisession), mirai, and
mizu — across four representative tidymodels workloads, against a
sequential `lapply()` baseline:

- **A: Bootstrap metrics** — regression bootstrap with yardstick metrics
- **B: Resampled model fits** — classification model over resamples
- **C: Permutation importance** — bootstrap permutation importance on a
  held-out validation set
- **D: Bootstrap prediction intervals** — bootstrap prediction intervals

## Requirements

- R with these packages:

  ```r
  install.packages(c(
    "tidymodels", "ggplot2", "parallelly",
    "future", "future.apply", "mirai"
  ))
  pak::pak("shikokuchuo/mizu")  # not on CRAN
  ```
- [Quarto CLI](https://quarto.org/docs/get-started/) for the report.

## Running the benchmark

```sh
Rscript benchmark.R                 # full run (~25 min)
BENCH_QUICK=1 Rscript benchmark.R   # smoke test (tiny sizes)
BENCH_SWEEP=1 Rscript benchmark.R   # worker-count sweep (2/4/8/12) for
                                    # workloads B and C (results/sweep.rds)
bash run_all.sh                     # resumable driver
```

`run_all.sh` runs each workload in a separate R process and skips
completed phases, so it is safe to relaunch after an interrupted run.
To run a single workload:

```sh
BENCH_ONLY=boot_metrics Rscript benchmark.R
```

(workload ids: `boot_metrics`, `resample_fit`, `perm_importance`,
`boot_pred_interval`)

Results are written to `results/` (`timings.rds`, `monitors.rds`,
`meta.rds`; the sweep writes `sweep.rds`).

## Generating the report

After the benchmark has produced `results/`:

```sh
quarto render benchmark.qmd
```

This renders `benchmark.html`, a self-contained report with timing
summaries, speedup and efficiency plots, and resource-monitor traces.

#!/bin/bash
# Mop-up driver: benchmark.R self-skips completed phases in BENCH_ONLY
# mode, so this simply invokes each workload in turn. Safe to relaunch
# repeatedly until everything is complete.
cd "$(dirname "$0")"
for w in boot_metrics resample_fit perm_importance boot_pred_interval; do
  echo "=== DRIVER2: invoke $w $(date) ==="
  BENCH_ONLY="$w" Rscript benchmark.R
  echo "=== DRIVER2: $w exited $? $(date) ==="
done
echo "=== DRIVER2: all done $(date) ==="

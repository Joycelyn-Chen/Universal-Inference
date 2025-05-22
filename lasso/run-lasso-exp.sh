#!/usr/bin/env bash
# run_experiments.sh  — launch overnight on your Ubuntu box

echo "Starting parameter sweep at $(date)"
Rscript lasso-exp.r \
  --vanilla \
  2>&1 | tee experiment_$(date +%Y%m%d_%H%M).log
echo "Done at $(date)"

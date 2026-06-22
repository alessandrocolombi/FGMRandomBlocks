#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DEFAULT="./Simulation_Study_parallel.R"   # default script if none given
TAIL_LINES="${TAIL_LINES:-200}"           # lines to show before follow

SCRIPT="${1:-$SCRIPT_DEFAULT}"            # allow: ./run_job.sh path/to/script.R
LOGDIR="logs"
mkdir -p "$LOGDIR"                        # ensure logs folder exists

if [[ ! -f "$SCRIPT" ]]; then
  echo "ERROR: script not found: $SCRIPT" # guard: missing script
  exit 1
fi

ts="$(date +%F_%H%M%S)"                   # timestamp
base="$(basename "$SCRIPT" .R)"           # script name without .R
log="$LOGDIR/${base}_${ts}.log"           # log file path
pidfile="$LOGDIR/${base}.pid"             # pid file path

# Limit threaded numeric libraries inside each parallel R worker.
# Without these, 15 R workers can each spawn many BLAS/OpenMP threads.
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OMP_THREAD_LIMIT="${OMP_THREAD_LIMIT:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export MKL_DOMAIN_NUM_THREADS="${MKL_DOMAIN_NUM_THREADS:-1}"
export BLIS_NUM_THREADS="${BLIS_NUM_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-1}"
export RCPP_PARALLEL_NUM_THREADS="${RCPP_PARALLEL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"
export GOTO_NUM_THREADS="${GOTO_NUM_THREADS:-1}"
export OMP_DYNAMIC="${OMP_DYNAMIC:-FALSE}"
export MKL_DYNAMIC="${MKL_DYNAMIC:-FALSE}"

# Run R in a NEW SESSION so Ctrl+C in this terminal won't kill the job
setsid bash -lc "Rscript '$SCRIPT' > '$log' 2>&1" </dev/null &  # detached job
pid=$!                                    # PID of the detached session leader
echo "$pid" > "$pidfile"                  # save PID for status/kill

echo "Started: $SCRIPT"                   # info
echo "PID: $pid  (saved in $pidfile)"     # info
echo "Log: $log"                          # info
echo "Stop: kill -- -$pid"                # kill process group: master + workers

# make it runnable typing: chmod +x run_job.sh

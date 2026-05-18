#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DEFAULT="Rscripts/KELM_Ale_new.R"   # default script if none given
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

# Run R in a NEW SESSION so Ctrl+C in this terminal won't kill the job
setsid bash -lc "Rscript '$SCRIPT' > '$log' 2>&1" </dev/null &  # detached job
pid=$!                                    # PID of the detached session leader
echo "$pid" > "$pidfile"                  # save PID for status/kill

echo "Started: $SCRIPT"                   # info
echo "PID: $pid  (saved in $pidfile)"     # info
echo "Log: $log"                          # info

# make it runnable typing: chmod +x run_job.sh
#!/bin/bash
# Wait for the Phase 5 sweep to exit, then run the mechanical close-out.
#
# Launched detached (setsid) so it survives the Claude session that started it.
# Bracket class in the pgrep pattern (trap 1): a bare "phase5_run.R" matches
# this script's own command line and would make the wait loop exit immediately.

cd /mnt/c/Rworking/negative-sgr || exit 1
LOG=analysis/overnight_finish.log

{
  echo "=== waiting for the sweep to exit; started $(date) ==="
  # Bounded wait: 6 h is well past the ~5 h a cell takes, so if this expires
  # something is wrong and the close-out should say so rather than hang silently.
  deadline=$(( $(date +%s) + 21600 ))
  while pgrep -f "phase5_ru[n].R" > /dev/null; do
    if [ "$(date +%s)" -gt "$deadline" ]; then
      echo "TIMED OUT after 6 h with the sweep still running. Close-out NOT run."
      echo "The sweep is untouched and still going; nothing has been damaged."
      exit 1
    fi
    sleep 60
  done
  echo "=== sweep process gone at $(date); starting close-out ==="
  Rscript analysis/overnight_finish.R
  echo "=== close-out exit status $? at $(date) ==="
} >> "$LOG" 2>&1

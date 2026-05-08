#!/usr/bin/env bash
# todo-lock.sh <dir> acquire|release
# Simple file-based lock using a .lock file in <dir>.
# Uses flock if available, falls back to mkdir-based locking on Windows/Git Bash.
set -euo pipefail

DIR="${1:-.}"
CMD="${2:-acquire}"
LOCKFILE="${DIR}/.todo.lock"
MAX_WAIT=30  # seconds

case "$CMD" in
  acquire)
    # Try flock first (Linux/macOS/Git Bash with util-linux)
    if command -v flock &>/dev/null; then
      exec 9>"$LOCKFILE"
      flock -w "$MAX_WAIT" 9 || { echo "dispatch: failed to acquire lock on $LOCKFILE" >&2; exit 1; }
      echo "$$" >> "$LOCKFILE"
    else
      # mkdir-based fallback (atomic on all POSIX filesystems)
      LOCKDIR="${LOCKFILE}.d"
      waited=0
      while ! mkdir "$LOCKDIR" 2>/dev/null; do
        sleep 0.2
        waited=$(echo "$waited + 0.2" | bc 2>/dev/null || echo $((waited + 1)))
        if (( ${waited%.*} >= MAX_WAIT )); then
          echo "dispatch: timeout waiting for lock $LOCKDIR" >&2
          exit 1
        fi
      done
      echo "$$" > "$LOCKDIR/pid"
    fi
    ;;
  release)
    if command -v flock &>/dev/null; then
      flock -u 9 2>/dev/null || true
      rm -f "$LOCKFILE"
    else
      LOCKDIR="${LOCKFILE}.d"
      rm -rf "$LOCKDIR"
    fi
    ;;
  *)
    echo "Usage: todo-lock.sh <dir> acquire|release" >&2
    exit 1
    ;;
esac

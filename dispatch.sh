#!/usr/bin/env bash
# dispatch.sh — Manual dispatcher for TreeSearch Claude Code agents
# Subcommands: allocate <budget>  task <T-ID> [budget]  list
#              checkin <agent-id> ...  reap  kill <agent-id>
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STATE_FILE="$REPO_ROOT/.dispatch/state.json"
TODO_FILE="$REPO_ROOT/to-do.md"
LOG_DIR="$REPO_ROOT/.dispatch/logs"
RANKER_PROMPT="$REPO_ROOT/dev/dispatch/ranker.txt"
AGENT_BRIEF="$REPO_ROOT/dev/dispatch/agent-brief.md"

err() { echo "dispatch: $*" >&2; }
die() { err "$*"; exit 1; }

# ---------------------------------------------------------------------------
# Budget parser: accepts 5%/5h  2%/wk  15m  30min  0.5h
# Sets globals BUDGET_MINUTES, BUDGET_WINDOW, BUDGET_PCT
# ---------------------------------------------------------------------------
parse_budget() {
  local raw="$1"
  BUDGET_PCT=0; BUDGET_WINDOW="fixed"; BUDGET_MINUTES=0

  if [[ "$raw" =~ ^([0-9]+(\.[0-9]*)?)%/(5h|wk|weekly)$ ]]; then
    BUDGET_PCT="${BASH_REMATCH[1]}"
    local window="${BASH_REMATCH[3]}"
    [[ "$window" == "5h" ]] && BUDGET_WINDOW="5h" || BUDGET_WINDOW="weekly"
    local total_min; [[ "$BUDGET_WINDOW" == "5h" ]] && total_min=300 || total_min=2520
    BUDGET_MINUTES=$(echo "$BUDGET_PCT * $total_min / 100" | bc 2>/dev/null || echo "0")
    BUDGET_MINUTES=${BUDGET_MINUTES%.*}
    [[ -z "$BUDGET_MINUTES" || "$BUDGET_MINUTES" -eq 0 ]] && BUDGET_MINUTES=1
    return 0
  fi
  if [[ "$raw" =~ ^([0-9]+)(m|min)$ ]]; then
    BUDGET_MINUTES="${BASH_REMATCH[1]}"; return 0
  fi
  if [[ "$raw" =~ ^([0-9]+(\.[0-9]*)?)h$ ]]; then
    BUDGET_MINUTES=$(echo "${BASH_REMATCH[1]} * 60" | bc 2>/dev/null || echo "0")
    BUDGET_MINUTES=${BUDGET_MINUTES%.*}
    [[ -z "$BUDGET_MINUTES" || "$BUDGET_MINUTES" -eq 0 ]] && BUDGET_MINUTES=1
    return 0
  fi
  die "unrecognised budget format '$raw'; use e.g. 5%/5h  2%/wk  15m  30min  0.5h"
}

# ---------------------------------------------------------------------------
# Lock helpers
# ---------------------------------------------------------------------------
acquire_lock() { bash "$REPO_ROOT/todo-lock.sh" "$REPO_ROOT" acquire; }
release_lock() { bash "$REPO_ROOT/todo-lock.sh" "$REPO_ROOT" release; }

# ---------------------------------------------------------------------------
# State helpers
# ---------------------------------------------------------------------------
read_state() {
  [[ -f "$STATE_FILE" ]] && cat "$STATE_FILE" \
    || echo '{"agents":[],"budget":{"5h_window_started":null,"5h_pct_committed":0,"weekly_pct_committed":0}}'
}
write_state() { echo "$1" > "$STATE_FILE"; }

# Return lowest dN id not in state.json
mint_id() {
  local state="$1" n=1
  while true; do
    local id="d${n}"
    echo "$state" | jq -e --arg id "$id" '.agents[] | select(.id == $id)' &>/dev/null \
      || { echo "$id"; return; }
    (( n++ ))
  done
}

# ---------------------------------------------------------------------------
# Parse per-task hints: [m:haiku e:low] etc.  Sets MODEL_HINT, EFFORT_HINT
# ---------------------------------------------------------------------------
parse_hints() {
  local notes="$1"
  MODEL_HINT=""; EFFORT_HINT=""
  [[ "$notes" =~ \[([^]]*)\] ]] || return 0
  local hint="${BASH_REMATCH[1]}"
  if [[ "$hint" =~ m:([a-z0-9_-]+) ]]; then
    case "${BASH_REMATCH[1]}" in
      haiku*)  MODEL_HINT="claude-haiku-4-5"  ;;
      sonnet*) MODEL_HINT="claude-sonnet-4-6" ;;
      opus*)   MODEL_HINT="claude-opus-4-7"   ;;
    esac
  fi
  [[ "$hint" =~ e:(low|medium|high) ]] && EFFORT_HINT="${BASH_REMATCH[1]}"
}

# ---------------------------------------------------------------------------
# Heuristic model/effort from task text
# ---------------------------------------------------------------------------
heuristic_model() {
  local text="${1,,}"
  if [[ "$text" =~ (housekeep|doc|triage|spelling|readme|comment|typo|format) ]]; then
    SELECTED_MODEL="claude-haiku-4-5"; SELECTED_EFFORT="low"; return
  fi
  if [[ "$text" =~ (architect|red.?team|refactor|design|overhaul|rewrite) ]]; then
    SELECTED_MODEL="claude-opus-4-7"; SELECTED_EFFORT="high"; return
  fi
  SELECTED_MODEL="claude-sonnet-4-6"; SELECTED_EFFORT="medium"
}

# ---------------------------------------------------------------------------
# Find first OPEN task not in-flight; echo "ID|pri|status|desc|notes"
# ---------------------------------------------------------------------------
find_open_task() {
  local inflight_ids="$1"
  while IFS='|' read -r id pri status blocks desc notes; do
    id="${id// /}"; status="${status// /}"
    [[ "$id" =~ ^[TS]-[0-9]+ ]] || continue
    [[ "$status" == "OPEN" ]] || continue
    local skip=0
    for fid in $inflight_ids; do [[ "$id" == "$fid" ]] && skip=1 && break; done
    [[ "$skip" -eq 1 ]] && continue
    echo "${id}|${pri// /}|${status}|${desc}|${notes}"; return 0
  done < <(grep '^\s*|' "$TODO_FILE")
  return 1
}

# ---------------------------------------------------------------------------
# Haiku ranker call; returns JSON or empty on failure
# ---------------------------------------------------------------------------
call_ranker() {
  local excerpt="$1" inflight="$2" budget_min="$3"
  [[ ! -f "$RANKER_PROMPT" ]] && return 1
  local prompt; prompt="$(cat "$RANKER_PROMPT")"$'\n\n'"TO-DO EXCERPT:\n${excerpt}\n\nIN-FLIGHT: ${inflight}\nBUDGET: ${budget_min} minutes"
  claude -p --model claude-haiku-4-5-20251001 --max-turns 1 "$prompt" 2>/dev/null || return 1
}

# ---------------------------------------------------------------------------
# Edit to-do.md: set the Status column (col 3) for a given task ID
# action: set_status  new_status  (OPEN / "ASSIGNED (d1)" / "PARKED (...)")
#         or: reset_open  (used by kill: ASSIGNED|PARKED → OPEN)
# Uses python3 (available on this Windows machine via Python 3.12)
# ---------------------------------------------------------------------------
_TODO_EDITOR='
import sys
fname, tid, new_st, *rest = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
reset = rest and rest[0] == "reset"
lines = open(fname, encoding="utf-8").readlines()
out = []
for line in lines:
    cols = line.split("|")
    if len(cols) > 3 and cols[1].strip() == tid:
        cur = cols[3].strip()
        if reset and ("ASSIGNED" in cur or "PARKED" in cur):
            cols[3] = " OPEN "; line = "|".join(cols)
        elif not reset and cur == "OPEN":
            cols[3] = f" {new_st} "; line = "|".join(cols)
        elif not reset and "ASSIGNED" in cur:
            cols[3] = f" {new_st} "; line = "|".join(cols)
    out.append(line)
open(fname, "w", encoding="utf-8").writelines(out)
'

todo_set_status() {
  local task_id="$1" new_status="$2" reset="${3:-}"
  python3 -c "$_TODO_EDITOR" "$TODO_FILE" "$task_id" "$new_status" ${reset:+reset} 2>/dev/null || \
    err "warning: could not update to-do.md for $task_id"
}

# ---------------------------------------------------------------------------
# Spawn local agent (detached), return PID
# ---------------------------------------------------------------------------
spawn_local() {
  local agent_id="$1" task_id="$2" model="$3" effort="$4" budget_min="$5"
  local logfile="$LOG_DIR/${agent_id}-${task_id}.log"
  mkdir -p "$LOG_DIR"
  local brief=""; [[ -f "$AGENT_BRIEF" ]] && brief="$(cat "$AGENT_BRIEF")"$'\n\n'
  local prompt="${brief}Agent ID: ${agent_id}
Task: ${task_id}  Model: ${model}  Effort: ${effort}  Budget: ${budget_min} minutes
Repo: ${REPO_ROOT}
When pausing: bash dispatch.sh checkin ${agent_id} --kind=<gha|hamilton|human|other> --ref=<id> --eta=<iso> --resume=\"<action>\"
When complete: bash dispatch.sh checkin ${agent_id} --done"
  nohup claude -p --model "$model" "$prompt" > "$logfile" 2>&1 &
  echo $!
}

# ---------------------------------------------------------------------------
# Spawn remote (stub)
# ---------------------------------------------------------------------------
spawn_remote() {
  local agent_id="$1" task_id="$2" model="$3" effort="$4" budget_min="$5"
  local logfile="$LOG_DIR/${agent_id}-${task_id}.log"
  mkdir -p "$LOG_DIR"
  { echo "TODO: wire RemoteTrigger / scheduled-agent"
    echo "Intended: agent=${agent_id} task=${task_id} model=${model} effort=${effort} budget=${budget_min}m"
  } | tee "$logfile" >&2
  echo "0"
}

# ===========================================================================
# SUBCOMMANDS
# ===========================================================================

cmd_list() {
  read_state | jq '{agents:[.agents[]|{id,task_id,model,effort,spawn_kind,spawned_at,next_checkin}]}'
}

cmd_allocate() {
  local budget_raw="${1:-}"
  [[ -z "$budget_raw" ]] && die "usage: dispatch.sh allocate <budget>"
  parse_budget "$budget_raw"

  acquire_lock; trap 'release_lock' EXIT
  local state; state=$(read_state)

  # Build in-flight set from state + to-do.md ASSIGNED/WORKTREE rows
  local inflight_ids; inflight_ids=$(echo "$state" | jq -r '.agents[].task_id' | tr '\n' ' ')
  while IFS='|' read -r id _ status _; do
    id="${id// /}"; status="${status// /}"
    [[ "$id" =~ ^[TS]-[0-9]+ ]] && [[ "$status" =~ ^(ASSIGNED|WORKTREE) ]] && inflight_ids="$inflight_ids $id"
  done < <(grep '^\s*|' "$TODO_FILE")

  local task_id="" task_desc="" task_notes=""
  local ranker_out; ranker_out=$(call_ranker "$(grep '^\s*|' "$TODO_FILE" | head -80)" \
    "$inflight_ids" "$BUDGET_MINUTES" 2>/dev/null) || ranker_out=""

  if [[ -n "$ranker_out" ]]; then
    task_id=$(echo "$ranker_out" | jq -r '.task_id // empty' 2>/dev/null) || task_id=""
    SELECTED_MODEL=$(echo "$ranker_out" | jq -r '.model // empty' 2>/dev/null) || SELECTED_MODEL=""
    SELECTED_EFFORT=$(echo "$ranker_out" | jq -r '.effort // empty' 2>/dev/null) || SELECTED_EFFORT=""
  fi

  if [[ -z "$task_id" ]]; then
    local task_row; task_row=$(find_open_task "$inflight_ids") || die "no OPEN tasks available"
    task_id=$(echo "$task_row" | cut -d'|' -f1)
    task_desc=$(echo "$task_row" | cut -d'|' -f4)
    task_notes=$(echo "$task_row" | cut -d'|' -f5)
  else
    local task_row; task_row=$(grep -E "^\s*\| ${task_id} \|" "$TODO_FILE" | head -1) || true
    task_desc=$(echo "$task_row" | cut -d'|' -f5)
    task_notes=$(echo "$task_row" | cut -d'|' -f6)
  fi

  parse_hints "$task_notes"
  [[ -n "$MODEL_HINT" ]]  && SELECTED_MODEL="$MODEL_HINT"
  [[ -n "$EFFORT_HINT" ]] && SELECTED_EFFORT="$EFFORT_HINT"
  [[ -z "${SELECTED_MODEL:-}" ]] && heuristic_model "$task_desc $task_notes"

  _do_dispatch "$task_id" "$task_desc" "$SELECTED_MODEL" "$SELECTED_EFFORT" \
    "$BUDGET_MINUTES" "$BUDGET_WINDOW" "${BUDGET_PCT:-0}" "$state"
}

cmd_task() {
  local task_id="${1:-}" budget_raw="${2:-}"
  [[ -z "$task_id" ]] && die "usage: dispatch.sh task <T-ID> [budget]"
  BUDGET_MINUTES=0; BUDGET_WINDOW="fixed"; BUDGET_PCT=0
  [[ -n "$budget_raw" ]] && parse_budget "$budget_raw"

  acquire_lock; trap 'release_lock' EXIT
  local state; state=$(read_state)

  local inflight_ids; inflight_ids=$(echo "$state" | jq -r '.agents[].task_id' | tr '\n' ' ')
  for fid in $inflight_ids; do [[ "$fid" == "$task_id" ]] && die "task $task_id is already in-flight"; done

  local task_row; task_row=$(grep -E "^\s*\| ${task_id} \|" "$TODO_FILE" | head -1) \
    || die "task $task_id not found in to-do.md"
  local task_desc; task_desc=$(echo "$task_row" | cut -d'|' -f5)
  local task_notes; task_notes=$(echo "$task_row" | cut -d'|' -f6)

  parse_hints "$task_notes"
  [[ -n "$MODEL_HINT" ]]  && SELECTED_MODEL="$MODEL_HINT"
  [[ -n "$EFFORT_HINT" ]] && SELECTED_EFFORT="$EFFORT_HINT"
  [[ -z "${SELECTED_MODEL:-}" ]] && heuristic_model "$task_desc $task_notes"

  _do_dispatch "$task_id" "$task_desc" "$SELECTED_MODEL" "$SELECTED_EFFORT" \
    "$BUDGET_MINUTES" "$BUDGET_WINDOW" "${BUDGET_PCT:-0}" "$state"
}

# Shared dispatch core (called under lock)
_do_dispatch() {
  local task_id="$1" task_desc="$2" model="$3" effort="$4"
  local budget_min="$5" budget_window="$6" budget_pct="$7" state="$8"

  local agent_id; agent_id=$(mint_id "$state")
  local now; now=$(date -u +%Y-%m-%dT%H:%M:%SZ)
  err "dispatching $agent_id → $task_id ($model / $effort / ${budget_min}m)"

  local pid=0
  pid=$(spawn_local "$agent_id" "$task_id" "$model" "$effort" "$budget_min")

  local new_entry
  new_entry=$(jq -n \
    --arg id "$agent_id" --arg tid "$task_id" --arg title "$task_desc" \
    --arg model "$model" --arg effort "$effort" --arg kind "local" \
    --argjson pid "$pid" --arg now "$now" \
    --argjson bmin "${budget_min:-0}" --arg bwin "$budget_window" \
    --argjson bpct "${budget_pct:-0}" \
    --arg log ".dispatch/logs/${agent_id}-${task_id}.log" \
    '{id:$id,task_id:$tid,task_title:$title,model:$model,effort:$effort,
      spawn_kind:$kind,pid:$pid,remote_id:null,spawned_at:$now,
      budget_slice:{window:$bwin,pct:$bpct,minutes_est:$bmin},
      next_checkin:null,log:$log,branch:null}')

  local new_state
  new_state=$(echo "$state" | jq \
    --argjson entry "$new_entry" --argjson bmin "${budget_min:-0}" --arg now "$now" \
    '.agents += [$entry]
     | if .budget["5h_window_started"]==null then .budget["5h_window_started"]=$now else . end
     | .budget["5h_pct_committed"] += ($bmin*100/300|floor)
     | .budget["weekly_pct_committed"] += ($bmin*100/2520|floor)')
  write_state "$new_state"

  todo_set_status "$task_id" "ASSIGNED (${agent_id})"

  err "spawned $agent_id (PID $pid) → log: .dispatch/logs/${agent_id}-${task_id}.log"
  echo "{\"agent_id\":\"${agent_id}\",\"task_id\":\"${task_id}\",\"model\":\"${model}\",\"effort\":\"${effort}\",\"pid\":${pid}}"
}

cmd_checkin() {
  local agent_id="${1:-}"; [[ -z "$agent_id" ]] && die "usage: dispatch.sh checkin <agent-id> [flags]"
  shift
  local kind="" ref="" eta="" resume="" done_flag=0
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --kind=*)   kind="${1#--kind=}"   ;;
      --ref=*)    ref="${1#--ref=}"     ;;
      --eta=*)    eta="${1#--eta=}"     ;;
      --resume=*) resume="${1#--resume=}" ;;
      --done)     done_flag=1           ;;
      *) die "unknown checkin flag: $1" ;;
    esac; shift
  done

  acquire_lock; trap 'release_lock' EXIT
  local state; state=$(read_state)

  local task_id
  task_id=$(echo "$state" | jq -r --arg id "$agent_id" '.agents[]|select(.id==$id)|.task_id')
  [[ -z "$task_id" ]] && die "agent $agent_id not found in state.json"

  if [[ "$done_flag" -eq 1 ]]; then
    write_state "$(echo "$state" | jq --arg id "$agent_id" '.agents=[.agents[]|select(.id!=$id)]')"
    err "agent $agent_id marked done"
    echo "{\"agent_id\":\"${agent_id}\",\"done\":true}"
  else
    local ci; ci=$(jq -n --arg k "$kind" --arg r "$ref" --arg e "$eta" --arg a "$resume" \
      '{kind:$k,ref:$r,eta:$e,resume_action:$a}')
    write_state "$(echo "$state" | jq --arg id "$agent_id" --argjson ci "$ci" \
      '.agents=[.agents[]|if .id==$id then .next_checkin=$ci else . end]')"
    todo_set_status "$task_id" "PARKED (${agent_id}, ${kind} ${ref})"
    err "agent $agent_id parked; ETA $eta"
    echo "{\"agent_id\":\"${agent_id}\",\"parked\":true,\"eta\":\"${eta}\"}"
  fi
}

cmd_reap() {
  local state; state=$(read_state)
  local now; now=$(date -u +%Y-%m-%dT%H:%M:%SZ)
  local ready_ids
  ready_ids=$(echo "$state" | jq -r --arg now "$now" \
    '.agents[]|select(.next_checkin!=null and .next_checkin.eta!=null and .next_checkin.eta!="" and .next_checkin.eta<=$now)|.id')
  if [[ -z "$ready_ids" ]]; then echo '{"ready":[]}'; return 0; fi

  acquire_lock; trap 'release_lock' EXIT
  state=$(read_state)
  local new_state
  new_state=$(echo "$state" | jq --arg now "$now" \
    '.agents=[.agents[]|if (.next_checkin!=null and .next_checkin.eta!=null and .next_checkin.eta!="" and .next_checkin.eta<=$now) then .+{ready:true} else . end]')
  write_state "$new_state"
  echo "$new_state" | jq '{ready:[.agents[]|select(.ready==true)|{id,task_id,next_checkin}]}'
}

cmd_kill() {
  local agent_id="${1:-}"; [[ -z "$agent_id" ]] && die "usage: dispatch.sh kill <agent-id>"

  acquire_lock; trap 'release_lock' EXIT
  local state; state=$(read_state)

  local pid task_id
  pid=$(echo "$state" | jq -r --arg id "$agent_id" '.agents[]|select(.id==$id)|.pid // 0')
  task_id=$(echo "$state" | jq -r --arg id "$agent_id" '.agents[]|select(.id==$id)|.task_id')
  [[ -z "$task_id" ]] && die "agent $agent_id not found"

  if [[ -n "$pid" ]] && (( pid > 0 )) 2>/dev/null; then
    kill -TERM "$pid" 2>/dev/null && err "sent SIGTERM to PID $pid" || err "PID $pid already gone"
  fi

  write_state "$(echo "$state" | jq --arg id "$agent_id" '.agents=[.agents[]|select(.id!=$id)]')"
  todo_set_status "$task_id" "OPEN" "reset"
  err "killed $agent_id (task $task_id reset to OPEN)"
  echo "{\"agent_id\":\"${agent_id}\",\"killed\":true,\"task_id\":\"${task_id}\"}"
}

# ===========================================================================
# Entry point
# ===========================================================================
SUBCMD="${1:-}"; shift || true
case "$SUBCMD" in
  allocate) cmd_allocate "$@" ;;
  task)     cmd_task     "$@" ;;
  list)     cmd_list          ;;
  checkin)  cmd_checkin  "$@" ;;
  reap)     cmd_reap          ;;
  kill)     cmd_kill     "$@" ;;
  "")       die "usage: dispatch.sh <allocate|task|list|checkin|reap|kill> ..." ;;
  *)        die "unknown subcommand: $SUBCMD" ;;
esac

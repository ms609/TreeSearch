#!/usr/bin/env bash
# Thin wrapper — engine lives in the shared dispatch skill.
exec bash "$HOME/.claude/skills/dispatch/todo-lock.sh" "$@"

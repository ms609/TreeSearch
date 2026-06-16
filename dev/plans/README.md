# dev/plans

Design and strategy documents for TreeSearch development, tracked in git.

This directory is the **Claude-convention home** for plans. It supersedes the
git-ignored `.positai/plans/` tree (the PositAI tool's plan store), which is
**retired** as of 2026-06-16. Historical PositAI plans remain available locally
under `.positai/plans/` but are no longer added to or maintained; their live
conclusions have been carried forward into the docs here and into the
file-based memory.

Conventions:
- One markdown file per plan, `YYYY-MM-DD-short-slug.md`.
- Plans are living specs: update status inline as work progresses.
- Durable cross-session facts (one fact each) go in the memory store, not here;
  link from a plan to memory by name where useful.

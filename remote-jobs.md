# Pending Remote Jobs

Track asynchronous jobs (Hamilton SLURM, long-running GHA, etc.) that
produce results an agent needs to retrieve later.

## How this works

- **Add a row** when you submit a remote job whose results won't be
  consumed in the same conversation turn.
- **Delete the row** once results have been retrieved and acted on
  (committed to repo, written up in coordination.md, etc.).
- Agents check this file at `/assign` time, after triaging `a.*` and
  `u.*` files but before claiming from `to-do.md`. If a retrievable
  job is listed, retrieving and processing it takes priority.
- Use the lock (`bash ../../todo-lock.sh . acquire/release`) if editing
  concurrently, same as `to-do.md`.

## Jobs

| Submitted | Type | Job ID | Branch | Description | Retrieve how | Owner |
|-----------|------|--------|--------|-------------|-------------|-------|

<!-- Example row:
| 2026-03-29 | SLURM | 16622483 | cpp-search | T-289f Stage 5: PR NNI polish benchmark (5 datasets, 131–206t) | `scp hamilton:scratch/ts_bench/t289f_*.csv dev/benchmarks/` | E |
-->

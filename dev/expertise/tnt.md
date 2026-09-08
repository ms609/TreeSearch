# TNT (Tree analysis using New Technology)

## Installation

TNT is installed at `C:\Programs\Phylogeny\tnt\`.

### Executables

| Path | Version | Notes |
|------|---------|-------|
| `tnt/tnt.exe` | older | **Do not use.** |
| `tnt/TNT-bin/tnt.exe` | 1.6 | **Use this one.** Console/script mode. |
| `tnt/TNT-bin/wTNT.exe` | 1.6 | Windows GUI version. |

Always use `C:\Programs\Phylogeny\tnt\TNT-bin\tnt.exe` (version 1.6).

### Invocation

**Never launch TNT without passing a script file.** TNT defaults to
interactive mode and will block waiting for keyboard input, hanging any
automated pipeline.

**Correct pattern** — pass a `.run` script as a positional argument with
trailing semicolon:

```bash
"C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe" "myscript.run;"
```

This launches TNT in PISH (batch) mode. It reads and executes the script,
then exits when it hits `quit;`.

**Critical: script files must use `.run` extension.** TNT interprets `.tnt`
files as data files. If you pass a `.tnt` script, TNT will try to parse it
as data and fail with "Can't open .tnt".

**Critical: script filenames must be purely alphabetic (no digits or
underscores).** TNT parses the filename as a command line — it splits on
digits and underscores, treating the first alphabetic token as a command.
`bench1.run` → command `bench`; `Vinther2008_EW.run` → command `vinther`.
Safe names: `tntbench.run`, `mytest.run`, `abc.run`.

**Piping via stdin does NOT work reliably** — `echo "..." | tnt.exe` launches
interactive mode (shows ASCII banner) and may hang.

**Encoding**: TNT stdout contains non-UTF8 progress bar characters. Use
`iconv(output, from = "", to = "UTF-8", sub = "")` to sanitize before
regex matching in R.

### TNT script basics

- Commands are terminated by `;`
- `mxram N;` — set memory (MB); must be first command
- `proc <file>;` — read data file (TNT `.tnt` or Nexus format)
- `xmult;` — heuristic search (new technology search)
- `xmult=hits N replic M;` — search with convergence/replicate limits
- `piwe = K;` — implied weights with concavity constant K
- `xpiwe = K;` — extended implied weights
- `rseed N;` — set random seed
- `timeout HH:MM:SS;` — set search time limit
- `best;` — report best score and tree count
- `length;` — print tree lengths
- `quit;` — exit TNT (essential for non-interactive use)

### Data format

TNT can read NEXUS (`.nex`) files and its own format (`.tnt`).
For NEXUS input, use `proc <file.nex>;`.

Export from R: `TreeTools::WriteTntCharacters(phyDat_obj, filepath)`.

### Output parsing

TNT stdout contains parseable lines:
- `"Best score: 78."` or `"Best score: 3.80000."` (IW) — best score
- `"N trees retained"` — number of trees found
- `"Best score hit N times."` — convergence hits
- `"Total rearrangements examined: N."` — total rearrangements

### Score comparability with TreeSearch

TNT standard Fitch treats inapplicable tokens as a regular character state
(column-based). TreeSearch uses Brazeau et al. (2019) three-pass algorithm.
For datasets with inapplicable characters, TNT EW scores will generally be
≤ TreeSearch EW scores. For IW, both use Goloboff's `e/(k+e)` formula.

Example: Vinther2008 — TNT EW = 78, TreeSearch EW = 79.

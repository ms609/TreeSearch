#!/bin/bash
# Diagnose why 64-bit TNT dies headless under SLURM/non-interactive ssh.
# Tests: plain (TERM=xterm) vs pseudo-TTY wrap (script -qec).
TNT=/nobackup/$USER/TreeSearch/tnt/TNT-bin/tnt
echo "=== TNT binary ==="; ls -la "$TNT" || { echo "NO TNT"; exit 1; }
W=/tmp/tnttest_$$; mkdir -p "$W"; cd "$W"
cat > d.tnt <<'DAT'
xread
2 4
A 01
B 01
C 10
D 10
;
DAT
cat > r.run <<'RUN'
log out.log;
proc d.tnt;
hold 100;
mult;
best;
log/;
quit;
RUN

echo "=== ATTEMPT 1: plain, TERM=xterm ==="
rm -f out.log
TERM=xterm "$TNT" "r.run;" >stdout1.txt 2>stderr1.txt
echo "exit=$?"
echo "--- out.log ---"; cat out.log 2>/dev/null || echo "(no log file)"
echo "--- stderr1 tail ---"; tail -6 stderr1.txt

echo "=== ATTEMPT 2: pseudo-TTY via script -qec, TERM=xterm ==="
rm -f out.log
TERM=xterm script -qec "$TNT r.run;" /dev/null >stdout2.txt 2>stderr2.txt
echo "exit=$?"
echo "--- out.log ---"; cat out.log 2>/dev/null || echo "(no log file)"
echo "--- stdout2 tail ---"; tail -10 stdout2.txt

echo "=== ATTEMPT 3: commands via stdin pipe, no run file, TERM=dumb ==="
rm -f out.log
printf 'log out.log;\nproc d.tnt;\nhold 100;\nmult;\nbest;\nlog/;\nquit;\n' | TERM=dumb "$TNT" >stdout3.txt 2>stderr3.txt
echo "exit=$?"
echo "--- out.log ---"; cat out.log 2>/dev/null || echo "(no log file)"
echo "--- stdout3 tail ---"; tail -10 stdout3.txt

rm -rf "$W"

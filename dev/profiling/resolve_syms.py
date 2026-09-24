import subprocess, re, bisect
dll = r"dev/profiling/.vtune-lib-resurvey/TreeSearch/libs/x64/TreeSearch.dll"
# nm output: "<hex addr> <type> <name>"
syms = []
out = subprocess.run(["nm", "-C", dll], capture_output=True, text=True).stdout
for ln in out.splitlines():
    m = re.match(r"^([0-9a-fA-F]{8,16})\s+\S\s+(.+)$", ln)
    if m:
        syms.append((int(m.group(1), 16), m.group(2)))
syms.sort()
addrs = [s[0] for s in syms]
# anchor: _pthread_self_lite runtime VA from VTune = 0x2cc2757c4
anchor_name = "_pthread_self_lite"
anchor_static = None
for a, n in syms:
    if anchor_name in n:
        anchor_static = a; break
slide = 0  # CONFIRMED: VTune addrs == nm static addrs (two exact matches)
print("slide =", hex(slide))
hot = {0x2cc1eeab0:1.653,0x2cc27ff80:0.761,0x2cc1e5980:0.594,0x2cc24b86f:0.232,
       0x2cc1e3d60:0.093,0x2cc249870:0.063,0x2cc1e47b0:0.047,0x2cc253280:0.031,
       0x2cc2444c0:0.031,0x2cc31db50:0.016,0x2cc1e5330:0.016}
print("--- hot TreeSearch.dll funcs (self s) -> symbol ---")
for rt, t in sorted(hot.items(), key=lambda x: -x[1]):
    st = rt - slide
    i = bisect.bisect_right(addrs, st) - 1
    name = syms[i][1] if i >= 0 else "?"
    off = st - syms[i][0] if i >= 0 else 0
    print(f"0x{rt:x}  self={t:5.3f}  +{off:#x}  {name}")
# also locate the structural fns of interest by name
print("--- key fns by name (RVA-ish static addr) ---")
for want in ["build_postorder","compute_insertion_edge_sets","tbr_search","fitch_indirect","any_hit_reduce","compute_from_above"]:
    hits = [(a,n) for a,n in syms if want in n]
    for a,n in hits[:2]:
        print(f"  static 0x{a:x}  runtime~0x{a+slide:x}  {n}")

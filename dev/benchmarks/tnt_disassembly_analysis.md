# TNT vs TreeSearch: Fitch Kernel Disassembly Comparison (T-250)

Date: 2026-03-26

## Scope limitation

**This analysis covers the native Windows TNT binary only.** The TNT
download page explicitly labels the Windows build as "[32 bits]". The
Mac, Linux, and Cygwin builds are compiled as 64-bit (Goloboff & Morales
2023). The 64-bit builds likely use wider registers and may include SIMD
or hardware `popcnt` — the "~4× throughput advantage" conclusion below
does **not** generalize to 64-bit TNT. Hamilton HPC benchmarks (T-249)
will run against the 64-bit Linux TNT and may show a very different
implementation-level gap.

## TNT Binary Profile (Windows, 32-bit)

- **File:** `C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe` (3.1 MB)
- **Format:** PE32 (32-bit i386), stripped (no symbols)
- **SIMD:** None. Zero xmm/ymm register references, zero popcnt instructions.
- **Code section:** `AUTO` — 2.4 MB, ~721K disassembly lines

## TreeSearch DLL Profile

- **File:** `.agent-E/TreeSearch/libs/x64/TreeSearch.dll` (1.8 MB)
- **Format:** PE32+ (64-bit x86-64), stripped
- **SIMD:** SSE2 (128-bit). 1281 integer SIMD ops (pand/por/pxor/pcmpeq),
  16472 xmm register references (includes scalar double FP). Zero ymm (no AVX2).
- **Popcount:** Software Hamming weight (0x5555.../0x3333.../0x0f0f... shift-mask
  pattern). No hardware `popcnt` instruction.

## Comparison Table

| Feature | TNT | TreeSearch |
|---------|-----|------------|
| Architecture | 32-bit i386 | 64-bit x86-64 |
| Word size | 32-bit | 64-bit |
| SIMD for Fitch | None | SSE2 (128-bit `pand`/`por`) |
| Popcount | 64KB lookup table (two 16-bit halves) | Software Hamming weight (shift+mask) |
| Hardware `popcnt` | No | No |
| AVX2 | No | No |
| Bits/inner-loop iteration | 32 | 128 (2 × uint64 via `movdqu`/`pand`/`por`) |

## TNT Fitch Kernel (0x420c04)

The main scoring loop at 0x420c04–0x420c77 is a single-pass design:

```
loop:
  dec counter; cmp -1; je exit         // iterate over character words
  mov 0x4(%esi),%eax                   // load left child state (32 bits)
  mov 0x4(%ebx),%ecx                   // load right child state
  add $0x4,%ebx; add $0x4,%esi        // advance pointers (stride 4 = 32 bits)
  not %eax; and %ecx,%eax             // ~left & right = "extra states in right"
  je skip                              // if zero, no extra states
  mov %ecx,%edx; xor %eax,%edx        // right XOR extra = intersection
  push %eax; call 0x5c9f30            // popcount(extra) via 64KB LUT
  mov %edx,(%ebx)                     // store intersection
  sub %eax,(%edx)                     // adjust score counter
skip:
  [symmetric check: ~result & left]
  jmp loop
```

The popcount function at 0x5c9f30 splits a 32-bit value into two 16-bit halves
and uses a 64KB lookup table at 0x718dbd:

```
mov 0x8(%ebp),%edx       // arg
mov %edx,%eax
and $0xffff,%eax         // low 16 bits
shr $0x10,%edx           // high 16 bits
mov 0x718dbd(%eax),%al   // table[low]
add 0x718dbd(%edx),%al   // + table[high]
movsbl %al,%eax
ret
```

**Key characteristics:**
- Processes one 32-bit word per iteration
- NOT+AND pattern (computes "extra states" directly) rather than AND then check-zero
- Includes a symmetric second check (right-to-left and left-to-right in the same loop body)
- Function call for popcount (not inlined)
- Branch per character word (`je skip`)

## TreeSearch Fitch Kernel

### Indirect scoring (TBR inner loop — 72% of wall time at 180 tips)

`any_hit_reduce3()` in `ts_simd.h` is the critical inner function:

```cpp
v128 acc = zero128();
for (; s + 2 <= n_states; s += 2) {
  v128 vc = loadu128(&clip[s]);       // 128-bit load (2 × uint64)
  v128 va = loadu128(&a[s]);
  v128 vb = loadu128(&b[s]);
  acc = or128(acc, and128(vc, or128(va, vb)));  // clip & (a | b)
}
```

Compiled to:
```
movdqu (%r8,%rax),%xmm0     // load 128 bits from clip
movdqu (%r9,%rax),%xmm2     // load 128 bits from a|b (pre-computed or inline)
add $0x10,%rax               // stride 16 = 128 bits
pand %xmm2,%xmm0            // 128-bit AND
por %xmm0,%xmm1             // OR accumulate
cmp %rax,%rdx
jne loop
```

**Key characteristics:**
- Processes 128 bits per iteration (2 × uint64)
- SSE2 `pand`/`por` for bit operations
- Branchless within the character loop (no per-word branching)
- `popcount64()` on the result mask (software Hamming weight)

### Downpass (`fitch_downpass_node`)

Two-pass design:
1. **Pass 1:** `any_hit_reduce()` — tight SSE2 `pand`+`por` loop to determine
   which characters have intersection (single 64-bit mask)
2. **Pass 2:** Broadcast mask + SSE2 select — no per-character branching

## Implications

### TNT's speed advantage is NOT implementation-level (Windows 32-bit)

On Windows, TreeSearch has a **~4× raw Fitch throughput advantage** (128-bit
SSE2 vs 32-bit scalar). Yet TNT converges 3–5× faster on the same datasets.
This means — at minimum on Windows:

1. **TNT's advantage is purely strategic** — fewer candidates evaluated,
   more effective heuristics, or both.
2. **T-246 (AVX2)** would double TreeSearch's throughput from 128→256 bits
   (and could add hardware `popcnt`). This is still worthwhile for absolute
   speed, but it won't close the strategic gap with TNT.
3. **T-251 (trajectory analysis) is the higher-priority investigation** —
   understanding *how many* candidates TNT evaluates per score improvement
   will reveal whether the gap is in candidate pruning, search ordering,
   or phase composition.

### Minor optimization opportunities

- **Hardware `popcnt`:** Neither program uses it. Adding `-mpopcnt` to
  TreeSearch's compile flags (or runtime dispatch) would replace the
  ~10-instruction software Hamming weight with a single `popcntq`. This
  affects step counting after each `any_hit_reduce`, not the inner loop
  itself, but could save ~5–10% of scoring time.
- **TNT's popcount is worse:** The 64KB LUT + function call overhead is
  significantly more expensive than TreeSearch's inlined shift-mask.
  This further confirms TNT's advantage is strategic.

### What to investigate next

The round 2 data shows TNT completing 50+ trees in 7–27s while TreeSearch
takes 45–110s for similar scores. If TreeSearch's per-candidate scoring is
faster, TNT must be evaluating far fewer candidates to achieve the same
result — either through better candidate pruning (e.g., more aggressive
clip skipping, smarter regraft ordering) or through phases that escape
local optima more efficiently (more effective ratchet/drift parameters).

T-249 (rerun comparison) and T-251 (trajectory analysis) should focus on
comparing **total candidates evaluated** and **score improvement per candidate**
rather than wall-clock timing.

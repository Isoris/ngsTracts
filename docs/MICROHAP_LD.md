# Microhap block coherence (STEP_TRC_07 / MicrohapLD)

Population-level **block coherence** layer. Independent of, and conceptually
upstream of, the pedigree tract-calling layer: it asks whether candidate LRRs
are held together as coherent long-range haplotype blocks across the cohort,
using local microhaplotypes as the unit of association.

## Submodule architecture (where this sits)

```
ngsTracts
├── MicrohapBuilder        (submodule 1) — phased variants -> microhap objects   [EXTERNAL]
├── MicrohapLD / Coherence (submodule 2) — THIS step (STEP_TRC_07)                [in-repo]
├── GHSL                   (submodule 3) — per-sample hap1-vs-hap2 divergence     [EXTERNAL]
├── TractCaller            — inherited-state tracts (STEP_TRC_01/02)              [in-repo]
├── EventCaller            — CO/DCO/NCO-like (STEP_TRC_01/04)                     [in-repo]
└── RegionClassifier       — overlay GHSL+LD+FIS+events+masks                     [join]
```

## The unit: a microhap is a haplotype object

Instead of "is SNP A in LD with SNP B?", the question is "is microhap object A
associated with microhap object B?". Each window's microhap allele is a local
phased allele string across 7-30 variants. Within one arrangement (U or V), the
arrangement-consistent microhaps (`M1_u, M2_u, M3_u ...`) tend to travel
together; recombination mixes them (`M1_u, M2_v, M3_u`), weakening association.

## The statistic

For each window pair, microhap LD = **bias-corrected Cramér's V** (Bergsma) of
the population's chromosome-resolved microhap alleles — a multi-allelic
association in [0,1] (`scripts/STEP_TRC_07_microhap_ld.py:60`). Chromosome
resolution uses the phased `hap1`/`hap2`: `hap1@A` pairs with `hap1@B`
(`scripts/STEP_TRC_07_microhap_ld.py:96`).

> **Phase assumption:** hap1/hap2 must be consistently phased across windows
> (same phase set), or this LD is not interpretable. Confirm against the phaser.

## Per-LRR test and background comparison

For each LRR: take its microhap windows, keep long-distance pairs
(`--min-distance-bp`, default 100 kb), and summarize `median_microhap_LD`,
`p90`, `frac_pairs_high_LD` (V ≥ `--high-ld-threshold`, default 0.5), and an
`LD_decay_slope_per_Mb`. Compare to a **genome background** of long-range pairs
whose windows lie outside all LRRs (`scripts/STEP_TRC_07_microhap_ld.py:163`).

`coherence_rank` (`scripts/STEP_TRC_07_microhap_ld.py:200`):

| rank | rule |
| --- | --- |
| `strong` | median LD ≥ background p90 **and** ≥ high-LD threshold |
| `moderate` | above background p90 **or** ≥ threshold (not both) |
| `weak` | neither (background-like or rapidly decaying) |
| `low_confidence` | fewer than `--min-pairs` long-range pairs |

Elevation is judged against the background **upper tail (p90)**, not its median
(which is ~0 for unlinked windows and would mark anything positive as elevated).

## Output `lrr_microhap_coherence.tsv`

Ranked, one row per LRR: `n_windows`, `n_longrange_pairs`,
`median_microhap_LD`, `p90_microhap_LD`, `frac_pairs_high_LD`,
`LD_decay_slope_per_Mb`, `background_median_LD`, `background_p90_LD`,
`delta_vs_background`, `coherence_rank`. Per-pair detail with `--emit-pairs`.

## Overlay (downstream, external join)

- high GHSL + high microhap-LD → divergent haplotypes preserved as a block
- high GHSL + low microhap-LD → divergent but not coherent
- high microhap-LD + negative FIS → heterozygote-excess block
- high microhap-LD + positive FIS → heterozygote-deficit block

GHSL and FIS are popstats / unified-ancestry outputs; STEP_TRC_07 supplies the
coherence side. MicrohapBuilder (phased Clair3 → microhap objects) is the
external input producer; its output contract is the `--microhaps` columns
(`sample, chrom, start, end, hap1, hap2`).

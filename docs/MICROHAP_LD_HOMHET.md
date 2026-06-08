# Microhap HOM/HET block coherence (STEP_TRC_08 / MicrohapLD_HOMHET)

A focused extension of the MicrohapLD coherence test (STEP_TRC_07): within each
candidate LRR, split samples into arrangement-**HOM** and arrangement-**HET**,
and compute population microhap-microhap LD *separately* in each group. This
asks whether the block's haplotype structure is held together equally in
homozygotes and heterozygotes.

## Exact test

1. Build microhap objects inside the LRR (windows of phased variants).
2. Split samples into HOM and HET for that LRR (karyotype, from upstream).
3. Project each microhap against every other microhap in the LRR.
4. Compute population microhap-LD (bias-corrected Cramér's V over
   chromosome-resolved hap1/hap2) separately in HOM and in HET samples
   (`scripts/STEP_TRC_08_microhap_hom_het.py:64`,
   `scripts/STEP_TRC_08_microhap_hom_het.py:90`).
5. Compare median LD_HOM vs median LD_HET.

## Interpretation (block level)

| HOM | HET | label | meaning |
| --- | --- | --- | --- |
| high | high | `coherent_in_both` | block coherent in both genotype classes |
| high | reduced (≥ low, < high) | `coherent_HOM_reduced_HET` | homozygotes coherent, heterozygotes partly mixed |
| high | very low (< low) | `coherent_HOM_disrupted_HET` | heterozygotes disrupted: recombination / NCO-like local disruption / bad boundary / noise |
| low | high | `coherent_HET_only` | unusual — check calls |
| low | low | `weak_both` | not a coherent block |

Thresholds: `--high-ld-threshold` (default 0.5), `--low-ld-threshold`
(default 0.2). HOM groups below `--min-group-samples` (default 4) or with no
computable pair are `low_confidence_HOM`; missing HET data → `*_no_HET_data`
(`scripts/STEP_TRC_08_microhap_hom_het.py:170`).

## Outputs

- `lrr_microhap_hom_het.tsv` (block level, ranked by Δ = HOM − HET):
  `n_hom_samples`, `n_het_samples`, `n_pairs_hom`, `n_pairs_het`,
  `median_microhapLD_HOM`, `median_microhapLD_HET`, `delta_HOM_minus_HET`,
  `interpretation`.
- `microhap_hom_het_profile.tsv` (per microhap): `mean_LD_to_others_HOM`,
  `mean_LD_to_others_HET` — the projection of each microhap onto all others.

## Provenance

- **Upstream (external):** MicrohapBuilder (phased Clair3 → microhap objects)
  and the per-sample × LRR karyotype call (HOM_REF/HET/HOM_INV). Consumed via
  `--microhaps` and `--karyotypes`.
- **In-repo:** this HOM/HET LD comparison.
- Same cross-window phase-consistency assumption as STEP_TRC_07.

> For each candidate LRR, local phased variants were collapsed into microhap
> objects. Each microhap was tested for association with all other microhaps in
> the same region, separately among arrangement-homozygous and
> arrangement-heterozygous samples. The resulting HOM and HET microhapLD
> profiles summarize whether microhap states remain coherently associated across
> the region in each genotype class.

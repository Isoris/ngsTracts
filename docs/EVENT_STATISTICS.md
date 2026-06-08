# Event-statistics catalog (STEP_TRC_04)

Reporting template borrowed from the pedigree-recombination literature — the
**structure**, not the biological values. Each statistic is emitted by
`scripts/STEP_TRC_04_event_statistics.py` into `event_statistics.tsv` with a
`status` of `computed`, `model:<...>`, or `requires:<input>`. `requires:` rows
carry NA, never a fabricated number.

## Statistic table

| statistic | definition | status / required input | inspired_by |
| --- | --- | --- | --- |
| `n_CO` / `n_NCO` / `n_DCO` | event counts (one-way switch / short bracket / long bracket) | computed | Smeds 2016; Korunes 2019; Navarro 1997 |
| `n_*` by `inside_inversion` | counts split inside/outside/partial | computed | Korunes 2019 |
| `NCO_CO_ratio` | gene-conversion-to-crossover leak | computed | ngsTracts (gene-flux) |
| `DCO_CO_ratio` | even-exchange recovery diagnostic | computed | inversion heterokaryotype |
| `gene_flux_index_count` / `_bp` | NCO + DCO, count and span-weighted | computed | Navarro 1997 (gene flux = GC + DCO) |
| `CO_resolution_median_bp` | median flanking-marker interval (or refined CI width) | computed (sharper with `--traversal-breakpoints`) | Smeds 2016; deCODE 2019 |
| `CO_resolution_lt{T}bp_pct` | % COs localized below T (default 10 kb) | computed | Smeds 2016 (86% < 10 kb) |
| `NCO_length_median_bp` (ALL/inside/outside) | median observed tract length | computed | Korunes 2019; Betrán 1997 |
| `NCO_length_mean_bp`, `NCO_length_geometric_p_hat` | tract-length model (p̂ = 1/mean) | computed | Betrán 1997; Halldorsson 2016 |
| `event_distance_to_breakpoint_median_bp` | per-class distance to nearest inversion breakpoint | computed (from Stage-3 `distance_to_nearest_inv_bp`) | Korunes 2019; Koury 2023 |
| `event_distance_to_telomere_median_bp` | per-class distance to chromosome end | computed (needs `.fai`) | Smeds 2016 |
| `{NCO,CO}_rate_*_inside_outside` | per-Mb rates inside vs outside, and ratio | `requires:inversion-atlas` (for inside-bp denominator) | Korunes 2019 |
| `NCO_preserved_CO_suppressed_flag` | true-inversion vs cold-region discriminator | `requires:inversion-atlas` | inversion-vs-cold-region |
| `{CO,NCO}_count_inside_{Het,Hom}` | event counts by transmitting-parent karyotype | `requires:parent-karyotypes + inversion-atlas` | classical inversion theory; Korunes 2019 |
| `event_rate_by_parent_sex`, `event_rate_male_over_female` | sex-specific event rate and ratio | `requires:sample-sex` | Smeds 2016; Halldorsson 2016 |
| `CO_interference_inter_event_gap_median_bp`, `CO_inter_event_gap_CV` | spacing between COs (CV<1 ⇒ positive interference) | computed (needs ≥2 COs per dyad×chrom) | Smeds 2016 (interference ≤ 14 Mb) |
| `NCO_CO_spatial_assoc_nearest_median_bp` | NCO-to-nearest-CO distance | computed | Smeds 2016 |
| `event_density_<name>_frac_overlap` | event overlap with a feature track | `requires:annotation(NAME=BED)` | Smeds 2016 (promoters/CpG) |
| `event_distance_to_centromere_median_bp` | per-class distance to centromere/cold region | `requires:centromeres-bed` | cold-region control |
| `switch_error_rate`, `expected_false_NCO_upper_bound` | low-coverage safeguard (E ≈ n_sites·ε²) | `model:` (needs `--switch-error-rate`) | low-coverage safeguard |
| `GC_bias_at_NCO` | gBGC transmission bias toward G/C | `requires:external(popstats allele/ancestral)` | Smeds 2016; Halldorsson 2016 |
| `rho_inside_outside_ratio` | LD-based historical recombination proxy | `requires:external(LD map: pyrho/LDhat)` | LDhat; pyrho |
| `IBD_NCO_support_frac` | IBD-supported NCO fraction | `requires:external(IBD)` | Browning & Browning 2024 |

## Provenance discipline

- **In-repo / computed:** everything derivable from `tract_classifications.tsv`
  (+ `.fai`, + optional `traversal_breakpoints.tsv`).
- **Optional in-repo (wired, degrade gracefully):** inside/outside rates
  (`--inversion-atlas`), Het/Hom split (`--parent-karyotypes`), sex effect
  (`--sample-sex`), feature density (`--annotation`), centromere distance
  (`--centromeres`), false-NCO bound (`--switch-error-rate`).
- **External (popstats / unified-ancestry, LD, IBD engines):** gBGC, rho, IBD —
  emitted as `requires:` so the report is the full intended atlas table with the
  unverified cells clearly marked.

`event_statistics.provenance.tsv` records which optional inputs were supplied
for a given run.

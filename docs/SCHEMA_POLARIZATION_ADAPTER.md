# Polarization adapter contract (ngsPedigree -> ngsTracts)

This fixes the boundary between **arrangement polarization** (upstream, in
ngsPedigree) and **transmission / drive testing + recombination tract calling**
(downstream, in ngsTracts). ngsTracts **plugs in at the end of polarization** and
does not re-implement it.

## Who owns what

| Layer | Computation | Owner | In this repo? |
| --- | --- | --- | --- |
| 1 | Karyotype calling: sample x inversion -> HOM_REF / HET / HOM_INV (from the arrangement-specific marker axis / local PCA band) | ngsPedigree | no |
| 2 | Dyad/triad compatibility: detect impossible parent-offspring genotype combinations (HOM_REF parent cannot yield HOM_INV offspring, etc.) | ngsPedigree | no |
| 3 | Transmission polarization: pick the REF/INV orientation that minimizes Mendelian contradictions; call which arrangement each HET parent transmitted | ngsPedigree | no |
| **adapter** | **read polarized transmissions** | **ngsTracts STEP_TRC_05** | **yes** |
| 4a | Mendelian transmission / meiotic-drive binomial test | ngsTracts STEP_TRC_05 | yes |
| 4b | Recombination tract calling (CO/NCO/DCO from marker-level transmitted state) | ngsTracts STEP_TRC_01/02 | yes |

Layers 1-3 are inferred from ngsPedigree's genotype/karyotype data and **must be
confirmed against the ngsPedigree source**; ngsTracts only consumes their output
through the tables below.

## Adapter input (required): `polarized_transmissions.tsv`

One row per (parent -> offspring) transmission per inversion, *after* the
polarity has been chosen upstream.

| col | name | type | semantics |
| --- | --- | --- | --- |
| 1 | `inversion_id` | string | candidate id, matches the atlas |
| 2 | `chrom` | string | |
| 3 | `parent_id` | string | the transmitting parent |
| 4 | `offspring_id` | string | |
| 5 | `parent_karyotype` | enum | `HOM_REF` \| `HET` \| `HOM_INV` (Layer 1) |
| 6 | `parent_role` | enum | `father` \| `mother` \| `unknown` (sex stratification) |
| 7 | `transmitted_arrangement` | enum | `REF` \| `INV` \| `ambiguous` (Layer 3) |
| 8 | `informative` | int | `1` if this transmission is usable for the drive test, else `0` |

Only **HET parents** with `informative == 1` and `transmitted_arrangement in
{REF, INV}` contribute to the drive test (a HET parent transmits REF or INV;
HOM parents are fixed and non-informative for 1:1 segregation). Triad-derived,
parent-specific transmissions are the strongest input; dyad-only rows where the
transmitting allele cannot be assigned should be marked `informative = 0` or
`transmitted_arrangement = ambiguous` upstream.

## Adapter input (optional): `inversion_polarity.tsv`

Per-inversion polarization summary computed upstream; passed through unchanged
into the report so the chosen orientation and contradiction counts travel with
the result.

| name | semantics |
| --- | --- |
| `inversion_id` | |
| `chosen_polarity` | retained REF/INV orientation |
| `polarization_score_REF_INV` | Mendelian-compatibility score under one orientation |
| `polarization_score_INV_REF` | score under the flipped orientation |
| `n_dyad_contradictions` | impossible dyad combinations |
| `n_triads_tested` | full trios checked |
| `n_triad_contradictions` | impossible trio combinations |

## What STEP_TRC_05 computes (downstream of the adapter)

Per inversion, into `transmission_drive.tsv`:

- `n_informative_het_parent_transmissions`, `n_REF_transmitted`, `n_INV_transmitted`
- `INV_transmission_rate = n_INV / (n_REF + n_INV)`
- `mendelian_drive_pvalue` — two-sided exact binomial test of `n_INV ~
  Binomial(n, 0.5)` (H0: 1:1 segregation; dependency-free `math.comb`)
- `paternal_*` / `maternal_*` — the same, stratified by `parent_role`
- pass-through polarity/contradiction fields (if `--polarity` supplied)
- `n_CO_candidates` / `n_NCO_candidates` / `n_DCO_candidates` — joined from
  STEP_TRC_01 by inversion overlap (if `--tract-classifications` +
  `--inversion-atlas` supplied)

## Method note (two-generation novelty)

> Inversion karyotype supplies the phase polarity normally obtained from
> grandparental haplotypes in classical three-generation recombination maps.
> Both possible polarizations are evaluated upstream against dyad and triad
> inheritance constraints, and the orientation minimizing Mendelian
> contradictions is retained. ngsTracts then tests the polarized HET-parent
> transmissions against 1:1 segregation and uses the polarized transmitted
> arrangement as the phase framework for CO/NCO/DCO classification.

Everything before the adapter line is external (ngsPedigree); everything after
is in this repository.

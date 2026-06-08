# LRR crossover-suppression ranking with population support (STEP_TRC_06)

Evaluates whether each candidate low-recombination region (LRR) behaves as a
coherent **inherited block**, by asking — across the whole cohort, not one
family — whether dyads that *had the opportunity* to show a crossover inside the
region actually do.

## Why a population score

A single family showing "no CO here" is weak: a region with no informative
transmissions trivially has zero COs. So the headline statistic is a
conservative population support score:

```
population_support_score = WilsonLowerBound_95%( k / n )
  n = dyads with OPPORTUNITY in the LRR interior
  k = those opportunity-dyads with NO interior CO
```

The Wilson lower bound penalizes thin support (one consistent family scores
~0.20, not 1.0) and rises toward 1 only with many consistent dyads. LRRs are
ranked by this score. (`scripts/STEP_TRC_06_lrr_suppression.py:60`,
`scripts/STEP_TRC_06_lrr_suppression.py:158`)

## Opportunity (the denominator that prevents false positives)

A dyad "has opportunity" in an LRR if it has at least
`--min-opportunity-sites` (default 10) informative interior sites
(from the Stage-3 per-site file), **or** it produced any observed interior
event there. Without `--per-site-file`, a weak proxy is used (any event on the
chromosome) and the `opportunity_source` column says so
(`scripts/STEP_TRC_06_lrr_suppression.py:131`). LRRs with fewer than
`--min-opportunity-dyads` (default 2) opportunity dyads are labelled
`low_confidence`, never "suppressed".

## Interior vs edge

Each LRR is split into edge zones (`--edge-frac` of length per end, default
0.10, or `--edge-bp`) and an interior. CO best-estimate switch points
(midpoints) are assigned to interior or edge, separating true internal
suppression from boundary signal (`scripts/STEP_TRC_06_lrr_suppression.py:177`).
NCO/DCO counts are reported but do **not** count against suppression — gene
conversion inside a block is expected gene flux, not a crossover.

## Pattern labels

| pattern | meaning |
| --- | --- |
| `no_CO_inside` | strong inherited-block candidate (no interior CO, none at edges) |
| `CO_at_boundary_only` | interior clean but CO at the edge → possible block boundary |
| `few_CO_inside` | some interior CO (< `--many-co-threshold`) → intermediate |
| `many_CO_inside` | interior CO ≥ threshold (default 2) → weak/fragmented LRR |
| `low_confidence` | too few opportunity dyads (or degenerate interior) |

## Output `lrr_co_suppression.tsv`

Ranked, one row per LRR: `population_support_score`, `suppression_consistency`
(raw k/n), `n_dyads_with_opportunity`, `n_dyads_no_interior_CO`,
`opportunity_source`, `mean_interior_sites_per_dyad`, interior CO/NCO/DCO and
edge CO counts, `pattern`, `flag_low_opportunity`.

## Where it sits

- **Upstream (external, ngsPedigree):** dyads/triads, mtDNA-supported maternal
  directionality, polarized transmissions (see
  `docs/SCHEMA_POLARIZATION_ADAPTER.md`).
- **In-repo:** CO/NCO/DCO calls (STEP_TRC_01/02) → this suppression ranking.
- **External (popstats / unified-ancestry):** GHSL, FIS and divergence. The
  overlay (`high GHSL + no_CO_inside + high support`, `+ negative FIS =
  heterozygote-excess block`) is a join the popstats engine does; STEP_TRC_06
  supplies the CO-suppression side.

> Because only adult broodstock are sampled, this does not separate meiotic
> drive from post-transmission viability. It asks the answerable question:
> across the cohort, do candidate LRRs behave as coherent inherited
> low-recombination blocks, supported by an absence of pedigree-backed CO-like
> events with adequate opportunity?

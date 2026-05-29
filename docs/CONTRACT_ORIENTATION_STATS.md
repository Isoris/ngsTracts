# Contract: `iv_subsample_orientation_stats()` — arrangement-diversity handoff

**Status: interface contract only. NOT implemented in this repository.**

The arrangement-specific population statistics (pi, DAF, dXY, polarization) are
computed by the **external popstats server / unified-ancestry engine**, not by
ngsTracts. This document fixes the boundary: what ngsTracts supplies, what the
popstats engine consumes and returns, so the implementation (in VS Code,
against that engine) can be written against a stable contract.

ngsTracts' role stops at producing the **exchange/gene-flux segment mask**
(`scripts/STEP_TRC_03_exchange_mask.py`) and the per-interval calls
(`STEP_TRC_01`). Everything below the dashed line in the data flow is external.

## Masking model for genotype likelihoods (per-sample, not global)

The popstats engine works on **genotype likelihoods**, not hard calls. A single
pooled BED is therefore the wrong masking unit: dropping a site globally would
discard *every* sample's GL there. The exchange tracts are observed in
parent->offspring transmission, and the **offspring inherits** the exchanged
sequence — so the carrier is the offspring, and masking must be **per sample**:
exclude only the carrier sample's GL contribution at its tract sites, keeping
all other samples at those sites.

STEP_TRC_03 therefore emits, in addition to the pooled BED:

- `by_sample/<sample_id>.exchange.bed` — one BED per carrier sample
  (`--per-sample`; carrier = offspring by default, `--sample-key`).
- `exchange_segments.by_sample.tsv` — long form `sample_id, role, chrom,
  bed_start, bed_end, class, interval_id, parent_id, offspring_id`.

The pooled / merged BED remains useful only as an *exchange-prone region*
annotation (cross-dyad pileup), not as the GL masking unit.

## Provenance

| Component | Provider | In this repo? |
| --- | --- | --- |
| Exchange-tract mask (`exchange_segments*.bed`) | ngsTracts STEP_TRC_03 | yes |
| Per-interval calls (`tract_classifications.tsv`) | ngsTracts STEP_TRC_01 | yes |
| Inversion core/flank/buffer partition | inversion atlas (external) | no |
| Genotypes + AA/AD/DD karyotype assignment | popstats / unified-ancestry | no |
| pi / DAF / dXY / bootstrap / polarization | popstats / unified-ancestry | no |
| `iv_subsample_orientation_stats()` itself | popstats / unified-ancestry | no (this file is the spec) |

## Signature (to be implemented in the popstats engine)

```python
def iv_subsample_orientation_stats(
    *,
    candidate_bed: Path,            # inversion candidates: chrom,start,end,inversion_id
    karyotype_table: Path,          # sample_id, inversion_id, karyotype in {AA, AD, DD}
    genotypes: Path | str,          # call set / handle the popstats server understands
    outgroup: Path | str | None,    # ancestral/outgroup allele source for polarization
    # --- masks supplied by ngsTracts + atlas (layers 1-3 of the strategy) ---
    breakpoint_buffer_bp: int | float = 0.10,   # int bp, or fraction of length if <1
    bad_window_bed: Path | None = None,         # low-mappability/repeat/depth/missingness
    # GL masking is PER SAMPLE: pass the by-sample table, not a single pooled BED.
    exchange_by_sample: Path | None = None,     # STEP_TRC_03 exchange_segments.by_sample.tsv
    exchange_prone_bed: Path | None = None,     # OPTIONAL annotation: merged pooled BED
    exchange_min_support: int = 1,              # min pileup depth to call a region exchange-prone
    # --- windowing & resampling ---
    window_bp: int = 20_000,
    step_bp: int | None = None,                 # default: non-overlapping (= window_bp)
    min_sites_per_window: int = 1,
    n_bootstrap: int = 500,
    seed: int | None = None,
) -> "OrientationStatsResult":
    """Compute arrangement-specific diversity per inversion candidate.

    Reports BOTH an INCLUSIVE core (breakpoint buffer + bad windows masked) and
    a MASKED core (additionally masking exchange/gene-flux tracts), so broad
    haplotype age is separable from localized sequence exchange.

    NOT IMPLEMENTED HERE — compute lives in the popstats / unified-ancestry
    engine. This signature is the contract only.
    """
    raise NotImplementedError("computed by external popstats / unified-ancestry engine")
```

## Inputs

| Arg | Source | Meaning |
| --- | --- | --- |
| `candidate_bed` | inversion atlas | one row per candidate; defines layer-1 geometry |
| `karyotype_table` | popstats / unified-ancestry | per-sample AA/AD/DD per inversion |
| `genotypes` | popstats / unified-ancestry | the call set diversity is computed from |
| `outgroup` | external | ancestral allele for derived-allele polarization |
| `breakpoint_buffer_bp` | strategy layer 1 | trim near breakpoints; bp, or fraction if `<1` |
| `bad_window_bed` | external QC | layer-3 mappability/repeat/depth/missingness mask |
| `exchange_by_sample` | **ngsTracts STEP_TRC_03** | layer-2, **per-sample** tracts (GL masking unit) |
| `exchange_prone_bed` | ngsTracts STEP_TRC_03 | optional pooled annotation only |
| `exchange_min_support` | — | merged BED col 4 (pileup depth) threshold for "exchange-prone" |
| `window_bp`/`step_bp` | strategy | subsampling unit is the window, not the raw SNP |
| `n_bootstrap`/`seed` | strategy | resample windows with replacement; CI + reproducibility |

The masking layers map exactly to the agreed strategy:
- **Layer 1 (core/flank/buffer):** `candidate_bed` + `breakpoint_buffer_bp`.
- **Layer 2 (exchange/meiosis tracts):** `exchange_by_sample` from STEP_TRC_03,
  applied **per carrier sample** (offspring); `exchange_prone_bed` +
  `exchange_min_support` give an optional pooled exchange-prone annotation.
- **Layer 3 (technical bad windows):** `bad_window_bed`.

## Output (`OrientationStatsResult`, one row per candidate)

```
inversion_id, chrom, core_start, core_end, core_length,
n_valid_windows,
n_masked_breakpoint, n_masked_low_quality, n_masked_exchange_like,
piAA_inclusive,  piDD_inclusive,  piDD_piAA_inclusive,
piAA_masked,     piDD_masked,     piDD_piAA_masked,
dXY_AA_DD, DAF_delta,
polarization_call, diversity_interpretation,
bootstrap_ci_piDD_piAA_low, bootstrap_ci_piDD_piAA_high,
support, warnings
```

Dual-estimate convention (decision rule, computed by the engine):
- inclusive ≈ masked → robust.
- inclusive high but masked low → localized gene-flux tracts (the exchange mask
  was load-bearing).
- both high → derived arrangement old / introgressed / mis-polarized.

Masking aggressiveness differs by test (per the strategy): **moderate** masking
for polarization, **strict** (include exchange mask) for founding-diversity/age.

## What ngsTracts guarantees to the engine

1. `by_sample/<sample_id>.exchange.bed` + `exchange_segments.by_sample.tsv`
   (`--per-sample`) — **the GL masking unit**: per carrier sample (offspring),
   BED6 0-based half-open, name `CLASS|interval_id|parent>offspring`.
2. `exchange_segments.bed` — pooled BED6 (all samples), for inspection.
3. `exchange_segments.merged.bed` — cross-dyad merged; column 4 = support
   (pileup depth) for `exchange_min_support` (exchange-prone annotation only).
4. `exchange_mask_summary.tsv` — per-class and (with `--inversion-atlas`)
   per-inversion tract counts and masked bp.

Coordinates are 0-based half-open in the BEDs (converted from the 1-based
inclusive interval calls); the engine must treat them as BED.

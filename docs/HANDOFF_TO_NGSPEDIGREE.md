# Handoff: what ngsPedigree Stage 3 needs to produce

This document is for the **ngsPedigree** repo, not ngsTracts. It's the
contract you implement on the Stage 3 side. Copy it into ngsPedigree's
docs when you start work on `STEP_PED_03_inheritance_map.py`.

## TL;DR

Stage 3 must produce four files in its output directory:

1. `departure_intervals.tsv` — runs of sites disagreeing with inferred background
2. `chromosome_background_intervals.tsv` — inferred background blocks per dyad
3. `stage3.args.tsv` — provenance + schema version
4. `per_site_haplotype_calls.tsv` — optional, only when `--emit-per-site`

Plus the dependencies it consumes:

- BEAGLE GLs (whole-genome merged, thin500) — same one used for ngsRelate
- ngsPedigree Stage 1 family roster (constraint vocab)
- Inversion atlas (optional, for `inside_inversion` tagging)
- Sample list

Schemas are normative — see `SCHEMA.md` for column-level details.

## Methodology summary (what Stage 3 actually computes)

For each (parent, offspring) dyad classified by Stage 1 as PO or candidate-PO,
and for each chromosome:

### A. Per-chromosome background inference

1. Restrict to **informative sites**: parent het (hard-called from GLs at
   threshold 0.90), offspring genotyped at threshold 0.95.
2. Use a per-chromosome Viterbi-style HMM with two hidden states (`hapA`,
   `hapB`) corresponding to which parental haplotype was transmitted.
3. Emit one row per inferred background block to
   `chromosome_background_intervals.tsv`.

The HMM transition probability sets the recombination prior; tune for
catfish-scale chromosomes (~30 Mb LGs, ~1-2 COs expected per chrom).

### B. Departure interval detection

A "departure" is a run of ≥3 consecutive informative sites whose
offspring genotype disagrees with the inferred background. Runs of 1-2
are dropped as singleton noise.

For each detected departure:
- `start`, `end` = first and last site position in the run
- `n_sites` = total informative sites in [start, end]
- `n_discordant` = subset disagreeing with background
- `flanking_left_state`, `flanking_right_state` = background state of the
  block immediately to the left/right
- `departure_state` = which haplotype the run agrees with (`hapA` |
  `hapB` | `neither` if mixed)

### C. Inversion tagging (if atlas supplied)

For each departure interval:
- `inside_inversion = yes` iff [start, end] is fully contained in an
  inversion region from the atlas
- `inside_inversion = partial` iff [start, end] crosses an inversion boundary
- `inside_inversion = no` otherwise
- `distance_to_nearest_inv_bp` = bp to nearest inversion boundary on
  the same chrom, or `-` if no inversions on this chrom

### D. Per-interval confidence

`confidence` should reflect:
- Informative-site density in/around the interval
- Hard-call rate (fraction of sites in the dyad that hit the 0.90/0.95
  thresholds vs. fell into ambiguous probability ranges)
- BAQ flag if available

Suggested mapping: `high` if all of (density >= 100 sites/Mb, hard-call
rate >= 0.95, no BAQ issues). `low` if any of (density < 50 sites/Mb,
hard-call rate < 0.85). `medium` otherwise.

### E. Provenance

`stage3.args.tsv` must contain at minimum:

```
schema_version          0.1
stage3_version          0.1
hard_call_threshold_hom 0.95
hard_call_threshold_het 0.90
min_run_length_sites    3
inversion_atlas_supplied yes|no
beagle_input            /path/to/input.beagle.gz
sample_list             /path/to/samples.txt
emit_per_site           yes|no
n_intervals             <count>
n_background_blocks     <count>
datetime                <ISO 8601 with timezone>
host                    <hostname>
```

## What ngsTracts does NOT need

Stage 3 does NOT need to:

- compute recombination rates (ngsTracts does this from CO/NCO counts)
- classify events as NCO/CO/DCO (that's the entire job of ngsTracts)
- decide what's "valid" — ngsTracts has its own confidence gate
- localize CO breakpoints to single-bp resolution (STEP_TRC_02 does that
  if Stage 3 emits per-site calls)

Keep Stage 3 focused on: "given these GLs and this pedigree, where does
the offspring's transmitted haplotype come from at each chromosome
position, and where does it diverge from that expectation?"

## Synthetic fixture as reference implementation

`tests/fixture_basic/` in this repo is a working example of all four
output files. If your Stage 3 implementation can produce files that
the ngsTracts smoke test accepts when pointed at them, you've got the
contract right.

Verify your Stage 3 outputs with:

```bash
python3 scripts/STEP_TRC_01_classify_intervals.py \
    --stage3-dir /path/to/your/stage3/out \
    --fai /path/to/ref.fa.fai \
    --outdir /tmp/test_handoff
```

If it runs to completion without errors, the handoff works. If it
errors, the error message points at the specific column or constraint
violation.

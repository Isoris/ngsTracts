# ngsTracts (v0.1.0)

Classifier for parent–offspring haplotype departures detected by
[ngsPedigree](https://github.com/...) Stage 3. Assigns each departure
interval to one of:

```
NCO              non-crossover (gene conversion)
CO               single crossover
DCO              double crossover (return-to-flank, 50–200kb)
MOSAIC_SHORT     gene conversion of unusual size (50–200kb inside inversion)
MOSAIC_LONG      suspicious long-span departure with matching flanks
AMBIG            data shape inconsistent with any single class
LOW_CONFIDENCE   below quality gate
```

Plus per-dyad and per-chrom event-rate summaries. Optional second pass
refines CO breakpoints to single-bp resolution via a sliding-window scan
over per-site haplotype calls.

This repo does NOT phase, call genotypes, or assign haplotype background.
That's ngsPedigree Stage 3. ngsTracts is the analytical layer that turns
Stage 3's intervals into biological events.

---

## Status

- **Stage 3 schema**: 0.1 (locked, see `docs/SCHEMA.md`)
- **ngsTracts version**: 0.1.0
- **STEP_TRC_01**: implemented, vectorized, tested on synthetic fixture
- **STEP_TRC_02**: implemented, tested on synthetic fixture (CO breakpoint
  refined to 500 bp of truth on the test case)
- **Real-data run**: pending Stage 3 producing real outputs on the 226-sample
  hatchery cohort

---

## Quickstart

```bash
# 1. Run smoke test (verifies install)
python3 tests/smoke_test.py

# 2. Classify intervals (production use)
python3 scripts/STEP_TRC_01_classify_intervals.py \
    --stage3-dir /scratch/.../ngsPedigree/stage3_out \
    --fai        /scratch/.../fClaHyb_Gar_LG.fa.fai \
    --outdir     /scratch/.../ngsTracts/out

# 3. Optional: refine CO breakpoints via per-site scan
python3 scripts/STEP_TRC_02_traversal_scan.py \
    --tract-classifications /scratch/.../ngsTracts/out/tract_classifications.tsv \
    --per-site-file         /scratch/.../ngsPedigree/stage3_out/per_site_haplotype_calls.tsv \
    --outdir                /scratch/.../ngsTracts/out
```

Both scripts run in minutes on a login node. No SLURM job needed.

---

## Outputs

Written to `--outdir`:

```
tract_classifications.tsv    per-interval classification (PRIMARY OUTPUT)
dyad_event_rates.tsv         per (parent, offspring, chrom) summary
chrom_summary.tsv            per-chrom aggregate across dyads
ngstracts.args.tsv           provenance + parameters
traversal_breakpoints.tsv    (if STEP_TRC_02 was run) refined CO breakpoints
```

For column-level schema of each output, see `docs/METHODOLOGY.md` §5.

---

## Repo layout

```
ngsTracts/
├── README.md
├── docs/
│   ├── METHODOLOGY.md            decision tree, parameters, biological reasoning
│   ├── SCHEMA.md                 column-level input/output contract
│   └── SCHEMA_COMPATIBILITY.md   Stage 3 ↔ ngsTracts version matrix
├── scripts/
│   ├── STEP_TRC_01_classify_intervals.py   main classifier
│   └── STEP_TRC_02_traversal_scan.py       optional CO breakpoint refiner
├── tests/
│   ├── make_fixture.py           synthetic Stage 3 fixture builder
│   ├── smoke_test.py             end-to-end smoke test
│   └── fixture_basic/            synthetic test data
├── examples/                     example invocation scripts
└── lanta/                        LANTA-specific launchers (optional)
```

---

## Input requirements (from ngsPedigree Stage 3)

Two required files in `--stage3-dir`:

1. `departure_intervals.tsv` — the PRIMARY input (intervals to classify)
2. `stage3.args.tsv` — provenance, must contain `schema_version` line

One optional file (only for STEP_TRC_02):

3. `per_site_haplotype_calls.tsv` — produced by Stage 3 when run with
   `--emit-per-site`

Full column specs in `docs/SCHEMA.md`. The schema is **normative** — Stage 3
must produce these columns with these names, types, and constraints, or
ngsTracts refuses to run.

---

## How the classifier works (5-line summary)

For each departure interval:
1. Confidence gate (drop low-quality / few-site intervals).
2. Inside-inversion check (short→NCO, mid→MOSAIC_SHORT, long→MOSAIC_LONG).
3. Outside-inversion short interval → NCO if flanks match, AMBIG otherwise.
4. Long interval with flank switch → CO.
5. Long interval with matching flanks → DCO (if 50-200kb, high disc frac) or MOSAIC_LONG.

Full decision tree in `docs/METHODOLOGY.md` §3.

---

## Dependencies

- Python 3.9+
- `numpy`
- `pandas`

No bioinformatics tools, no SLURM, no compiled extensions.

---

## Schema compatibility

| ngsPedigree Stage 3 schema | ngsTracts version | OK |
| -------------------------- | ----------------- | -- |
| 0.1                        | 0.1.x             | ✓  |

ngsTracts reads `schema_version` from `stage3.args.tsv` on startup and
refuses to run on a mismatch. See `docs/SCHEMA_COMPATIBILITY.md` for
upgrade rules.

---

## Handoff to ngsPedigree

If you maintain ngsPedigree and are implementing Stage 3, the
authoritative input spec is `docs/SCHEMA.md`. Specifically:

- Required columns in `departure_intervals.tsv` (15 of them)
- Allowed enum values per column
- Required constraints (`n_discordant <= n_sites`, etc.)
- `interval_id` format: `DEP_NNNNNN` with NNNNNN zero-padded to width ≥6
- `stage3.args.tsv` must contain `schema_version`

The `tests/fixture_basic/` directory is a working example of the expected
output format. If your Stage 3 implementation produces files that pass
`STEP_TRC_01` without errors on the schema check, you've got the contract
right.

# ngsPedigree Stage 3 ↔ ngsTracts handoff schema (v0.1)

This document is the **machine-readable contract**. Anything in this file
is part of the API boundary between the two repos.

## File: `chromosome_background_intervals.tsv`

| col | name                  | type     | required | semantics                                                  |
| --- | --------------------- | -------- | -------- | ---------------------------------------------------------- |
|   1 | `parent_id`           | string   | yes      | sample ID, matches ngsRelate `ida`/`idb`                   |
|   2 | `offspring_id`        | string   | yes      | sample ID                                                  |
|   3 | `chrom`               | string   | yes      | chromosome name, e.g. `C_gar_LG01`                         |
|   4 | `start`               | int      | yes      | 1-based inclusive start of background block (bp)           |
|   5 | `end`                 | int      | yes      | 1-based inclusive end of background block (bp)             |
|   6 | `background_state`    | enum     | yes      | `hapA` \| `hapB` \| `unphased`                             |
|   7 | `n_inf_sites`         | int      | yes      | informative sites in block (parent het, offspring hard)    |
|   8 | `concordant_frac`     | float    | yes      | sites matching background ÷ n_inf_sites (excl. departures) |
|   9 | `notes`               | string   | no       | free text; `-` if none                                     |

Sort order: `parent_id, offspring_id, chrom, start`.
Coordinates: 1-based, inclusive. Same convention as ANGSD/samtools.

## File: `departure_intervals.tsv`

The primary input to ngsTracts. Schema is normative.

| col | name                          | type     | required | semantics                                                       |
| --- | ----------------------------- | -------- | -------- | --------------------------------------------------------------- |
|   1 | `interval_id`                 | string   | yes      | format `DEP_NNNNNN` (6+ digits, stable across runs)             |
|   2 | `parent_id`                   | string   | yes      |                                                                 |
|   3 | `offspring_id`                | string   | yes      |                                                                 |
|   4 | `chrom`                       | string   | yes      |                                                                 |
|   5 | `start`                       | int      | yes      | 1-based inclusive                                               |
|   6 | `end`                         | int      | yes      | 1-based inclusive                                               |
|   7 | `n_sites`                     | int      | yes      | informative sites in interval                                   |
|   8 | `n_discordant`                | int      | yes      | sites disagreeing with background                               |
|   9 | `span_bp`                     | int      | yes      | `end - start + 1`                                               |
|  10 | `flanking_left_state`         | enum     | yes      | `hapA` \| `hapB` \| `boundary`                                  |
|  11 | `flanking_right_state`        | enum     | yes      | `hapA` \| `hapB` \| `boundary`                                  |
|  12 | `departure_state`             | enum     | yes      | `hapA` \| `hapB` \| `neither`                                   |
|  13 | `distance_to_nearest_inv_bp`  | int / -  | no       | bp to nearest inversion boundary; `-` if no inversions on chrom |
|  14 | `inside_inversion`            | enum     | yes      | `yes` \| `partial` \| `no`                                      |
|  15 | `confidence`                  | enum     | yes      | `high` \| `medium` \| `low`                                     |

Constraints (Stage 3 must enforce):
- `start <= end`
- `span_bp == end - start + 1`
- `n_discordant <= n_sites`
- `n_sites >= 3` (intervals below this are not emitted; runs of 1-2
  discordant sites are filtered as singleton noise)
- `inside_inversion = yes` only when the entire interval `[start, end]`
  is contained in an annotated inversion; `partial` when it crosses a
  boundary; `no` otherwise (including when no atlas is supplied)
- When no inversion atlas is supplied: `inside_inversion` always `no`,
  `distance_to_nearest_inv_bp` always `-`
- `interval_id` globally unique across all dyads and chromosomes in a
  single Stage 3 run; format `DEP_NNNNNN` with `NNNNNN` zero-padded to
  width >=6

## File: `per_site_haplotype_calls.tsv` (optional)

Large file (one row per informative site per dyad). Emitted only when
Stage 3 is invoked with `--emit-per-site`. ngsTracts only reads this
when `STEP_TRC_02_traversal_scan.py` is run.

| col | name              | type   | required | semantics                                       |
| --- | ----------------- | ------ | -------- | ----------------------------------------------- |
|   1 | `parent_id`       | string | yes      |                                                 |
|   2 | `offspring_id`    | string | yes      |                                                 |
|   3 | `chrom`           | string | yes      |                                                 |
|   4 | `pos`             | int    | yes      | 1-based                                         |
|   5 | `parent_gt`       | enum   | yes      | `0/0` \| `0/1` \| `1/1`                         |
|   6 | `offspring_gt`    | enum   | yes      | `0/0` \| `0/1` \| `1/1`                         |
|   7 | `inferred_state`  | enum   | yes      | `hapA` \| `hapB` (locally inferred background)  |
|   8 | `departure`       | enum   | yes      | `0` \| `1`                                      |
|   9 | `interval_id`     | string | no       | `DEP_NNNNNN` if site is in an interval; `-` else |

Sort order: `parent_id, offspring_id, chrom, pos`.

## File: `stage3.args.tsv` (provenance)

Stage 3 must emit this alongside the data files. ngsTracts reads it on
startup to verify schema compatibility.

```
key                            value
schema_version                 0.1
stage3_version                 0.1
n_dyads                        12
n_chroms                       28
inversion_atlas_supplied       yes
inversion_atlas_path           /scratch/.../inversion_atlas.tsv
hard_call_threshold_hom        0.95
hard_call_threshold_het        0.90
min_run_length_sites           3
beagle_input                   /scratch/.../catfish.wholegenome.byRF.thin_500.beagle.gz
sample_list                    /scratch/.../list_of_samples_one_per_line_same_bamfile_list.tsv
emit_per_site                  no
n_intervals                    8421
n_background_blocks            142
datetime                       2026-05-12T08:31:00+07:00
host                           x3002c0s7b0n0
```

`schema_version` is the line ngsTracts checks. If it doesn't match
ngsTracts' supported set (see `SCHEMA_COMPATIBILITY.md`), ngsTracts
refuses to run with a clear error message.

## File: `inversion_atlas.tsv` (optional, supplied externally)

Provided to BOTH Stage 3 and ngsTracts. Source: MS_Inversions
supplementary atlas.

| col | name              | type   | required | semantics                                          |
| --- | ----------------- | ------ | -------- | -------------------------------------------------- |
|   1 | `inversion_id`    | string | yes      | e.g. `INV_028`                                     |
|   2 | `chrom`           | string | yes      |                                                    |
|   3 | `start`           | int    | yes      | 1-based inclusive                                  |
|   4 | `end`             | int    | yes      | 1-based inclusive                                  |
|   5 | `karyotype_0`     | int    | yes      | n samples homozygous for state 0 (reference)       |
|   6 | `karyotype_1`     | int    | yes      | n heterozygous                                     |
|   7 | `karyotype_2`     | int    | yes      | n homozygous for state 1 (alternate)               |
|   8 | `fst_hom1_hom2`   | float  | yes      | Fst between karyotype groups 0 and 2               |

The atlas is shared by both ngsPedigree (for tagging `inside_inversion`
in `departure_intervals.tsv`) and ngsTracts (for the inversion-aware
priors in `STEP_TRC_01`).

## Per-parent karyotype lookup (optional)

If a per-sample karyotype table is supplied to ngsTracts via
`--parent-karyotypes`, the inside-inversion NCO rule (METHODOLOGY §3.2)
can be applied more precisely — it only fires when the **parent is known
heterozygous** for the inversion containing the interval.

Format: `inversion_id\tsample_id\tkaryotype` where karyotype is `0`/`1`/`2`.
If not supplied, ngsTracts assumes het for any parent with an interval
inside an inversion (less conservative — sets `manual_review_flag` for
those calls).

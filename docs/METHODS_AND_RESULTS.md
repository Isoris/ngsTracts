# Classification of parent–offspring haplotype-background departures into recombination and gene-conversion events: methods and results

## Provenance note (read first)

This manuscript describes a two-stage analytical chain. The two stages have
different provenance, and the reader must keep them separate:

- **Upstream / external engine (NOT in this repository).** All raw-sequence
  encoding, genotype-likelihood computation, hard genotype calling,
  informative-site selection, per-chromosome haplotype-background inference,
  departure-interval detection, inversion tagging, and per-interval confidence
  assignment are performed by an upstream caller ("Stage 3",
  `STEP_PED_03_inheritance_map.py`) whose source is **not present here**. This
  repository only consumes its tabular output through a locked schema
  (`docs/SCHEMA.md`). Where I state the estimator that the upstream stage
  *must* be using, that estimator is **inferred from the data contract and the
  handoff document** (`docs/HANDOFF_TO_NGSPEDIGREE.md`) and **must be confirmed
  against the upstream source**, which I could not read.

- **In-repository compute (verifiable here).** Everything from the integer
  interval encoding onward — the discordance statistic, the comparison
  primitives, the deterministic decision tree that assigns each interval to a
  recombination/gene-conversion class, the per-call confidence and review
  flags, the per-group event-rate statistics, the cross-group aggregation, and
  the sliding-window crossover-breakpoint refinement — is implemented in
  `scripts/STEP_TRC_01_classify_intervals.py` and
  `scripts/STEP_TRC_02_traversal_scan.py`. Every formula, threshold, and
  decision rule below cites its implementing `path:line`.

Three documented behaviours are **specified but not wired into the executed
code**, and are flagged inline and in Appendix B so they are not mistaken for
verified method:

1. The karyotype-aware refinement of the inside-inversion gene-conversion rule
   is loaded and validated but its return value is **discarded** and never
   reaches the classifier (`scripts/STEP_TRC_01_classify_intervals.py:386`).
2. The position-dependent prior is computed but is **informational only** and
   never changes a call (`scripts/STEP_TRC_01_classify_intervals.py:278`;
   `docs/METHODOLOGY.md:113`).
3. The genome-wide summary row (`chrom == "GENOME"`) promised in
   `docs/METHODOLOGY.md:169` is **not produced** by the aggregation code
   (`scripts/STEP_TRC_01_classify_intervals.py:290`).

Finally: the headline cohort (a 226-sample hatchery cohort, named in the
launcher `lanta/run_on_lanta.sh:9`) **has not been run** — the upstream stage
has not yet produced real outputs for it (`README.md:35`). The only data with
real result artifacts on disk is the synthetic conformance fixture
(`tests/fixture_basic/`), which I executed for the Results section. Every
number in Results is taken verbatim from those artifacts and is explicitly
labelled as synthetic.

---

## Abstract

We classify runs of sites at which a transmitted offspring haplotype departs
from its inferred single-parent haplotype background into the biological event
that produced them: non-crossover gene conversion, single crossover, double
crossover, or one of several mosaic/ambiguous residual classes. The input is a
per-interval table emitted by an upstream parent–offspring inheritance caller;
each interval carries its genomic span, its informative-site count, its count
of background-discordant sites, the haplotype state of the flanking background
blocks, and an inside-inversion flag. The chain computes a single per-interval
discordance fraction, derives two boolean comparison primitives over the
flanking states, and runs a deterministic first-match decision tree whose
thresholds are length scales (50 kb, 200 kb, 1 Mb) and a strong-discordance
cut (0.9); it then aggregates calls per (parent, offspring, chromosome) dyad
into crossover and non-crossover rates per megabase and a crossover-to-gene-
conversion ratio, and across dyads into per-chromosome means. An optional
second pass refines crossover breakpoints to the position of maximum
haplotype-concordance gradient in a 21-site sliding window. On the only data
with artifacts on disk — an 8-interval synthetic conformance fixture across two
dyads and two linkage groups — the classifier returns 2 non-crossover, 1
crossover, 1 double-crossover, 1 short-mosaic, 1 long-mosaic, 1 ambiguous and 1
low-confidence call, flags 4/8 intervals for manual review, and the refinement
pass localizes the single synthetic crossover to 19,949,500 bp against a
ground truth of 19,950,000 bp (500 bp error). **No result exists for the
226-sample cohort; the cohort run is pending upstream output.**

---

## 1. Introduction

The analytical problem is to read meaning out of *departures*: contiguous runs
of informative sites where an offspring's genotype is inconsistent with the
parental haplotype it was otherwise inferred to have inherited. A departure is
a symptom; its cause may be (i) a true crossover, where the transmitted
chromosome switches from one parental haplotype to the other and stays
switched; (ii) a double crossover, where it switches and switches back over a
short span; (iii) a non-crossover gene-conversion tract, where a short patch is
copied from the homolog without an exchange of flanks; or (iv) an artefact —
mis-phasing, a mis-assigned background block, or a structurally suppressed
region such as the interior of a segregating inversion, where ordinary
crossing-over does not occur and a short departure is far more likely to be
gene conversion than recombination.

This is hard for three reasons. First, the classes are distinguished not by the
departure itself but by its *context* — the haplotype identity of the blocks on
either side, and whether they match. A run with mismatched flanks is a
crossover; the same-length run with matching flanks is gene conversion or a
double crossover, and which one depends on length and discordance density.
Second, length scales are biological priors, not measurements: gene-conversion
tracts are short (sub-50 kb), double crossovers that "return to flank" are
plausible only over an intermediate window, and very long matched-flank
departures are more likely to be background-inference errors than real events.
Third, recombination suppression inside inversions inverts the default
interpretation, so the same data shape must be read differently depending on an
externally supplied structural annotation. The chain therefore encodes its
biology as an ordered set of thresholded rules over a handful of robust
per-interval statistics, keeps the calling deterministic, and quarantines every
soft or model-based quantity (priors, karyotype refinement) out of the calling
path so that a reported call is always reproducible from the input row alone.

The remainder follows the data flow. Section 2 covers the input encoding,
windowing, and the schema gate that admits data. Section 3 defines the core
per-interval detection statistic. Section 4 defines the reusable comparison
primitives over flanking state. Section 5 is the candidate-construction
decision tree (the calling step). Section 6 covers the per-call confidence and
manual-review logic. Section 7 is the informational position prior. Section 8
is the per-dyad grouping and population statistics; Section 9 the cross-dyad
per-chromosome aggregation. Section 10 is the inheritance/relatedness context
and the (unwired) karyotype refinement. Section 11 is the crossover-breakpoint
refinement pass. Section 12 covers serialization and reproducibility.

---

## 2. Input encoding, windowing, and the schema gate

### 2.1 Upstream encoding (external; confirm against engine)

The biological signal originates as per-site genotype likelihoods that the
upstream stage hard-calls and then reduces to intervals. The handoff document
states the upstream estimator chain as: restrict to informative sites (parent
heterozygous, offspring genotyped), hard-call from genotype likelihoods at a
het threshold of 0.90 and a homozygote threshold of 0.95, infer the
transmitted background with a **two-state Viterbi HMM** over {hapA, hapB} per
chromosome, and emit a departure wherever ≥3 consecutive informative sites
disagree with that background (`docs/HANDOFF_TO_NGSPEDIGREE.md:31`,
`docs/HANDOFF_TO_NGSPEDIGREE.md:44`, `docs/SCHEMA.md:96`). **These estimators
and thresholds live in the upstream `STEP_PED_03` source, which is not in this
repository; treat them as the contract-implied standard and confirm against
that source.** This repository never sees a likelihood or a base — it begins at
the integer interval table.

### 2.2 In-repo encoding

Each departure is encoded as one row with 15 required fields: identity
(`interval_id`, `parent_id`, `offspring_id`, `chrom`), genomic extent (`start`,
`end`, `span_bp`, all 1-based inclusive), site counts (`n_sites`,
`n_discordant`), flanking context (`flanking_left_state`,
`flanking_right_state` ∈ {hapA, hapB, boundary}), the departing haplotype
(`departure_state`), and structural context (`distance_to_nearest_inv_bp`,
`inside_inversion` ∈ {yes, partial, no}, `confidence` ∈ {high, medium, low})
(`scripts/STEP_TRC_01_classify_intervals.py:88`; column contract
`docs/SCHEMA.md:27`).

On load, span, count, and ordering fields are cast to 64-bit integers and the
inversion distance is coerced to a float with the sentinel `"-"` mapped to NaN
(`scripts/STEP_TRC_01_classify_intervals.py:102`). Three hard integrity
constraints are enforced and abort the run on violation: `span_bp == end −
start + 1`, `n_discordant ≤ n_sites`, and `start ≤ end`
(`scripts/STEP_TRC_01_classify_intervals.py:112`). Five categorical fields are
validated against closed vocabularies, and `interval_id` must be unique
(`scripts/STEP_TRC_01_classify_intervals.py:120`,
`scripts/STEP_TRC_01_classify_intervals.py:132`). Chromosome lengths used by
the rate statistics are read from a FASTA index
(`scripts/STEP_TRC_01_classify_intervals.py:139`).

### 2.3 Schema gate

Before any data file is opened, the provenance file's `schema_version` is
parsed to MAJOR.MINOR and checked against the supported set `{"0.1"}`; a
mismatch aborts with both versions reported
(`scripts/STEP_TRC_01_classify_intervals.py:26`,
`scripts/STEP_TRC_01_classify_intervals.py:64`,
`scripts/STEP_TRC_01_classify_intervals.py:73`). This makes the data contract
load-bearing: the analysis refuses to run on an unrecognized upstream encoding
rather than silently misinterpreting columns.

### 2.4 Windowing

There is no windowing in the calling stage; intervals are the atomic unit.
Windowing appears only in the optional breakpoint-refinement pass (Section 11),
which slides a fixed-width window of informative sites along the per-site
sequence.

---

## 3. Core detection statistic: the discordance fraction

The single per-interval statistic that drives confidence and the
double-crossover gate is the discordance fraction

```
disc_frac = n_discordant / n_sites      (0 if n_sites == 0)
```

computed vectorized over all intervals with a guarded division
(`scripts/STEP_TRC_01_classify_intervals.py:187`). It is the fraction of
informative sites within the interval that actually disagree with the inherited
background. A high value (near 1) means the departure is "clean" — almost every
informative site in the span supports the alternative haplotype — and is used
both as a high-confidence criterion for short tracts and as a necessary
condition for a double-crossover call (Section 5).

---

## 4. Reusable comparison primitives over flanking state

Two boolean vectors over the flanking states are computed once and reused by
every branch of the decision tree
(`scripts/STEP_TRC_01_classify_intervals.py:189`):

```
flank_same  = (flanking_left_state == flanking_right_state)
flanks_real = (flanking_left_state != "boundary") AND
              (flanking_right_state != "boundary")
```

`flank_same` is the topological discriminator between an *exchange* (crossover:
the chromosome enters in one haplotype and leaves in the other) and a *return*
(gene conversion or double crossover: it enters and leaves in the same
haplotype). `flanks_real` guards against calling a crossover when one flank is a
chromosome/assembly boundary rather than a genuine haplotype block, in which
case the switch cannot be interpreted (`scripts/STEP_TRC_01_classify_intervals.py:233`).

---

## 5. Candidate construction: the deterministic decision tree (calling)

Calling assigns each interval exactly one class from {NCO, CO, DCO,
MOSAIC_SHORT, MOSAIC_LONG, AMBIG, LOW_CONFIDENCE}
(`scripts/STEP_TRC_01_classify_intervals.py:48`). The tree is **first-match**:
each rule writes only into rows not yet finalized, tracked by an `assigned`
mask, and the default before any rule is `AMBIG`
(`scripts/STEP_TRC_01_classify_intervals.py:167`,
`scripts/STEP_TRC_01_classify_intervals.py:174`). The ordered rules are:

**(0) Confidence/coverage gate.** If upstream `confidence == low` or
`n_sites < min_n_sites_per_interval` (default 3) → `LOW_CONFIDENCE`
(`scripts/STEP_TRC_01_classify_intervals.py:193`).

**(1) Inside inversion** (`inside_inversion == yes`), where ordinary
crossing-over is suppressed so a departure is read as gene conversion, sized by
span (`scripts/STEP_TRC_01_classify_intervals.py:199`):
- `span_bp < short_tract_max_bp` (default 50 kb) → `NCO`;
- `short_tract_max_bp ≤ span_bp < mosaic_short_max_bp` (50 kb–200 kb) →
  `MOSAIC_SHORT` (review-flagged);
- `span_bp ≥ mosaic_short_max_bp` (≥200 kb) → `MOSAIC_LONG` (review-flagged).

**(2) Short tract outside an inversion**
(`inside_inversion == no` and `span_bp < short_tract_max_bp`)
(`scripts/STEP_TRC_01_classify_intervals.py:219`):
- if `flank_same` → `NCO` (classical short gene-conversion tract);
- else → `AMBIG` (a short span cannot host a real crossover, so mismatched
  flanks are internally inconsistent).

**(3) Flank switch → crossover.** If `not flank_same` and `flanks_real` →
`CO`, with the best-estimate breakpoint taken as the interval midpoint (the
interval is the resolution limit unless Section 11 is run)
(`scripts/STEP_TRC_01_classify_intervals.py:233`;
`docs/METHODOLOGY.md:91`).

**(4) Long matched-flank return.** If `flank_same` and
`span_bp ≥ short_tract_max_bp`
(`scripts/STEP_TRC_01_classify_intervals.py:241`):
- `DCO` iff `span_bp ≤ mosaic_short_max_bp` **and**
  `disc_frac ≥ discordant_frac_strong` (default 0.9) — an intermediate-length,
  clean return-to-flank consistent with two nearby crossovers;
- otherwise `MOSAIC_LONG` (review-flagged).

**(5) Default.** Anything still unassigned remains `AMBIG`, low confidence,
review-flagged (`scripts/STEP_TRC_01_classify_intervals.py:254`).

The calling depends only on the input row and the length/fraction thresholds;
it is deterministic and uses no random number generator, no model fit, and no
position prior.

---

## 6. Per-call confidence and manual-review flagging

**Confidence** is assigned per branch as a categorical {high, medium, low}
(`scripts/STEP_TRC_01_classify_intervals.py:176` and within each rule):
inside-inversion NCO is `high` when `n_sites ≥ n_sites_high_conf_min` (default
5) else `medium` (`scripts/STEP_TRC_01_classify_intervals.py:203`); a short
outside-inversion NCO is `high` only when it is *strong*, i.e.
`disc_frac ≥ discordant_frac_strong` and `n_sites ≥ n_sites_high_conf_min`
(`scripts/STEP_TRC_01_classify_intervals.py:222`); a crossover is graded by how
finely the interval localizes it — `high` for `span_bp ≤ 200 kb`, `medium` up
to `co_breakpoint_max_resolution` (default 1 Mb), `low` beyond
(`scripts/STEP_TRC_01_classify_intervals.py:235`); DCO is fixed at `medium`
(intrinsically rare; `scripts/STEP_TRC_01_classify_intervals.py:244`); mosaics
and the unmatched default are `low`.

**Manual-review flag** (`manual_review_flag` ∈ {0,1}) is set when any of:
DCO with `span_bp > implausible_dco_min_bp` (default 5 Mb); DCO inside an
inversion (crossing-over is suppressed there, so a double crossover is
impossible); NCO with `span_bp > implausible_nco_min_bp` (default 100 kb,
longer than a typical gene-conversion tract); any `MOSAIC_*`; any low-confidence
call; or `n_sites < n_sites_high_conf_min`
(`scripts/STEP_TRC_01_classify_intervals.py:259`). The flag is advisory and
does not alter summary statistics; downstream rate computations are expected to
filter on `class ∈ {CO, NCO} AND confidence != low AND manual_review_flag == 0`
(`docs/METHODOLOGY.md:202`).

---

## 7. Position-dependent prior (informational only — not used in calling)

A per-interval log prior-ratio favouring crossover over gene conversion as a
function of fractional chromosome position is emitted as a column but, by
design, **never changes a call**
(`scripts/STEP_TRC_01_classify_intervals.py:278`; `docs/METHODOLOGY.md:113`).
With `f = midpoint / chrom_len`:

```
prior_dco(f)      = prior_max_dco_rate * 4 * f * (1 - f)
prior_log_ratio   = log( max(prior_dco, 1e-12) / prior_gc_constant )
```

(`scripts/STEP_TRC_01_classify_intervals.py:266`). Defaults are
`prior_max_dco_rate = 0.01` and `prior_gc_constant = 0.001`
(`scripts/STEP_TRC_01_classify_intervals.py:43`). Note that the implemented
form `4·f·(1−f)` is a **parabola peaking at the chromosome midpoint** (value =
`prior_max_dco_rate` at `f = 0.5`, falling to 0 at the telomeres); the prose
label "logistic" in `docs/METHODOLOGY.md:117` is a misnomer — the source is
quadratic. At the midpoint this gives `log(0.01/0.001) = log 10 ≈ 2.3026`,
matching the executed output (Section 13). Because it is informational, it is
not part of the verified calling method and is recorded as such in Appendix B.

---

## 8. Grouping and per-group population statistics (per dyad × chromosome)

Calls are grouped by (`parent_id`, `offspring_id`, `chrom`) and cross-tabulated
by class into per-dyad-per-chromosome counts, with every class column
guaranteed present even when zero
(`scripts/STEP_TRC_01_classify_intervals.py:290`). From these counts the
following per-group statistics are computed
(`scripts/STEP_TRC_01_classify_intervals.py:302`):

```
n_intervals_total   = sum over the 7 classes
fraction_classified = (n_total - n_AMBIG - n_LOW_CONFIDENCE) / n_total
co_per_mb           = (n_CO + n_DCO) / (chrom_len_bp / 1e6)
nco_per_mb          = n_NCO       / (chrom_len_bp / 1e6)
co_nco_ratio        = (n_CO + n_DCO) / n_NCO        (NaN if n_NCO == 0)
fraction_in_inversions = mean over the group of [inside_inversion == "yes"]
```

`co_per_mb` and `nco_per_mb` are crossover- and gene-conversion event
densities — the recombination-landscape quantities of interest — normalized by
physical chromosome length; `co_nco_ratio` is their balance, with a zero
denominator mapped to NaN rather than infinity
(`scripts/STEP_TRC_01_classify_intervals.py:309`). `fraction_classified` is the
yield of the calling step (the complement of the ambiguous + low-confidence
fraction), and `fraction_in_inversions` is the share of the group's intervals
that fell in suppressed regions
(`scripts/STEP_TRC_01_classify_intervals.py:315`). Note that double crossovers
are pooled with single crossovers in the crossover rate (each is at least one
exchange event for landscape purposes).

> **Flag.** `docs/METHODOLOGY.md:169` states a genome-wide `chrom == "GENOME"`
> row is also emitted. The aggregation code produces only per-(dyad,chrom)
> rows and **no GENOME row**
> (`scripts/STEP_TRC_01_classify_intervals.py:290`). Documented, not
> implemented.

---

## 9. Cross-dyad per-chromosome aggregation (regime construction)

A second aggregation collapses the per-dyad table over dyads to one row per
chromosome, giving the count of contributing dyads, summed crossover /
double-crossover / non-crossover totals, and the **mean and standard deviation
across dyads** of the per-Mb crossover and non-crossover rates
(`scripts/STEP_TRC_01_classify_intervals.py:340`):

```
n_dyads, total_co, total_dco, total_nco,
mean_co_per_mb, sd_co_per_mb, mean_nco_per_mb, sd_nco_per_mb
```

This is the per-chromosome recombination regime: a central tendency and a
between-individual dispersion of the rate. The standard deviation is the pandas
sample SD and is therefore undefined (blank) for any chromosome with a single
contributing dyad — a property visible in the Results.

---

## 10. Relatedness / inheritance context and the karyotype refinement (unwired)

The entire analysis is conditioned on parent–offspring inheritance: the dyad
keys (`parent_id`, `offspring_id`) and the haplotype backgrounds against which
departures are measured are upstream products of a pedigree-aware HMM
(Section 2.1; `docs/HANDOFF_TO_NGSPEDIGREE.md:31`). Inversion karyotype enters
the calling only through the upstream `inside_inversion` flag; this repository
does **not** read the inversion atlas itself in the calling script.

A more precise inside-inversion rule is documented: supply per-sample inversion
karyotypes so the gene-conversion interpretation fires only when the parent is
a confirmed inversion heterozygote (`docs/SCHEMA.md:129`). The loader for this
table exists and validates its columns
(`scripts/STEP_TRC_01_classify_intervals.py:150`), **but its return value is
assigned to `_` and discarded** — it is never passed to the classifier
(`scripts/STEP_TRC_01_classify_intervals.py:386`). The karyotype-aware
refinement is therefore **method only / not wired**, and the executed
inside-inversion rule uses the upstream flag unconditionally. This is recorded
in Appendix B.

---

## 11. Crossover-breakpoint refinement (optional second pass)

The refinement pass localizes crossover breakpoints below interval resolution
and runs **only on `CO` calls** — gene conversions and double crossovers have
no single breakpoint to refine
(`scripts/STEP_TRC_02_traversal_scan.py:103`). Per-site inferred states are
encoded `0 = hapA`, `1 = hapB`, restricted to confidently phased sites, and
windowed around each crossover interval with a ±100 kb pad
(`scripts/STEP_TRC_02_traversal_scan.py:125`,
`scripts/STEP_TRC_02_traversal_scan.py:148`).

Within the local site sequence, a centered window of `W` informative sites
(`window_size`, default 21; `scripts/STEP_TRC_02_traversal_scan.py:93`) is slid
along and the within-window fraction of hapB sites is computed by a cumulative
sum (`scripts/STEP_TRC_02_traversal_scan.py:53`):

```
frac_B[i] = ( cumsum(state)[i + half + 1] - cumsum(state)[i - half] ) / W
grad[i]   = frac_B[i+1] - frac_B[i]
breakpoint index = argmax_i |grad[i]|
```

The breakpoint is reported as the midpoint of the genomic positions of the two
windows bracketing the maximum-gradient step
(`scripts/STEP_TRC_02_traversal_scan.py:62`,
`scripts/STEP_TRC_02_traversal_scan.py:71`). The confidence interval spans the
first and last window whose absolute gradient is within `0.1` of the maximum,
i.e. `|grad| ≥ max(|grad|_max − 0.1, 0)`
(`scripts/STEP_TRC_02_traversal_scan.py:73`). At least `W + 1` sites are
required or the interval is skipped
(`scripts/STEP_TRC_02_traversal_scan.py:48`,
`scripts/STEP_TRC_02_traversal_scan.py:152`). The output carries the refined
breakpoint, its CI, the gradient magnitude, and the window-site count
(`scripts/STEP_TRC_02_traversal_scan.py:165`); it is intended to be joined back
to the calling table on `interval_id` (`docs/METHODOLOGY.md:143`).

---

## 12. Serialization and reproducibility

Three TSVs are written by the calling stage — the per-interval classifications,
the per-dyad-per-chromosome rates, and the per-chromosome aggregate — plus a
provenance file `ngstracts.args.tsv` recording the tool version, the upstream
schema and version, the input paths, the interval count, and **every effective
parameter value** (`scripts/STEP_TRC_01_classify_intervals.py:404`,
`scripts/STEP_TRC_01_classify_intervals.py:421`). The refinement stage writes a
fourth TSV, and writes a well-formed header-only file when there are no
crossovers to refine (`scripts/STEP_TRC_02_traversal_scan.py:104`). All
parameters surface as CLI flags generated directly from the `DEFAULTS`
dictionary, so a run is fully specified by its argument vector and the recorded
`.args` file (`scripts/STEP_TRC_01_classify_intervals.py:371`). The compute is
deterministic (no RNG in the calling or rate paths) and depends only on
`numpy` and `pandas`.

---

## 13. Results

> **Scope of this section.** The headline 226-sample hatchery cohort has **not
> been run**: the upstream stage has not produced real outputs for it
> (`README.md:35`), and no cohort artifacts exist on disk (verified — the
> repository contains no `tract_classifications.tsv`, `dyad_event_rates.tsv`,
> `chrom_summary.tsv`, or `traversal_breakpoints.tsv` outside what is generated
> from the fixture). The only data with real result artifacts is the
> **synthetic conformance fixture** `tests/fixture_basic/`
> (`tests/make_fixture.py:1`). I executed the full chain on it; every number
> below is read verbatim from the generated artifacts. These numbers validate
> the *implementation*; they are **not biological results** and must not be
> read as cohort findings.

### 13.1 Worked example on the synthetic fixture

**Input.** 8 departure intervals across 2 dyads (P1→O1 on `C_gar_LG01`,
40,000,000 bp; P2→O2 on `C_gar_LG28`, 45,000,000 bp)
(`tests/fixture_basic/departure_intervals.tsv`;
`tests/fixture_basic/ref.fa.fai`). The schema gate passed
(`stage3 = 0.1`, `ngstracts = 0.1.0`).

**Per-interval calls** (`tract_classifications.tsv`, 8 rows):

| interval_id | dyad | chrom | start–end | span_bp | class | confidence | review |
| --- | --- | --- | --- | --- | --- | --- | --- |
| DEP_000001 | P1→O1 | C_gar_LG01 | 5,000,000–5,000,500 | 501 | NCO | high | 0 |
| DEP_000002 | P2→O2 | C_gar_LG28 | 16,000,000–16,000,200 | 201 | NCO | high | 0 |
| DEP_000003 | P1→O1 | C_gar_LG01 | 19,900,000–20,100,000 | 200,001 | CO | medium | 0 |
| DEP_000004 | P2→O2 | C_gar_LG28 | 5,000,000–5,100,000 | 100,001 | DCO | medium | 0 |
| DEP_000005 | P2→O2 | C_gar_LG28 | 10,000,000–10,600,000 | 600,001 | MOSAIC_LONG | low | 1 |
| DEP_000006 | P1→O1 | C_gar_LG01 | 2,000,000–2,000,300 | 301 | AMBIG | low | 1 |
| DEP_000007 | P1→O1 | C_gar_LG01 | 7,000,000–7,000,300 | 301 | LOW_CONFIDENCE | low | 1 |
| DEP_000008 | P2→O2 | C_gar_LG28 | 16,500,000–16,600,000 | 100,001 | MOSAIC_SHORT | medium | 1 |

Class breakdown: **NCO 2, CO 1, DCO 1, MOSAIC_SHORT 1, MOSAIC_LONG 1, AMBIG 1,
LOW_CONFIDENCE 1**; **`manual_review_flag` set on 4/8 (50.0%)**. Each call
exercises a distinct branch of Section 5: DEP_000002 and DEP_000008 the
inside-inversion size ladder (short→NCO, mid→MOSAIC_SHORT); DEP_000001 the
short matched-flank NCO; DEP_000003 the flank-switch CO; DEP_000004 the clean
intermediate matched-flank DCO (disc_frac = 68/70 = 0.971 ≥ 0.9); DEP_000005
the long matched-flank MOSAIC_LONG; DEP_000006 the short mismatched-flank AMBIG;
DEP_000007 the `n_sites = 2 < 3` gate.

The position prior column behaves as specified: DEP_000003 sits at the
chromosome midpoint (`f = 0.5`) and carries `prior_log_ratio_co_nco = 2.302585`
= `log 10`, confirming Section 7.

**Per-dyad rates** (`dyad_event_rates.tsv`, 2 rows):

| dyad | chrom | chrom_len | n_co | n_dco | n_nco | co_per_mb | nco_per_mb | co_nco_ratio | frac_classified | frac_in_inv |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| P1→O1 | C_gar_LG01 | 40,000,000 | 1 | 0 | 1 | 0.025 | 0.025 | 1.0 | 0.5 | 0.0 |
| P2→O2 | C_gar_LG28 | 45,000,000 | 0 | 1 | 1 | 0.022222 | 0.022222 | 1.0 | 1.0 | 0.5 |

These match the formulas of Section 8 exactly: e.g. P1→O1 `co_per_mb =
(1+0)/(40e6/1e6) = 0.025`; P2→O2 `co_per_mb = (0+1)/45 = 0.022222`;
`fraction_classified` for P1→O1 = `(4 − 1 − 1)/4 = 0.5` (its AMBIG and
LOW_CONFIDENCE intervals); `fraction_in_inversions` for P2→O2 = `2/4 = 0.5` (its
two inside-inversion intervals).

**Per-chromosome aggregate** (`chrom_summary.tsv`, 2 rows): `C_gar_LG01`
(40 Mb, n_dyads = 1, total_co = 1, total_nco = 1, mean_co_per_mb = 0.025) and
`C_gar_LG28` (45 Mb, n_dyads = 1, total_dco = 1, total_nco = 1, mean_co_per_mb =
0.022222). As predicted in Section 9, both `sd_*` columns are **blank** because
each chromosome has only one contributing dyad (sample SD of one value is
undefined).

**Breakpoint refinement** (`traversal_breakpoints.tsv`, 1 row): the single CO
(DEP_000003) was refined; 0 skipped. From 401 windowed sites the refined
breakpoint is **19,949,500 bp** with CI **[19,810,000, 20,190,000]** and
gradient magnitude **0.047619** (= 1/21, the per-step fraction change of a clean
single-site-resolution transition spread across the 21-site window). The
synthetic ground-truth crossover is at 19,950,000 bp
(`tests/make_fixture.py:133`), so the localization error is **500 bp**, well
inside the smoke test's 50 kb tolerance (`tests/smoke_test.py:100`). The wide CI
is the honest behaviour of the `max − 0.1` rule on a fixture whose gradient is
uniform across the transition region.

### 13.2 What was NOT produced on this data

- **No cohort result.** The 226-sample hatchery cohort
  (`lanta/run_on_lanta.sh:9`) was not run; upstream Stage 3 outputs for it do
  not exist (`README.md:35`). No per-interval, per-dyad, per-chromosome, or
  breakpoint artifacts exist for any real cohort.
- **No genome-wide (`GENOME`) summary row** was produced for any data — the
  feature is documented but unimplemented (Section 8;
  `scripts/STEP_TRC_01_classify_intervals.py:290`).
- **No karyotype-aware NCO refinement** was applied — the loader output is
  discarded (Section 10; `scripts/STEP_TRC_01_classify_intervals.py:386`). The
  fixture run used no `--parent-karyotypes` table.
- **No inversion-aware prior was applied to any call** — the position prior is
  informational and was emitted but never used (Section 7).
- **No between-dyad dispersion** could be estimated: with one dyad per
  chromosome the SD columns are empty.
- **No multi-CO refinement statistics**: only one interval (DEP_000003) was a
  crossover, so the refinement output contains a single row.

---

## Appendix A — Default parameters

All defaults are defined in the `DEFAULTS` dictionary of the calling script and
surfaced as CLI flags (`scripts/STEP_TRC_01_classify_intervals.py:371`); the
window size is a flag on the refinement script.

| Parameter | Default | Role | Source (path:line) |
| --- | --- | --- | --- |
| `min_n_sites_per_interval` | 3 | confidence/coverage gate (rule 0) | `scripts/STEP_TRC_01_classify_intervals.py:34` |
| `short_tract_max_bp` | 50,000 | NCO upper bound / mosaic lower bound | `scripts/STEP_TRC_01_classify_intervals.py:35` |
| `mosaic_short_max_bp` | 200,000 | MOSAIC_SHORT / DCO upper bound | `scripts/STEP_TRC_01_classify_intervals.py:36` |
| `co_breakpoint_max_resolution` | 1,000,000 | CO confidence downgrade threshold | `scripts/STEP_TRC_01_classify_intervals.py:37` |
| `implausible_dco_min_bp` | 5,000,000 | DCO manual-review trigger | `scripts/STEP_TRC_01_classify_intervals.py:38` |
| `implausible_nco_min_bp` | 100,000 | NCO manual-review trigger | `scripts/STEP_TRC_01_classify_intervals.py:39` |
| `discordant_frac_strong` | 0.9 | high-confidence / DCO discordance cut | `scripts/STEP_TRC_01_classify_intervals.py:40` |
| `n_sites_high_conf_min` | 5 | high-confidence site-count floor | `scripts/STEP_TRC_01_classify_intervals.py:41` |
| `chrom_end_proximity_bp` | 2,000,000 | "near chromosome end" (declared; not used in calling) | `scripts/STEP_TRC_01_classify_intervals.py:42` |
| `prior_max_dco_rate` | 0.01 | prior peak at chrom midpoint (informational) | `scripts/STEP_TRC_01_classify_intervals.py:43` |
| `prior_gc_constant` | 0.001 | prior GC baseline (informational) | `scripts/STEP_TRC_01_classify_intervals.py:44` |
| `window_size` (TRC_02) | 21 | sliding-window width in informative sites | `scripts/STEP_TRC_02_traversal_scan.py:93` |

(Note: `chrom_end_proximity_bp` is a declared default with no consumer in the
calling code; it is documented in `docs/METHODOLOGY.md:218` as a "near
chromosome end" scale but no rule references it.)

---

## Appendix B — Statistic provenance and run status

"Run on this data" = executed on the synthetic fixture `tests/fixture_basic/`
with a real artifact on disk (Section 13). No statistic has been run on the
226-sample cohort.

| Statistic / stage | In-repo | External | Run on this data |
| --- | --- | --- | --- |
| Genotype-likelihood encoding & hard calling (thresholds 0.90/0.95) | no | yes (upstream Stage 3; confirm) | no |
| Informative-site selection | no | yes (upstream; confirm) | no |
| Two-state Viterbi HMM haplotype-background inference | no | yes (upstream; confirm) | no |
| Departure detection (run ≥3 discordant sites) | no | yes (upstream; confirm) | no |
| Inversion tagging / `inside_inversion` | no | yes (upstream + atlas; confirm) | no |
| Per-interval upstream `confidence` | no | yes (upstream; confirm) | no |
| Schema-version gate | yes (`STEP_TRC_01:64`) | no | yes |
| Input integrity constraints (span, counts, enums, uniqueness) | yes (`STEP_TRC_01:112`) | no | yes |
| Discordance fraction `disc_frac` | yes (`STEP_TRC_01:187`) | no | yes |
| Flank comparison primitives (`flank_same`, `flanks_real`) | yes (`STEP_TRC_01:189`) | no | yes |
| Decision-tree classification (7 classes) | yes (`STEP_TRC_01:167`) | no | yes |
| Per-call confidence assignment | yes (`STEP_TRC_01:203`) | no | yes |
| Manual-review flagging | yes (`STEP_TRC_01:259`) | no | yes |
| Position prior `prior_log_ratio_co_nco` (informational; not in calling) | yes (`STEP_TRC_01:266`) | no | yes (emitted, unused) |
| Per-dyad rates (`co_per_mb`, `nco_per_mb`, `co_nco_ratio`, fractions) | yes (`STEP_TRC_01:302`) | no | yes |
| Per-chromosome aggregate (means, SDs) | yes (`STEP_TRC_01:340`) | no | yes |
| Genome-wide `GENOME` summary row | documented only (`METHODOLOGY.md:169`) | no | no (not implemented) |
| Karyotype-aware NCO refinement | loader only; **return discarded** (`STEP_TRC_01:386`) | no | no (not wired) |
| CO breakpoint refinement (max-gradient, 21-site window) | yes (`STEP_TRC_02:37`) | no | yes |

---

*End of draft — Methods + Results, ngsTracts analysis chain. All numbers in
Section 13 derive from the synthetic fixture `tests/fixture_basic/`; the
226-sample cohort has not been run.*

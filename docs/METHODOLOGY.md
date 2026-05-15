# ngsTracts methodology specification (v0.1)

Downstream of `ngsPedigree` Stage 3. Classifies departure intervals from
expected single-parent haplotype background as **non-crossover (NCO,
gene conversion)**, **crossover (CO)**, **double crossover (DCO)**, or
**ambiguous**. Per-parent–offspring dyad, per chromosome.

This is the analytical layer that turns Stage 3's "site disagrees with
inherited background" events into biological calls. ngsPedigree owns the
sites; ngsTracts owns what those sites mean.

---

## 1. Scope and non-goals

**Scope.**
1. Classify each Stage 3 `departure_interval` (a run of sites disagreeing
   with the inherited haplotype background) into one of: `NCO`, `CO`,
   `DCO`, `MOSAIC_SHORT`, `MOSAIC_LONG`, `AMBIG`, `LOW_CONFIDENCE`.
2. Per-dyad, per-chromosome event-rate summary (events/Mb, CO:NCO ratio).
3. Recombination-aware priors: the prior on each class depends on
   position relative to known inversions, recombination map (if available),
   and chromosome-end proximity.
4. Output a per-call confidence so downstream filtering is principled.

**Non-goals.**
- ngsTracts does NOT re-do haplotype assignment. That's Stage 3.
- ngsTracts does NOT re-call genotypes. Stage 3 emits hard calls per site;
  ngsTracts trusts them (gated by Stage 3's per-site confidence flag).
- ngsTracts does NOT estimate genome-wide recombination rate from scratch.
  It uses Stage 3 background + inversion atlas as input priors.
- ngsTracts is NOT a phasing tool. Stage 3 has already established which
  parental haplotype block each region descends from.

---

## 2. Input contract (from ngsPedigree Stage 3)

ngsPedigree Stage 3 (`STEP_PED_03_inheritance_map.py`) must emit the four
files described in `SCHEMA.md`. Summary:

- `chromosome_background_intervals.tsv` — per-dyad, per-chrom inferred
  background blocks.
- `departure_intervals.tsv` — the PRIMARY input to ngsTracts.
- `per_site_haplotype_calls.tsv` — optional; only read by STEP_TRC_02.
- `stage3.args.tsv` — provenance + schema version.

The inversion atlas (`inversion_atlas.tsv`) is supplied externally to
both Stage 3 and ngsTracts.

See `SCHEMA.md` for column-level details and constraints.

---

## 3. Decision logic — interval classifier (`STEP_TRC_01`)

For each interval in `departure_intervals.tsv`, apply this decision tree.
The first matching rule wins.

### 3.1 Pre-filter: confidence gate

If `confidence == low` or `n_sites < 3`: emit `LOW_CONFIDENCE`, skip.

### 3.2 Inside-inversion check (priority 1)

If `inside_inversion == yes`:
- Recombination is **suppressed** in inversion heterozygotes (loop-pairing
  doesn't form chiasmata in standard models).
- A short departure inside the inversion almost certainly = gene conversion
  (non-allelic homologous recombination resolved without crossover).
- `span_bp < 50_000` → **`NCO`** (confidence: `high` if `n_sites >= 5`,
  `medium` otherwise).
- `50_000 <= span_bp < 200_000` → **`MOSAIC_SHORT`** (rare; gene
  conversion of unusual size; manual review flag set).
- `span_bp >= 200_000` → **`MOSAIC_LONG`** (suspicious; possible
  recurrent NCO tract OR misclassified background; manual review).

### 3.3 Short-interval rule (priority 2)

If `inside_inversion == no` AND `span_bp < 50_000`:
- **`NCO`** (gene conversion, classical short tract).
- Confidence: `high` if `n_discordant / n_sites >= 0.9` AND `n_sites >= 5`,
  else `medium`.
- Sanity check: flanking states must be the same. If not, downgrade
  to `AMBIG`.

### 3.4 Long-interval-with-flanking-switch rule (priority 3)

If `flanking_left_state != flanking_right_state` (and both are real states,
not `boundary`):
- Single crossover within the interval. **`CO`**.
- Best-estimate breakpoint = midpoint of `[start, end]` (interval is the
  resolution limit; finer resolution requires STEP_TRC_02).
- Confidence: `high` if `span_bp <= 200_000`, `medium` if `200_000 <
  span_bp <= 1_000_000`, `low` otherwise (interval too coarse to localize).

### 3.5 Long-interval-with-matching-flanks rule (priority 4)

If `flanking_left_state == flanking_right_state` AND `span_bp >= 50_000`:
- The departure "comes back" to the original haplotype = double crossover.
- **`DCO`** if `50_000 <= span_bp <= 200_000` AND
  `n_discordant / n_sites >= 0.9`.
- **`MOSAIC_LONG`** otherwise.
- DCO calls always emitted with confidence `medium` (intrinsically rare;
  downstream validation recommended).

### 3.6 Default

Anything unmatched: `AMBIG` with `n_sites`, `span_bp`, and flank states
preserved for manual review.

### 3.7 Position-dependent priors (informational only)

Emitted as a separate column `prior_log_ratio_co_nco`. NOT used to change
the classification (the decision tree above is deterministic on the data).

Logistic position prior:
```
prior_dco(pos, chrom_len) = max_dco_rate * 4 * (pos/chrom_len) * (1 - pos/chrom_len)
```
with `max_dco_rate = 0.01` (at chromosome midpoint), `prior_gc = 0.001`
constant. Used as a confidence modifier and to flag implausible DCOs
at chrom ends.

---

## 4. Traversal scan — `STEP_TRC_02` (refinement layer)

Optional second pass. Only when high-precision crossover breakpoint
localization is required (e.g. for the manuscript's CO landscape figure).

Per (parent, offspring, chrom):
1. Read `per_site_haplotype_calls.tsv` (Stage 3, with `--emit-per-site`).
2. Walk sites in position order.
3. At each site, compute a 21-site sliding concordance with hapA vs hapB.
4. Crossover breakpoint = position of maximum concordance gradient.
5. Confidence interval = leftmost and rightmost site within
   `max_gradient - 0.1`.

Output: `traversal_breakpoints.tsv` with `interval_id`, refined
`breakpoint_bp`, `ci_left_bp`, `ci_right_bp`, `gradient_magnitude`.

Joined back to `STEP_TRC_01`'s table on `interval_id` to give each CO
call a fine-mapped breakpoint when per-site data is available.

---

## 5. Outputs

### 5.1 `tract_classifications.tsv` (PRIMARY)

One row per interval. Columns:

```
interval_id, parent_id, offspring_id, chrom, start, end, span_bp,
class, confidence,
flanking_left_state, flanking_right_state, departure_state,
n_sites, n_discordant, inside_inversion, distance_to_nearest_inv_bp,
prior_log_ratio_co_nco,
refined_breakpoint_bp, refined_ci_left, refined_ci_right,
manual_review_flag, notes
```

`class`: `NCO` | `CO` | `DCO` | `MOSAIC_SHORT` | `MOSAIC_LONG` | `AMBIG` | `LOW_CONFIDENCE`.
`refined_*` columns populated only if STEP_TRC_02 was run.

### 5.2 `dyad_event_rates.tsv` (SUMMARY)

One row per (parent, offspring, chrom). Genome-wide totals also emitted
with `chrom == "GENOME"`.

```
parent_id, offspring_id, chrom, chrom_len_bp,
n_co, n_dco, n_nco, n_mosaic_short, n_mosaic_long, n_ambig, n_low_conf,
co_per_mb, nco_per_mb, co_nco_ratio,
n_intervals_total, fraction_classified, fraction_in_inversions
```

### 5.3 `chrom_summary.tsv`

Per-chromosome aggregated across all dyads:

```
chrom, chrom_len_bp, n_dyads,
total_co, total_dco, total_nco,
mean_co_per_mb, sd_co_per_mb, mean_nco_per_mb, sd_nco_per_mb
```

---

## 6. Quality flags

`manual_review_flag` is set to `1` (otherwise `0`) when ANY of:
- `class == DCO` AND `span_bp > 5_000_000` (implausibly long DCO)
- `class == DCO` AND `inside_inversion == yes` (CO suppressed → DCO impossible)
- `class == NCO` AND `span_bp > 100_000` (longer than typical GC tract)
- `class == MOSAIC_LONG` (always)
- `confidence == low`
- `n_sites < 5`

The flag does not affect downstream summary stats by default; filter
explicitly when computing recombination rates:
`class in {CO, NCO} AND confidence != low AND manual_review_flag == 0`.

---

## 7. Parameter table (defaults)

```
min_n_sites_per_interval       = 3
short_tract_max_bp             = 50_000      # NCO upper bound
mosaic_short_max_bp            = 200_000     # MOSAIC_SHORT upper bound
co_breakpoint_max_resolution   = 1_000_000   # CO confidence downgrade threshold
implausible_dco_min_bp         = 5_000_000   # manual review trigger
implausible_nco_min_bp         = 100_000     # manual review trigger
discordant_frac_strong         = 0.9         # "high" confidence threshold
n_sites_high_conf_min          = 5
chrom_end_proximity_bp         = 2_000_000   # "near chromosome end"
prior_max_dco_rate             = 0.01        # at chrom midpoint
prior_gc_constant              = 0.001
sliding_window_n_sites         = 21          # STEP_TRC_02
```

All parameters surface as CLI flags. Defaults applied when not specified.

---

## 8. Provenance and dependencies

- **Input**: ngsPedigree Stage 3 outputs (locked schema, see `SCHEMA.md`).
- **Inversion atlas**: optional, from MS_Inversions supplementary data.
  If absent, ngsTracts runs in degraded mode (no inversion-aware priors,
  `inside_inversion` assumed `no` everywhere).
- **Reference**: chromosome lengths from `${REF}.fai`.
- **Python deps**: `numpy`, `pandas`. No SLURM dependencies — runs in
  minutes on a login node.

---

## 9. Versioning

This is the v0.1 spec. Schema-breaking changes bump the major version and
require ngsPedigree Stage 3 to emit a matching `schema_version` line in
its `stage3.args.tsv`. ngsTracts refuses to run if schema versions don't
match. See `SCHEMA_COMPATIBILITY.md`.

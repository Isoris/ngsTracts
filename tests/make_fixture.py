#!/usr/bin/env python3
"""
Generate a synthetic Stage 3 output set for testing ngsTracts.

Builds tiny but valid TSVs covering each classification path:
  - one NCO outside inversion (short, matching flanks)
  - one NCO inside inversion (very short)
  - one CO (long, flank switch)
  - one DCO (medium long, matching flanks, high disc frac)
  - one MOSAIC_LONG (very long, matching flanks)
  - one AMBIG (short but flanks differ — impossible NCO)
  - one LOW_CONFIDENCE (n_sites=2 or confidence=low)
  - one inside-inv MOSAIC_SHORT (50-200kb inside inversion)

Plus a matching stage3.args.tsv with schema_version=0.1.

Writes to ./fixture_basic/ relative to this script's parent.
"""
import os
from pathlib import Path

FIXTURE_DIR = Path(__file__).parent / "fixture_basic"
FIXTURE_DIR.mkdir(parents=True, exist_ok=True)


# --- stage3.args.tsv ------------------------------------------------------
(FIXTURE_DIR / "stage3.args.tsv").write_text(
    "key\tvalue\n"
    "schema_version\t0.1\n"
    "stage3_version\t0.1\n"
    "n_dyads\t2\n"
    "n_chroms\t2\n"
    "inversion_atlas_supplied\tyes\n"
    "inversion_atlas_path\tNA\n"
    "hard_call_threshold_hom\t0.95\n"
    "hard_call_threshold_het\t0.90\n"
    "min_run_length_sites\t3\n"
    "beagle_input\tNA\n"
    "sample_list\tNA\n"
    "emit_per_site\tno\n"
    "n_intervals\t8\n"
    "n_background_blocks\t4\n"
    "datetime\t2026-05-12T08:00:00+07:00\n"
    "host\ttest\n"
)


# --- chromosome_background_intervals.tsv ----------------------------------
bg = [
    # parent_id offspring_id chrom start end background_state n_inf_sites concordant_frac notes
    ("P1", "O1", "C_gar_LG01",        1, 20_000_000, "hapA", 1500, "0.992", "-"),
    ("P1", "O1", "C_gar_LG01", 20_000_001, 30_000_000, "hapB",  800, "0.988", "crossover_event_1"),
    ("P2", "O2", "C_gar_LG28",        1, 30_000_000, "hapA", 2400, "0.985", "-"),
    ("P2", "O2", "C_gar_LG28", 30_000_001, 40_000_000, "hapA", 1100, "0.989", "-"),
]
with open(FIXTURE_DIR / "chromosome_background_intervals.tsv", "w") as f:
    f.write("parent_id\toffspring_id\tchrom\tstart\tend\tbackground_state\tn_inf_sites\tconcordant_frac\tnotes\n")
    for row in bg:
        f.write("\t".join(map(str, row)) + "\n")


# --- departure_intervals.tsv ---------------------------------------------
def row(iid, p, o, ch, s, e, ns, nd, fl, fr, ds, dist, ininv, conf):
    span = e - s + 1
    return (iid, p, o, ch, s, e, ns, nd, span, fl, fr, ds, dist, ininv, conf)

intervals = [
    # 1) NCO outside inversion (short, matching flanks, high disc frac)
    row("DEP_000001", "P1", "O1", "C_gar_LG01",
        5_000_000, 5_000_500,  8, 8,
        "hapA", "hapA", "hapB", "1245000", "no", "high"),

    # 2) NCO inside inversion (very short)
    row("DEP_000002", "P2", "O2", "C_gar_LG28",
        16_000_000, 16_000_200, 6, 6,
        "hapA", "hapA", "hapB", "0", "yes", "high"),

    # 3) CO — long, flank switch
    row("DEP_000003", "P1", "O1", "C_gar_LG01",
        19_900_000, 20_100_000, 80, 60,
        "hapA", "hapB", "hapB", "-", "no", "high"),

    # 4) DCO — 100kb span, matching flanks, high disc
    row("DEP_000004", "P2", "O2", "C_gar_LG28",
        5_000_000, 5_100_000, 70, 68,
        "hapA", "hapA", "hapB", "-", "no", "high"),

    # 5) MOSAIC_LONG — 600kb span, matching flanks
    row("DEP_000005", "P2", "O2", "C_gar_LG28",
        10_000_000, 10_600_000, 200, 180,
        "hapA", "hapA", "hapB", "-", "no", "medium"),

    # 6) AMBIG — short but flanks differ
    row("DEP_000006", "P1", "O1", "C_gar_LG01",
        2_000_000, 2_000_300,  5, 5,
        "hapA", "hapB", "hapB", "-", "no", "high"),

    # 7) LOW_CONFIDENCE — n_sites=2
    row("DEP_000007", "P1", "O1", "C_gar_LG01",
        7_000_000, 7_000_300,  2, 2,
        "hapA", "hapA", "hapB", "-", "no", "high"),

    # 8) MOSAIC_SHORT — 100kb inside inversion
    row("DEP_000008", "P2", "O2", "C_gar_LG28",
        16_500_000, 16_600_000, 60, 58,
        "hapA", "hapA", "hapB", "0", "yes", "high"),
]

with open(FIXTURE_DIR / "departure_intervals.tsv", "w") as f:
    f.write("\t".join([
        "interval_id", "parent_id", "offspring_id", "chrom",
        "start", "end", "n_sites", "n_discordant", "span_bp",
        "flanking_left_state", "flanking_right_state", "departure_state",
        "distance_to_nearest_inv_bp", "inside_inversion", "confidence",
    ]) + "\n")
    for r in intervals:
        f.write("\t".join(map(str, r)) + "\n")


# --- chrom lengths (.fai) -------------------------------------------------
fai_path = FIXTURE_DIR / "ref.fa.fai"
with open(fai_path, "w") as f:
    f.write("C_gar_LG01\t40000000\t0\t60\t61\n")
    f.write("C_gar_LG28\t45000000\t40000061\t60\t61\n")


# --- per_site_haplotype_calls.tsv (small, for STEP_TRC_02) ---------------
# Just for the CO interval (DEP_000003): hapA before 19.95M, hapB after
ps_rows = []
import random
random.seed(42)
for pos in range(19_800_000, 20_200_001, 1000):
    state = "hapA" if pos < 19_950_000 else "hapB"
    ps_rows.append(("P1", "O1", "C_gar_LG01", pos, "0/1", "0/0", state, "0",
                   "DEP_000003" if 19_900_000 <= pos <= 20_100_000 else "-"))

with open(FIXTURE_DIR / "per_site_haplotype_calls.tsv", "w") as f:
    f.write("\t".join([
        "parent_id", "offspring_id", "chrom", "pos",
        "parent_gt", "offspring_gt", "inferred_state", "departure", "interval_id",
    ]) + "\n")
    for r in ps_rows:
        f.write("\t".join(map(str, r)) + "\n")


print(f"Wrote fixture to {FIXTURE_DIR}")
print(f"  stage3.args.tsv")
print(f"  chromosome_background_intervals.tsv ({len(bg)} blocks)")
print(f"  departure_intervals.tsv             ({len(intervals)} intervals)")
print(f"  per_site_haplotype_calls.tsv        ({len(ps_rows)} sites)")
print(f"  ref.fa.fai                          (2 chroms)")

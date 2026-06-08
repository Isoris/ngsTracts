#!/usr/bin/env python3
"""
Smoke test for STEP_TRC_06 (LRR CO-suppression ranking + population support).

Synthesizes a 5-dyad cohort with informative per-site coverage and three
candidate LRRs:
  LRR_A  suppressed: 5 opportunity dyads, none with an interior CO  -> no_CO_inside
  LRR_B  not suppressed: 3 of 5 dyads have an interior CO            -> many_CO_inside
  LRR_C  boundary: interior clean but an edge CO present             -> CO_at_boundary_only
Asserts patterns, opportunity counts, and that the population support score
(Wilson lower bound) ranks A and C above B.
"""
import subprocess
import sys
from math import sqrt
from pathlib import Path

import pandas as pd


HERE = Path(__file__).parent
ROOT = HERE.parent
OUTDIR = HERE / "tmp_smoke_trc06"
CHROM = "C_gar_LG01"
DYADS = [(f"P{i}", f"O{i}") for i in range(1, 6)]  # 5 dyads


def wilson_lower(k, n, z=1.96):
    phat = k / n
    denom = 1 + z * z / n
    center = phat + z * z / (2 * n)
    margin = z * sqrt(phat * (1 - phat) / n + z * z / (4 * n * n))
    return max(0.0, (center - margin) / denom)


def build_inputs(d: Path):
    # LRRs (BED, 0-based half-open)
    (d / "lrrs.bed").write_text(
        f"{CHROM}\t999999\t2000000\tLRR_A\n"      # interior 1.1M-1.9M
        f"{CHROM}\t4999999\t6000000\tLRR_B\n"      # interior 5.1M-5.9M
        f"{CHROM}\t7999999\t9000000\tLRR_C\n")     # interior 8.1M-8.9M

    # per-site informative coverage: every dyad, every 20kb across 1M-9M
    ps_cols = ["parent_id", "offspring_id", "chrom", "pos",
               "parent_gt", "offspring_gt", "inferred_state", "departure", "interval_id"]
    ps_rows = []
    for (pa, of) in DYADS:
        for pos in range(1_000_000, 9_000_001, 20_000):
            ps_rows.append((pa, of, CHROM, pos, "0/1", "0/0", "hapA", "0", "-"))
    pd.DataFrame(ps_rows, columns=ps_cols).to_csv(d / "per_site.tsv", sep="\t", index=False)

    # tract_classifications: CO/NCO events
    tc_cols = ["interval_id", "parent_id", "offspring_id", "chrom", "start", "end", "class"]
    tc_rows = []
    n = 0

    def add(pa, of, s, e, cls):
        nonlocal n
        n += 1
        tc_rows.append((f"DEP_{n:06d}", pa, of, CHROM, s, e, cls))

    # LRR_B: 3 dyads with an interior CO near 5.5M
    for (pa, of) in DYADS[:3]:
        add(pa, of, 5_490_000, 5_510_000, "CO")
    # LRR_A: one interior NCO (gene conversion is allowed inside a block), no CO
    add("P1", "O1", 1_500_000, 1_500_400, "NCO")
    # LRR_C: an edge CO (~8.05M, inside the 0.1*1M = 100kb edge), no interior CO
    add("P2", "O2", 8_045_000, 8_055_000, "CO")
    pd.DataFrame(tc_rows, columns=tc_cols).to_csv(d / "tracts.tsv", sep="\t", index=False)


def main():
    errors = 0
    OUTDIR.mkdir(parents=True, exist_ok=True)
    build_inputs(OUTDIR)

    print("[smoke-trc06] running STEP_TRC_06...")
    r = subprocess.run([
        sys.executable, str(ROOT / "scripts" / "STEP_TRC_06_lrr_suppression.py"),
        "--lrr-bed", str(OUTDIR / "lrrs.bed"),
        "--tract-classifications", str(OUTDIR / "tracts.tsv"),
        "--per-site-file", str(OUTDIR / "per_site.tsv"),
        "--outdir", str(OUTDIR),
    ], check=False, capture_output=True, text=True)
    print(r.stdout); print(r.stderr)
    if r.returncode != 0:
        sys.exit("STEP_TRC_06 failed")

    df = pd.read_csv(OUTDIR / "lrr_co_suppression.tsv", sep="\t").set_index("lrr_id")

    def chk(cond, msg):
        nonlocal errors
        print(f"  [{'OK  ' if cond else 'FAIL'}] {msg}")
        errors += not cond

    A, B, C = df.loc["LRR_A"], df.loc["LRR_B"], df.loc["LRR_C"]

    chk(A["pattern"] == "no_CO_inside", f"LRR_A pattern = {A['pattern']} (want no_CO_inside)")
    chk(A["n_dyads_with_opportunity"] == 5, f"LRR_A opportunity dyads = {A['n_dyads_with_opportunity']} (want 5)")
    chk(A["n_interior_CO"] == 0 and A["n_interior_NCO"] == 1,
        f"LRR_A interior CO/NCO = {A['n_interior_CO']}/{A['n_interior_NCO']} (want 0/1)")
    chk(abs(A["population_support_score"] - wilson_lower(5, 5)) < 1e-6,
        f"LRR_A score = {A['population_support_score']:.4f} (want {wilson_lower(5,5):.4f})")

    chk(B["pattern"] == "many_CO_inside", f"LRR_B pattern = {B['pattern']} (want many_CO_inside)")
    chk(B["n_interior_CO"] == 3, f"LRR_B interior CO = {B['n_interior_CO']} (want 3)")
    chk(abs(B["suppression_consistency"] - 0.4) < 1e-9,
        f"LRR_B consistency = {B['suppression_consistency']} (want 0.4)")

    chk(C["pattern"] == "CO_at_boundary_only", f"LRR_C pattern = {C['pattern']} (want CO_at_boundary_only)")
    chk(C["n_edge_CO"] == 1 and C["n_interior_CO"] == 0,
        f"LRR_C edge/interior CO = {C['n_edge_CO']}/{C['n_interior_CO']} (want 1/0)")

    chk(A["population_support_score"] > B["population_support_score"]
        and C["population_support_score"] > B["population_support_score"],
        "population support ranks suppressed A & boundary C above recombining B")

    if errors:
        print(f"\n[smoke-trc06] {errors} failures"); sys.exit(1)
    print("\n[smoke-trc06] ALL PASS")


if __name__ == "__main__":
    main()

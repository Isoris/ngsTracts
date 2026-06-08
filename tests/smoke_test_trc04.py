#!/usr/bin/env python3
"""
Smoke test for STEP_TRC_04 (event-statistics report).

Builds the fixture, runs STEP_TRC_01 + STEP_TRC_02, then STEP_TRC_04 and
asserts the computed statistics and the requires/model status flags.
"""
import subprocess
import sys
from pathlib import Path

import pandas as pd


HERE = Path(__file__).parent
ROOT = HERE.parent
FIXTURE = HERE / "fixture_basic"
OUTDIR = HERE / "tmp_smoke_trc04"
TC = OUTDIR / "trc01"
STATS = OUTDIR / "stats"


def run(cmd):
    r = subprocess.run(cmd, check=False, capture_output=True, text=True)
    print(r.stdout); print(r.stderr)
    return r.returncode


def main():
    errors = 0
    OUTDIR.mkdir(parents=True, exist_ok=True)

    if run([sys.executable, str(HERE / "make_fixture.py")]) != 0:
        sys.exit("fixture build failed")
    if run([
        sys.executable, str(ROOT / "scripts" / "STEP_TRC_01_classify_intervals.py"),
        "--stage3-dir", str(FIXTURE), "--fai", str(FIXTURE / "ref.fa.fai"),
        "--outdir", str(TC),
    ]) != 0:
        sys.exit("STEP_TRC_01 failed")
    run([
        sys.executable, str(ROOT / "scripts" / "STEP_TRC_02_traversal_scan.py"),
        "--tract-classifications", str(TC / "tract_classifications.tsv"),
        "--per-site-file", str(FIXTURE / "per_site_haplotype_calls.tsv"),
        "--outdir", str(TC),
    ])

    print("[smoke-trc04] running STEP_TRC_04...")
    if run([
        sys.executable, str(ROOT / "scripts" / "STEP_TRC_04_event_statistics.py"),
        "--tract-classifications", str(TC / "tract_classifications.tsv"),
        "--fai", str(FIXTURE / "ref.fa.fai"),
        "--traversal-breakpoints", str(TC / "traversal_breakpoints.tsv"),
        "--switch-error-rate", "0.001",
        "--outdir", str(STATS),
    ]) != 0:
        sys.exit("STEP_TRC_04 failed")

    st = pd.read_csv(STATS / "event_statistics.tsv", sep="\t")

    def val(stat, stratum="ALL"):
        r = st[(st["statistic"] == stat) & (st["stratum"] == stratum)]
        return None if r.empty else float(r["value"].iloc[0])

    def status(stat, stratum="ALL"):
        r = st[(st["statistic"] == stat) & (st["stratum"] == stratum)]
        return None if r.empty else r["status"].iloc[0]

    checks = [
        ("n_CO", "ALL", 1.0), ("n_NCO", "ALL", 2.0), ("n_DCO", "ALL", 1.0),
        ("NCO_CO_ratio", "ALL", 2.0), ("DCO_CO_ratio", "ALL", 1.0),
        ("gene_flux_index_count", "ALL", 3.0),
        ("NCO_length_median_bp", "ALL", 351.0),
        ("NCO_length_median_bp", "inside_inversion=yes", 201.0),
        ("NCO_length_median_bp", "inside_inversion=no", 501.0),
    ]
    for stat, strat, want in checks:
        got = val(stat, strat)
        ok = got == want
        print(f"  [{'OK  ' if ok else 'FAIL'}] {stat}[{strat}] = {got} (want {want})")
        errors += not ok

    # status flags: external/optional must be flagged, not fabricated
    status_checks = [
        ("GC_bias_at_NCO", "ALL", "requires:external(popstats allele/ancestral)"),
        ("rho_inside_outside_ratio", "ALL", "requires:external(LD map: pyrho/LDhat)"),
        ("NCO_rate_ratio_inside_outside", "ALL", "requires:inversion-atlas"),
        ("CO_count_inside_Het", "karyotype=Het", "requires:parent-karyotypes+inversion-atlas"),
        ("event_rate_by_parent_sex", "ALL", "requires:sample-sex"),
        ("expected_false_NCO_upper_bound", "ALL", "model:n_sites*eps^2"),
    ]
    for stat, strat, want in status_checks:
        got = status(stat, strat)
        ok = got == want
        print(f"  [{'OK  ' if ok else 'FAIL'}] status {stat}[{strat}] = {got}")
        errors += not ok

    # No fabricated numbers: every requires:/external row must have empty value
    bad = st[st["status"].str.startswith("requires:") & st["value"].notna()]
    if len(bad):
        print(f"  [FAIL] {len(bad)} requires: rows carry a value (should be NA)"); errors += 1
    else:
        print("  [OK  ] all requires: rows are NA (no fabricated numbers)")

    if errors:
        print(f"\n[smoke-trc04] {errors} failures"); sys.exit(1)
    print("\n[smoke-trc04] ALL PASS")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Smoke test for ngsTracts. Builds the fixture, runs STEP_TRC_01, and asserts
each test interval gets its expected class.
"""
import subprocess
import sys
from pathlib import Path

import pandas as pd


HERE = Path(__file__).parent
ROOT = HERE.parent
FIXTURE = HERE / "fixture_basic"
OUTDIR = HERE / "tmp_smoke_out"

EXPECTED = {
    "DEP_000001": "NCO",
    "DEP_000002": "NCO",            # inside inversion, short
    "DEP_000003": "CO",
    "DEP_000004": "DCO",
    "DEP_000005": "MOSAIC_LONG",
    "DEP_000006": "AMBIG",
    "DEP_000007": "LOW_CONFIDENCE",
    "DEP_000008": "MOSAIC_SHORT",   # inside inversion, 50-200kb
}


def main():
    # Build fixture
    print("[smoke] building fixture...")
    r = subprocess.run([sys.executable, str(HERE / "make_fixture.py")],
                       check=False, capture_output=True, text=True)
    print(r.stdout)
    if r.returncode != 0:
        print(r.stderr, file=sys.stderr)
        sys.exit(1)

    # Run STEP_TRC_01
    OUTDIR.mkdir(parents=True, exist_ok=True)
    print("[smoke] running STEP_TRC_01...")
    cmd = [
        sys.executable, str(ROOT / "scripts" / "STEP_TRC_01_classify_intervals.py"),
        "--stage3-dir", str(FIXTURE),
        "--fai", str(FIXTURE / "ref.fa.fai"),
        "--outdir", str(OUTDIR),
    ]
    r = subprocess.run(cmd, check=False, capture_output=True, text=True)
    print(r.stdout)
    print(r.stderr)
    if r.returncode != 0:
        print(f"[smoke] STEP_TRC_01 failed (exit {r.returncode})", file=sys.stderr)
        sys.exit(1)

    # Verify each interval's class
    tc = pd.read_csv(OUTDIR / "tract_classifications.tsv", sep="\t")
    got = dict(zip(tc["interval_id"], tc["class"]))

    errors = 0
    for iid, expected_cls in EXPECTED.items():
        actual = got.get(iid)
        ok = actual == expected_cls
        marker = "OK  " if ok else "FAIL"
        print(f"  [{marker}] {iid}: expected={expected_cls}  got={actual}")
        if not ok:
            errors += 1

    # Verify dyad and chrom summaries exist and have rows
    dyad = pd.read_csv(OUTDIR / "dyad_event_rates.tsv", sep="\t")
    chrom = pd.read_csv(OUTDIR / "chrom_summary.tsv", sep="\t")
    print(f"\n[smoke] dyad_event_rates.tsv: {len(dyad)} rows")
    print(f"[smoke] chrom_summary.tsv:    {len(chrom)} rows")
    if len(dyad) == 0:
        print("[smoke] dyad summary is empty!"); errors += 1
    if len(chrom) == 0:
        print("[smoke] chrom summary is empty!"); errors += 1

    # Run STEP_TRC_02 on the CO call
    print("\n[smoke] running STEP_TRC_02...")
    cmd = [
        sys.executable, str(ROOT / "scripts" / "STEP_TRC_02_traversal_scan.py"),
        "--tract-classifications", str(OUTDIR / "tract_classifications.tsv"),
        "--per-site-file", str(FIXTURE / "per_site_haplotype_calls.tsv"),
        "--outdir", str(OUTDIR),
    ]
    r = subprocess.run(cmd, check=False, capture_output=True, text=True)
    print(r.stdout)
    print(r.stderr)
    if r.returncode != 0:
        print("[smoke] STEP_TRC_02 failed", file=sys.stderr)
        errors += 1
    else:
        tb = pd.read_csv(OUTDIR / "traversal_breakpoints.tsv", sep="\t")
        print(f"[smoke] traversal_breakpoints.tsv: {len(tb)} rows")
        if len(tb) >= 1:
            # The synthetic CO breakpoint is at ~19_950_000
            bp = int(tb.iloc[0]["refined_breakpoint_bp"])
            print(f"[smoke] refined breakpoint = {bp:,} (truth ~19,950,000)")
            if abs(bp - 19_950_000) > 50_000:
                print(f"[smoke] FAIL: breakpoint off by {abs(bp - 19_950_000):,} bp")
                errors += 1
        else:
            print("[smoke] FAIL: STEP_TRC_02 produced no rows")
            errors += 1

    if errors:
        print(f"\n[smoke] {errors} failures")
        sys.exit(1)
    print("\n[smoke] ALL PASS")


if __name__ == "__main__":
    main()

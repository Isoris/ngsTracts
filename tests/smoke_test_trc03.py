#!/usr/bin/env python3
"""
Smoke test for STEP_TRC_03 (exchange/gene-flux segment mask).

Builds the fixture, runs STEP_TRC_01 to produce tract_classifications.tsv,
then runs STEP_TRC_03 in two modes and asserts the emitted segments.
"""
import subprocess
import sys
from pathlib import Path

import pandas as pd


HERE = Path(__file__).parent
ROOT = HERE.parent
FIXTURE = HERE / "fixture_basic"
OUTDIR = HERE / "tmp_smoke_trc03"
TC_DIR = OUTDIR / "trc01"
MASK_ALL = OUTDIR / "mask_all"
MASK_INV = OUTDIR / "mask_inv"

BED_COLS = ["chrom", "start", "end", "name", "score", "strand"]

# All exchange classes from the fixture (NCO/DCO/MOSAIC_*), interval -> class.
EXPECTED_ALL = {
    "DEP_000001": "NCO",
    "DEP_000002": "NCO",
    "DEP_000004": "DCO",
    "DEP_000005": "MOSAIC_LONG",
    "DEP_000008": "MOSAIC_SHORT",
}
# inside_inversion in {yes, partial} only.
EXPECTED_INV = {"DEP_000002", "DEP_000008"}


def run(cmd):
    r = subprocess.run(cmd, check=False, capture_output=True, text=True)
    print(r.stdout); print(r.stderr)
    return r.returncode


def iid_of(name):  # "CLASS|DEP_xxxx|P>O" -> "DEP_xxxx"
    return name.split("|")[1]


def cls_of(name):
    return name.split("|")[0]


def main():
    errors = 0
    OUTDIR.mkdir(parents=True, exist_ok=True)

    print("[smoke-trc03] building fixture...")
    if run([sys.executable, str(HERE / "make_fixture.py")]) != 0:
        sys.exit("fixture build failed")

    print("[smoke-trc03] running STEP_TRC_01...")
    if run([
        sys.executable, str(ROOT / "scripts" / "STEP_TRC_01_classify_intervals.py"),
        "--stage3-dir", str(FIXTURE),
        "--fai", str(FIXTURE / "ref.fa.fai"),
        "--outdir", str(TC_DIR),
    ]) != 0:
        sys.exit("STEP_TRC_01 failed")
    tc_path = TC_DIR / "tract_classifications.tsv"

    # ---- (a) default: all exchange classes ----
    print("[smoke-trc03] running STEP_TRC_03 (all exchange tracts)...")
    if run([
        sys.executable, str(ROOT / "scripts" / "STEP_TRC_03_exchange_mask.py"),
        "--tract-classifications", str(tc_path),
        "--outdir", str(MASK_ALL),
    ]) != 0:
        sys.exit("STEP_TRC_03 (all) failed")

    bed = pd.read_csv(MASK_ALL / "exchange_segments.bed", sep="\t", header=None, names=BED_COLS)
    got = {iid_of(n): cls_of(n) for n in bed["name"]}
    for iid, cls in EXPECTED_ALL.items():
        ok = got.get(iid) == cls
        print(f"  [{'OK  ' if ok else 'FAIL'}] {iid}: expected={cls} got={got.get(iid)}")
        errors += not ok
    if len(bed) != len(EXPECTED_ALL):
        print(f"  [FAIL] expected {len(EXPECTED_ALL)} segments, got {len(bed)}"); errors += 1
    # CO / AMBIG / LOW_CONFIDENCE must NOT appear
    for forbidden in ("DEP_000003", "DEP_000006", "DEP_000007"):
        if forbidden in got:
            print(f"  [FAIL] non-exchange interval leaked into mask: {forbidden}"); errors += 1

    # BED coords are 0-based half-open: DEP_000001 is 1-based 5000000..5000500
    row1 = bed[bed["name"].str.contains("DEP_000001")].iloc[0]
    if not (row1["start"] == 4999999 and row1["end"] == 5000500):
        print(f"  [FAIL] DEP_000001 BED coords {row1['start']}-{row1['end']} != 4999999-5000500")
        errors += 1

    summ = pd.read_csv(MASK_ALL / "exchange_mask_summary.tsv", sep="\t")
    n_all = int(summ[(summ["group"] == "ALL") & (summ["key"] == "exchange_tracts")]["n_tracts"].iloc[0])
    if n_all != 5:
        print(f"  [FAIL] summary exchange_tracts={n_all} != 5"); errors += 1

    # ---- (b) inside-only with breakpoint pad ----
    print("[smoke-trc03] running STEP_TRC_03 (inside-only, pad 25kb)...")
    if run([
        sys.executable, str(ROOT / "scripts" / "STEP_TRC_03_exchange_mask.py"),
        "--tract-classifications", str(tc_path),
        "--outdir", str(MASK_INV), "--inside-only", "--pad", "25000",
    ]) != 0:
        sys.exit("STEP_TRC_03 (inside-only) failed")

    bed_inv = pd.read_csv(MASK_INV / "exchange_segments.bed", sep="\t", header=None, names=BED_COLS)
    got_inv = {iid_of(n) for n in bed_inv["name"]}
    if got_inv != EXPECTED_INV:
        print(f"  [FAIL] inside-only set {got_inv} != {EXPECTED_INV}"); errors += 1
    else:
        print(f"  [OK  ] inside-only segments = {sorted(got_inv)}")
    # pad applied: DEP_000002 1-based start 16000000 -> 16000000-1-25000 = 15974999
    r2 = bed_inv[bed_inv["name"].str.contains("DEP_000002")].iloc[0]
    if not (r2["start"] == 15974999 and r2["end"] == 16025200):
        print(f"  [FAIL] DEP_000002 padded coords {r2['start']}-{r2['end']} != 15974999-16025200")
        errors += 1
    else:
        print("  [OK  ] breakpoint pad applied correctly")

    # ---- (c) per-sample BEDs (GL-based masking; keyed on offspring carrier) ----
    print("[smoke-trc03] running STEP_TRC_03 (--per-sample)...")
    MASK_PS = OUTDIR / "mask_ps"
    if run([
        sys.executable, str(ROOT / "scripts" / "STEP_TRC_03_exchange_mask.py"),
        "--tract-classifications", str(tc_path),
        "--outdir", str(MASK_PS), "--per-sample",
    ]) != 0:
        sys.exit("STEP_TRC_03 (per-sample) failed")

    by = pd.read_csv(MASK_PS / "exchange_segments.by_sample.tsv", sep="\t")
    # O1 carries DEP_000001 only; O2 carries the four LG28 tracts.
    per_sample_counts = by.groupby("sample_id").size().to_dict()
    if per_sample_counts != {"O1": 1, "O2": 4}:
        print(f"  [FAIL] per-sample counts {per_sample_counts} != {{'O1':1,'O2':4}}"); errors += 1
    else:
        print(f"  [OK  ] per-sample counts = {per_sample_counts}")
    if set(by["role"].unique()) != {"offspring"}:
        print(f"  [FAIL] expected carrier role=offspring, got {set(by['role'].unique())}"); errors += 1
    o1 = pd.read_csv(MASK_PS / "by_sample" / "O1.exchange.bed", sep="\t",
                     header=None, names=BED_COLS)
    o2 = pd.read_csv(MASK_PS / "by_sample" / "O2.exchange.bed", sep="\t",
                     header=None, names=BED_COLS)
    if not (len(o1) == 1 and len(o2) == 4):
        print(f"  [FAIL] O1.bed={len(o1)} (want 1), O2.bed={len(o2)} (want 4)"); errors += 1
    else:
        print("  [OK  ] per-sample BED files written with correct row counts")

    if errors:
        print(f"\n[smoke-trc03] {errors} failures"); sys.exit(1)
    print("\n[smoke-trc03] ALL PASS")


if __name__ == "__main__":
    main()

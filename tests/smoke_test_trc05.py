#!/usr/bin/env python3
"""
Smoke test for STEP_TRC_05 (transmission / meiotic-drive).

Builds a tiny synthetic end-of-polarization adapter table (polarized
transmissions) with one driven inversion and one balanced inversion, runs
STEP_TRC_05, and asserts the transmission counts and binomial drive p-values.
"""
import subprocess
import sys
from pathlib import Path

import pandas as pd


HERE = Path(__file__).parent
ROOT = HERE.parent
OUTDIR = HERE / "tmp_smoke_trc05"


def write_adapter(path: Path):
    rows = []
    # INV_001: 8 informative HET transmissions, 7 INV / 1 REF  -> drive
    arr = ["INV"] * 7 + ["REF"] * 1
    for i, a in enumerate(arr):
        rows.append(("INV_001", "C_gar_LG01", f"PF{i}", f"O{i}", "HET", "father", a, "1"))
    # INV_002: 6 informative HET transmissions, 3 INV / 3 REF  -> balanced
    arr2 = ["INV", "REF"] * 3
    for i, a in enumerate(arr2):
        rows.append(("INV_002", "C_gar_LG28", f"PM{i}", f"Q{i}", "HET", "mother", a, "1"))
    # noise that must be excluded: HOM parents + ambiguous + non-informative
    rows.append(("INV_001", "C_gar_LG01", "PH", "OH", "HOM_INV", "father", "INV", "1"))
    rows.append(("INV_001", "C_gar_LG01", "PA", "OA", "HET", "father", "ambiguous", "1"))
    rows.append(("INV_002", "C_gar_LG28", "PN", "ON", "HET", "mother", "INV", "0"))
    cols = ["inversion_id", "chrom", "parent_id", "offspring_id",
            "parent_karyotype", "parent_role", "transmitted_arrangement", "informative"]
    pd.DataFrame(rows, columns=cols).to_csv(path, sep="\t", index=False)


def main():
    errors = 0
    OUTDIR.mkdir(parents=True, exist_ok=True)
    adapter = OUTDIR / "polarized_transmissions.tsv"
    write_adapter(adapter)

    print("[smoke-trc05] running STEP_TRC_05...")
    r = subprocess.run([
        sys.executable, str(ROOT / "scripts" / "STEP_TRC_05_transmission_drive.py"),
        "--transmissions", str(adapter), "--outdir", str(OUTDIR),
    ], check=False, capture_output=True, text=True)
    print(r.stdout); print(r.stderr)
    if r.returncode != 0:
        sys.exit("STEP_TRC_05 failed")

    td = pd.read_csv(OUTDIR / "transmission_drive.tsv", sep="\t").set_index("inversion_id")

    def chk(cond, msg):
        nonlocal errors
        print(f"  [{'OK  ' if cond else 'FAIL'}] {msg}")
        errors += not cond

    a = td.loc["INV_001"]
    chk(a["n_INV_transmitted"] == 7 and a["n_REF_transmitted"] == 1,
        f"INV_001 counts INV/REF = {a['n_INV_transmitted']}/{a['n_REF_transmitted']} (want 7/1)")
    chk(abs(a["INV_transmission_rate"] - 0.875) < 1e-9,
        f"INV_001 rate = {a['INV_transmission_rate']} (want 0.875)")
    # exact two-sided binomial p for k=7,n=8: 18/256 = 0.0703125
    chk(abs(a["mendelian_drive_pvalue"] - 0.0703125) < 1e-6,
        f"INV_001 drive p = {a['mendelian_drive_pvalue']} (want 0.0703125)")
    chk(abs(a["paternal_INV_transmission_rate"] - 0.875) < 1e-9,
        "INV_001 paternal rate routed from parent_role=father")

    b = td.loc["INV_002"]
    chk(b["n_INV_transmitted"] == 3 and b["n_REF_transmitted"] == 3,
        f"INV_002 counts INV/REF = {b['n_INV_transmitted']}/{b['n_REF_transmitted']} (want 3/3)")
    chk(abs(b["mendelian_drive_pvalue"] - 1.0) < 1e-9,
        f"INV_002 drive p = {b['mendelian_drive_pvalue']} (want 1.0)")
    chk(b["n_informative_het_parent_transmissions"] == 6,
        "INV_002 excluded the non-informative (informative=0) row")

    if errors:
        print(f"\n[smoke-trc05] {errors} failures"); sys.exit(1)
    print("\n[smoke-trc05] ALL PASS")


if __name__ == "__main__":
    main()

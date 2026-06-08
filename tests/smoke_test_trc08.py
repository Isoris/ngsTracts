#!/usr/bin/env python3
"""
Smoke test for STEP_TRC_08 (MicrohapLD_HOMHET).

One LRR, 5 windows, 26 samples:
  HOM group (16): 8 UU + 8 VV -> within HOM each chromosome is pure U or pure V
                  across windows -> near-perfect cross-window LD (high).
  HET group (10): each window's hap1/hap2 shuffled independently -> within HET
                  windows are uncorrelated -> ~no LD (disrupted).
Expect: median LD HOM high, HET very low, interpretation coherent_HOM_disrupted_HET.
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).parent
ROOT = HERE.parent
OUTDIR = HERE / "tmp_smoke_trc08"
CHROM = "C_gar_LG01"


def build(d: Path):
    rng = np.random.default_rng(8)
    mh_rows, kt_rows = [], []
    starts = list(range(1_000_000, 2_300_000, 300_000))  # 5 windows

    hom_uu = [f"hom_uu{i:02d}" for i in range(8)]
    hom_vv = [f"hom_vv{i:02d}" for i in range(8)]
    het = [f"het{i:02d}" for i in range(10)]

    for st in starts:
        en = st + 5_000
        for s in hom_uu:
            mh_rows.append((s, CHROM, st, en, "U", "U"))
        for s in hom_vv:
            mh_rows.append((s, CHROM, st, en, "V", "V"))
        for s in het:
            h1, h2 = rng.choice(["U", "V"], 2)  # independent per window -> mixed
            mh_rows.append((s, CHROM, st, en, h1, h2))

    for s in hom_uu:
        kt_rows.append((s, "LRR_TEST", "0"))   # HOM_REF
    for s in hom_vv:
        kt_rows.append((s, "LRR_TEST", "2"))   # HOM_INV
    for s in het:
        kt_rows.append((s, "LRR_TEST", "1"))   # HET

    pd.DataFrame(mh_rows, columns=["sample", "chrom", "start", "end", "hap1", "hap2"]) \
        .to_csv(d / "microhaps.tsv", sep="\t", index=False)
    pd.DataFrame(kt_rows, columns=["sample", "lrr_id", "karyotype"]) \
        .to_csv(d / "karyotypes.tsv", sep="\t", index=False)
    (d / "lrrs.bed").write_text(f"{CHROM}\t900000\t2400000\tLRR_TEST\n")


def main():
    errors = 0
    OUTDIR.mkdir(parents=True, exist_ok=True)
    build(OUTDIR)

    print("[smoke-trc08] running STEP_TRC_08...")
    r = subprocess.run([
        sys.executable, str(ROOT / "scripts" / "STEP_TRC_08_microhap_hom_het.py"),
        "--microhaps", str(OUTDIR / "microhaps.tsv"),
        "--lrr-bed", str(OUTDIR / "lrrs.bed"),
        "--karyotypes", str(OUTDIR / "karyotypes.tsv"),
        "--outdir", str(OUTDIR),
    ], check=False, capture_output=True, text=True)
    print(r.stdout); print(r.stderr)
    if r.returncode != 0:
        sys.exit("STEP_TRC_08 failed")

    blk = pd.read_csv(OUTDIR / "lrr_microhap_hom_het.tsv", sep="\t").set_index("lrr_id")
    prof = pd.read_csv(OUTDIR / "microhap_hom_het_profile.tsv", sep="\t")
    L = blk.loc["LRR_TEST"]

    def chk(cond, msg):
        nonlocal errors
        print(f"  [{'OK  ' if cond else 'FAIL'}] {msg}")
        errors += not cond

    chk(L["n_hom_samples"] == 16 and L["n_het_samples"] == 10,
        f"group sizes HOM/HET = {L['n_hom_samples']}/{L['n_het_samples']} (want 16/10)")
    chk(L["median_microhapLD_HOM"] > 0.8,
        f"median LD HOM = {L['median_microhapLD_HOM']:.3f} (want > 0.8)")
    chk(L["median_microhapLD_HET"] < 0.3,
        f"median LD HET = {L['median_microhapLD_HET']:.3f} (want < 0.3)")
    chk(L["delta_HOM_minus_HET"] > 0.5,
        f"delta = {L['delta_HOM_minus_HET']:.3f} (want > 0.5)")
    chk(L["interpretation"] == "coherent_HOM_disrupted_HET",
        f"interpretation = {L['interpretation']} (want coherent_HOM_disrupted_HET)")
    chk((prof["mean_LD_to_others_HOM"] > 0.8).all(),
        "every microhap has high HOM mean-LD-to-others (projection works)")

    if errors:
        print(f"\n[smoke-trc08] {errors} failures"); sys.exit(1)
    print("\n[smoke-trc08] ALL PASS")


if __name__ == "__main__":
    main()

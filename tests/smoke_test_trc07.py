#!/usr/bin/env python3
"""
Smoke test for STEP_TRC_07 (microhap block coherence).

Synthesizes 30 samples and two candidate LRRs plus background windows:
  LRR_COHERENT  : every chromosome is purely U-microhaps or purely V-microhaps
                  across all windows -> near-perfect long-distance association
                  -> strong coherence.
  LRR_FRAGMENTED: microhap alleles drawn independently per window per
                  chromosome -> ~no association -> weak coherence.
Background windows are also independent/random -> low background LD.
Asserts coherent >> fragmented and rank labels.
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).parent
ROOT = HERE.parent
OUTDIR = HERE / "tmp_smoke_trc07"
CHROM = "C_gar_LG01"
N = 30


def build(d: Path):
    rng = np.random.default_rng(7)
    rows = []  # sample, chrom, start, end, hap1, hap2

    # arrangement per sample: 0=UU, 1=VV, 2=het
    arr = np.array([0] * 12 + [1] * 12 + [2] * 6)

    # LRR_COHERENT windows 1.0..2.2 Mb, 300 kb spacing
    coh_starts = list(range(1_000_000, 2_300_000, 300_000))
    for st in coh_starts:
        en = st + 5_000
        for i in range(N):
            if arr[i] == 0:
                h1, h2 = "U", "U"
            elif arr[i] == 1:
                h1, h2 = "V", "V"
            else:
                h1, h2 = "U", "V"
            rows.append((f"fish{i:02d}", CHROM, st, en, h1, h2))

    # LRR_FRAGMENTED windows 5.0..6.2 Mb: independent random alleles per chrom
    frag_starts = list(range(5_000_000, 6_300_000, 300_000))
    for st in frag_starts:
        en = st + 5_000
        for i in range(N):
            h1, h2 = rng.choice(["A", "B", "C"], 2)
            rows.append((f"fish{i:02d}", CHROM, st, en, h1, h2))

    # background windows 8.0..12.0 Mb, random
    bg_starts = list(range(8_000_000, 12_100_000, 400_000))
    for st in bg_starts:
        en = st + 5_000
        for i in range(N):
            h1, h2 = rng.choice(["A", "B", "C"], 2)
            rows.append((f"fish{i:02d}", CHROM, st, en, h1, h2))

    pd.DataFrame(rows, columns=["sample", "chrom", "start", "end", "hap1", "hap2"]) \
        .to_csv(d / "microhaps.tsv", sep="\t", index=False)

    (d / "lrrs.bed").write_text(
        f"{CHROM}\t900000\t2400000\tLRR_COHERENT\n"
        f"{CHROM}\t4900000\t6400000\tLRR_FRAGMENTED\n")


def main():
    errors = 0
    OUTDIR.mkdir(parents=True, exist_ok=True)
    build(OUTDIR)

    print("[smoke-trc07] running STEP_TRC_07...")
    r = subprocess.run([
        sys.executable, str(ROOT / "scripts" / "STEP_TRC_07_microhap_ld.py"),
        "--microhaps", str(OUTDIR / "microhaps.tsv"),
        "--lrr-bed", str(OUTDIR / "lrrs.bed"),
        "--min-distance-bp", "100000", "--min-samples", "10",
        "--outdir", str(OUTDIR),
    ], check=False, capture_output=True, text=True)
    print(r.stdout); print(r.stderr)
    if r.returncode != 0:
        sys.exit("STEP_TRC_07 failed")

    df = pd.read_csv(OUTDIR / "lrr_microhap_coherence.tsv", sep="\t").set_index("lrr_id")

    def chk(cond, msg):
        nonlocal errors
        print(f"  [{'OK  ' if cond else 'FAIL'}] {msg}")
        errors += not cond

    coh = df.loc["LRR_COHERENT"]
    frag = df.loc["LRR_FRAGMENTED"]

    chk(coh["median_microhap_LD"] > 0.8,
        f"coherent median LD = {coh['median_microhap_LD']:.3f} (want > 0.8)")
    chk(coh["coherence_rank"] == "strong",
        f"coherent rank = {coh['coherence_rank']} (want strong)")
    chk(frag["median_microhap_LD"] < 0.3,
        f"fragmented median LD = {frag['median_microhap_LD']:.3f} (want < 0.3)")
    chk(frag["coherence_rank"] in ("weak", "low_confidence"),
        f"fragmented rank = {frag['coherence_rank']} (want weak/low_confidence)")
    chk(coh["median_microhap_LD"] > frag["median_microhap_LD"] + 0.5,
        "coherent block ranks well above fragmented region")
    chk(coh["n_longrange_pairs"] >= 3,
        f"coherent has {coh['n_longrange_pairs']} long-range pairs (want >= 3)")

    if errors:
        print(f"\n[smoke-trc07] {errors} failures"); sys.exit(1)
    print("\n[smoke-trc07] ALL PASS")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
STEP_TRC_07_microhap_ld.py — population microhap block-coherence (MicrohapLD)

Submodule role: population-level BLOCK COHERENCE, before/independent of the
pedigree tract-calling layer. Treats each microhap window as a local
multi-allelic haplotype OBJECT and asks, across the cohort:

  Within a candidate LRR, do arrangement-consistent microhaps show stronger
  long-distance association than the genomic background?

For each window pair it computes a multi-allelic LD between the two microhap
loci as bias-corrected Cramer's V over the population's chromosome-resolved
microhap alleles. Per LRR it summarizes long-distance (>= --min-distance-bp)
pairs and compares the block to a matched genome background of pairs lying
outside all LRRs.

PROVENANCE:
  - Upstream (external): MicrohapBuilder (submodule 1) turns phased Clair3
    variants into per-sample per-window microhap objects. NOT in this repo;
    its output is the input contract below.
  - In-repo: this LD / coherence computation.
  - External (popstats / unified-ancestry): GHSL, FIS. The overlay
    (high GHSL + high microhap-LD; +/- FIS) is a downstream join.

PHASE ASSUMPTION (important): hap1/hap2 must be CONSISTENTLY phased across
windows within a region (same phase set), so hap1@A and hap1@B lie on the same
physical chromosome. If phasing is not cross-window consistent this LD is not
interpretable; confirm against the phaser.

Input (--microhaps), MicrohapBuilder output contract:
  sample, chrom, start, end, hap1, hap2   (hap1/hap2 = microhap allele strings)
Input (--lrr-bed): candidate LRRs (chrom,start,end[,name]).
Output (--outdir): lrr_microhap_coherence.tsv (one ranked row per LRR);
  microhap_pair_ld.tsv (per-pair detail) if --emit-pairs.
"""
from __future__ import annotations

import argparse
import sys
import time
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd


def log(msg: str) -> None:
    sys.stderr.write(f"[{time.strftime('%F %T')}] [ngsTracts:trc07] {msg}\n")
    sys.stderr.flush()


def die(msg: str, code: int = 1) -> None:
    log(f"[ERROR] {msg}")
    sys.exit(code)


def cramers_v_corrected(a: np.ndarray, b: np.ndarray) -> tuple[float, int]:
    """Bias-corrected Cramer's V (Bergsma) between two categorical vectors.

    Returns (V, n_obs). a, b are aligned arrays of category labels over the
    population's chromosome-resolved microhap alleles for two windows.
    """
    n = len(a)
    if n < 2:
        return float("nan"), n
    tab = pd.crosstab(pd.Series(a), pd.Series(b)).to_numpy(dtype=float)
    r, c = tab.shape
    if r < 2 or c < 2:
        return 0.0, n  # one window monomorphic -> no association measurable
    row = tab.sum(1, keepdims=True)
    col = tab.sum(0, keepdims=True)
    expected = row @ col / n
    chi2 = float(np.sum((tab - expected) ** 2 / expected))
    phi2 = chi2 / n
    phi2corr = max(0.0, phi2 - (r - 1) * (c - 1) / (n - 1))
    rcorr = r - (r - 1) ** 2 / (n - 1)
    ccorr = c - (c - 1) ** 2 / (n - 1)
    denom = min(rcorr - 1, ccorr - 1)
    if denom <= 0:
        return 0.0, n
    return float(np.sqrt(phi2corr / denom)), n


def pair_ld(win_alleles: dict, wa, wb, min_samples: int):
    """Chromosome-resolved alleles for two windows -> (V, n_obs)."""
    sa = win_alleles[wa]
    sb = win_alleles[wb]
    shared = sa.keys() & sb.keys()
    if len(shared) < min_samples:
        return None
    A, B = [], []
    for s in shared:
        h1a, h2a = sa[s]
        h1b, h2b = sb[s]
        A.append(h1a); B.append(h1b)   # chromosome 1 (hap1 across windows)
        A.append(h2a); B.append(h2b)   # chromosome 2 (hap2 across windows)
    v, n = cramers_v_corrected(np.array(A), np.array(B))
    return v, n


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="ngsTracts STEP_TRC_07 — microhap block coherence")
    p.add_argument("--microhaps", required=True, type=Path)
    p.add_argument("--lrr-bed", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--min-distance-bp", type=int, default=100_000,
                   help="only long-distance window pairs >= this (default 100kb)")
    p.add_argument("--min-samples", type=int, default=10,
                   help="min shared samples to compute a pair LD (default 10)")
    p.add_argument("--high-ld-threshold", type=float, default=0.5,
                   help="V at/above which a pair counts as 'high LD' (default 0.5)")
    p.add_argument("--min-pairs", type=int, default=3,
                   help="min long-range pairs or LRR is low_confidence (default 3)")
    p.add_argument("--max-pairs-per-lrr", type=int, default=5000)
    p.add_argument("--background-max-pairs", type=int, default=2000)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--emit-pairs", action="store_true")
    args = p.parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    mh = pd.read_csv(args.microhaps, sep="\t", dtype=str)
    need = {"sample", "chrom", "start", "end", "hap1", "hap2"}
    if not need.issubset(mh.columns):
        die(f"--microhaps missing columns; need {sorted(need)}")
    mh["start"] = pd.to_numeric(mh["start"]).astype(np.int64)
    mh["end"] = pd.to_numeric(mh["end"]).astype(np.int64)
    mh["mid"] = (mh["start"] + mh["end"]) // 2

    # window -> {sample: (hap1, hap2)}
    win_alleles: dict[tuple, dict] = {}
    win_mid: dict[tuple, int] = {}
    win_chrom: dict[tuple, str] = {}
    for row in mh.itertuples():
        key = (row.chrom, row.start, row.end)
        win_alleles.setdefault(key, {})[row.sample] = (row.hap1, row.hap2)
        win_mid[key] = int(row.mid)
        win_chrom[key] = row.chrom
    log(f"{len(win_alleles):,} microhap windows; {mh['sample'].nunique():,} samples")

    lrr = pd.read_csv(args.lrr_bed, sep="\t", header=None, comment="#", dtype=str)
    lrr = lrr.rename(columns={0: "chrom", 1: "bed_start", 2: "bed_end"})
    lrr["bed_start"] = pd.to_numeric(lrr["bed_start"]).astype(np.int64)
    lrr["bed_end"] = pd.to_numeric(lrr["bed_end"]).astype(np.int64)
    lrr["lrr_id"] = lrr[3].astype(str) if lrr.shape[1] >= 4 else \
        [f"LRR_{i+1:06d}" for i in range(len(lrr))]

    def windows_in(chrom, s, e):
        return [w for w in win_alleles
                if win_chrom[w] == chrom and s <= win_mid[w] < e]

    def lrr_of(w):
        for L in lrr.itertuples():
            if win_chrom[w] == L.chrom and L.bed_start <= win_mid[w] < L.bed_end:
                return L.lrr_id
        return None

    # ---- genome background: long-range pairs with neither window in any LRR ----
    bg_windows = [w for w in win_alleles if lrr_of(w) is None]
    bg_by_chrom: dict[str, list] = {}
    for w in bg_windows:
        bg_by_chrom.setdefault(win_chrom[w], []).append(w)
    bg_pairs = []
    for ch, ws in bg_by_chrom.items():
        for wa, wb in combinations(ws, 2):
            if abs(win_mid[wa] - win_mid[wb]) >= args.min_distance_bp:
                bg_pairs.append((wa, wb))
    if len(bg_pairs) > args.background_max_pairs:
        idx = rng.choice(len(bg_pairs), args.background_max_pairs, replace=False)
        bg_pairs = [bg_pairs[i] for i in idx]
    bg_v = []
    for wa, wb in bg_pairs:
        res = pair_ld(win_alleles, wa, wb, args.min_samples)
        if res is not None and not np.isnan(res[0]):
            bg_v.append(res[0])
    bg_median = float(np.median(bg_v)) if bg_v else float("nan")
    bg_p90 = float(np.percentile(bg_v, 90)) if bg_v else float("nan")
    log(f"background: {len(bg_v):,} long-range pairs, median V={bg_median:.3f}")

    pair_rows = []
    out_rows = []
    for L in lrr.itertuples():
        ws = windows_in(L.chrom, L.bed_start, L.bed_end)
        pairs = [(wa, wb) for wa, wb in combinations(ws, 2)
                 if abs(win_mid[wa] - win_mid[wb]) >= args.min_distance_bp]
        if len(pairs) > args.max_pairs_per_lrr:
            idx = rng.choice(len(pairs), args.max_pairs_per_lrr, replace=False)
            pairs = [pairs[i] for i in idx]
        vs, dists = [], []
        for wa, wb in pairs:
            res = pair_ld(win_alleles, wa, wb, args.min_samples)
            if res is None or np.isnan(res[0]):
                continue
            v, n = res
            d = abs(win_mid[wa] - win_mid[wb])
            vs.append(v); dists.append(d)
            if args.emit_pairs:
                pair_rows.append({"lrr_id": L.lrr_id, "chrom": L.chrom,
                                  "mid_a": win_mid[wa], "mid_b": win_mid[wb],
                                  "distance_bp": d, "n_obs": n, "microhap_LD": v})
        n_pairs = len(vs)
        median_v = float(np.median(vs)) if vs else float("nan")
        p90_v = float(np.percentile(vs, 90)) if vs else float("nan")
        frac_high = float(np.mean(np.array(vs) >= args.high_ld_threshold)) if vs else float("nan")
        slope = float("nan")
        if n_pairs >= 3 and len(set(dists)) >= 2:
            slope = float(np.polyfit(np.array(dists) / 1e6, np.array(vs), 1)[0])  # V per Mb
        delta = (median_v - bg_median) if not np.isnan(bg_median) else float("nan")

        # Elevated = above the background UPPER tail (p90), not its median (which
        # is ~0 for unlinked windows and would make anything positive "moderate").
        if n_pairs < args.min_pairs or len(ws) < 2:
            rank = "low_confidence"
        else:
            bar_bg = bg_p90 if not np.isnan(bg_p90) else args.high_ld_threshold
            elevated = median_v >= bar_bg
            high_abs = median_v >= args.high_ld_threshold
            if elevated and high_abs:
                rank = "strong"
            elif elevated or high_abs:
                rank = "moderate"
            else:
                rank = "weak"

        out_rows.append({
            "lrr_id": L.lrr_id, "chrom": L.chrom, "start": L.bed_start, "end": L.bed_end,
            "n_windows": len(ws), "n_longrange_pairs": n_pairs,
            "median_microhap_LD": median_v, "p90_microhap_LD": p90_v,
            "frac_pairs_high_LD": frac_high, "LD_decay_slope_per_Mb": slope,
            "background_median_LD": bg_median, "background_p90_LD": bg_p90,
            "delta_vs_background": delta, "coherence_rank": rank,
        })

    out = pd.DataFrame(out_rows).sort_values(
        ["median_microhap_LD", "n_longrange_pairs"], ascending=[False, False],
        na_position="last").reset_index(drop=True)
    out_path = args.outdir / "lrr_microhap_coherence.tsv"
    out.to_csv(out_path, sep="\t", index=False)
    log(f"wrote {out_path}  ({len(out)} LRRs)")
    if args.emit_pairs:
        pd.DataFrame(pair_rows).to_csv(args.outdir / "microhap_pair_ld.tsv", sep="\t", index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())

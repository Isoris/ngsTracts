#!/usr/bin/env python3
"""
STEP_TRC_08_microhap_hom_het.py — MicrohapLD_HOMHET

Within each candidate LRR, split samples into arrangement-HOMOZYGOUS and
arrangement-HETEROZYGOUS (karyotype from upstream), then compute population
microhap-microhap LD SEPARATELY in each group and compare. Each microhap is
projected against all other microhaps in the same LRR; the per-microhap mean
LD and the block-level median LD are reported for HOM and HET.

Interpretation (block level):
  HOM high & HET high          -> coherent in both genotype classes
  HOM high & HET reduced       -> homozygous arrangements coherent, heterozygotes
                                  partially mixed
  HOM high & HET very low      -> heterozygotes disrupted (recombination /
                                  NCO-like local disruption / bad boundary / noise)

PROVENANCE:
  - Upstream (external): MicrohapBuilder (phased Clair3 -> microhap objects) and
    the per-sample x LRR karyotype call (HOM_REF/HET/HOM_INV). NOT in this repo.
  - In-repo: this HOM/HET LD comparison.
Phase assumption is the same as STEP_TRC_07: hap1/hap2 consistently phased
across windows.

Inputs:
  --microhaps    sample, chrom, start, end, hap1, hap2 [, microhap_id]
  --lrr-bed      candidate LRRs (chrom,start,end[,name])
  --karyotypes   sample, lrr_id, karyotype   (0/1/2 or HOM_REF/HET/HOM_INV)
Output (--outdir):
  lrr_microhap_hom_het.tsv         block-level HOM vs HET
  microhap_hom_het_profile.tsv     per-microhap mean LD to others, HOM vs HET
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
    sys.stderr.write(f"[{time.strftime('%F %T')}] [ngsTracts:trc08] {msg}\n")
    sys.stderr.flush()


def die(msg: str, code: int = 1) -> None:
    log(f"[ERROR] {msg}")
    sys.exit(code)


def cramers_v_corrected(a: np.ndarray, b: np.ndarray) -> float:
    n = len(a)
    if n < 2:
        return float("nan")
    tab = pd.crosstab(pd.Series(a), pd.Series(b)).to_numpy(dtype=float)
    r, c = tab.shape
    if r < 2 or c < 2:
        return 0.0
    expected = tab.sum(1, keepdims=True) @ tab.sum(0, keepdims=True) / n
    chi2 = float(np.sum((tab - expected) ** 2 / expected))
    phi2 = chi2 / n
    phi2corr = max(0.0, phi2 - (r - 1) * (c - 1) / (n - 1))
    rcorr = r - (r - 1) ** 2 / (n - 1)
    ccorr = c - (c - 1) ** 2 / (n - 1)
    denom = min(rcorr - 1, ccorr - 1)
    return 0.0 if denom <= 0 else float(np.sqrt(phi2corr / denom))


def pair_ld_subset(win_alleles, wa, wb, samples, min_samples):
    """Chromosome-resolved Cramer's V for two windows over a sample subset."""
    sa, sb = win_alleles[wa], win_alleles[wb]
    shared = (sa.keys() & sb.keys()) & samples
    if len(shared) < min_samples:
        return float("nan")
    A, B = [], []
    for s in shared:
        h1a, h2a = sa[s]; h1b, h2b = sb[s]
        A.append(h1a); B.append(h1b)
        A.append(h2a); B.append(h2b)
    return cramers_v_corrected(np.array(A), np.array(B))


def load_karyotypes(path: Path) -> dict[tuple, str]:
    """(sample, lrr_id) -> 'HOM' | 'HET'."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    need = {"sample", "lrr_id", "karyotype"}
    if not need.issubset(df.columns):
        die(f"--karyotypes needs columns {sorted(need)}")
    grp = {}
    for r in df.itertuples():
        k = str(r.karyotype).strip()
        if k in {"1", "HET"}:
            g = "HET"
        elif k in {"0", "2", "HOM_REF", "HOM_INV", "HOM"}:
            g = "HOM"
        else:
            die(f"unrecognized karyotype value: {k!r} (use 0/1/2 or HOM_*/HET)")
        grp[(r.sample, r.lrr_id)] = g
    return grp


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="ngsTracts STEP_TRC_08 — microhap HOM/HET LD")
    p.add_argument("--microhaps", required=True, type=Path)
    p.add_argument("--lrr-bed", required=True, type=Path)
    p.add_argument("--karyotypes", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--min-distance-bp", type=int, default=0,
                   help="minimum window-pair distance (default 0 = all within-LRR pairs)")
    p.add_argument("--min-samples", type=int, default=6,
                   help="min group samples sharing a pair to compute its LD (default 6)")
    p.add_argument("--min-group-samples", type=int, default=4,
                   help="min HOM/HET samples or that group is NA (default 4)")
    p.add_argument("--high-ld-threshold", type=float, default=0.5)
    p.add_argument("--low-ld-threshold", type=float, default=0.2)
    args = p.parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)

    mh = pd.read_csv(args.microhaps, sep="\t", dtype=str)
    need = {"sample", "chrom", "start", "end", "hap1", "hap2"}
    if not need.issubset(mh.columns):
        die(f"--microhaps missing columns; need {sorted(need)}")
    mh["start"] = pd.to_numeric(mh["start"]).astype(np.int64)
    mh["end"] = pd.to_numeric(mh["end"]).astype(np.int64)
    mh["mid"] = (mh["start"] + mh["end"]) // 2
    has_id = "microhap_id" in mh.columns

    win_alleles, win_mid, win_chrom, win_id = {}, {}, {}, {}
    for row in mh.itertuples():
        key = (row.chrom, row.start, row.end)
        win_alleles.setdefault(key, {})[row.sample] = (row.hap1, row.hap2)
        win_mid[key] = int(row.mid)
        win_chrom[key] = row.chrom
        win_id[key] = getattr(row, "microhap_id") if has_id else f"{row.chrom}:{row.start}-{row.end}"

    lrr = pd.read_csv(args.lrr_bed, sep="\t", header=None, comment="#", dtype=str)
    lrr = lrr.rename(columns={0: "chrom", 1: "bed_start", 2: "bed_end"})
    lrr["bed_start"] = pd.to_numeric(lrr["bed_start"]).astype(np.int64)
    lrr["bed_end"] = pd.to_numeric(lrr["bed_end"]).astype(np.int64)
    lrr["lrr_id"] = lrr[3].astype(str) if lrr.shape[1] >= 4 else \
        [f"LRR_{i+1:06d}" for i in range(len(lrr))]

    karyo = load_karyotypes(args.karyotypes)

    block_rows, profile_rows = [], []
    for L in lrr.itertuples():
        ws = [w for w in win_alleles
              if win_chrom[w] == L.chrom and L.bed_start <= win_mid[w] < L.bed_end]
        hom = {s for (s, lid), g in karyo.items() if lid == L.lrr_id and g == "HOM"}
        het = {s for (s, lid), g in karyo.items() if lid == L.lrr_id and g == "HET"}

        pairs = [(wa, wb) for wa, wb in combinations(ws, 2)
                 if abs(win_mid[wa] - win_mid[wb]) >= args.min_distance_bp]

        v_hom, v_het = {}, {}   # window-pair -> V, per group
        for wa, wb in pairs:
            if len(hom) >= args.min_group_samples:
                v_hom[(wa, wb)] = pair_ld_subset(win_alleles, wa, wb, hom, args.min_samples)
            if len(het) >= args.min_group_samples:
                v_het[(wa, wb)] = pair_ld_subset(win_alleles, wa, wb, het, args.min_samples)

        def med(d):
            vals = [v for v in d.values() if not np.isnan(v)]
            return float(np.median(vals)) if vals else float("nan")

        m_hom, m_het = med(v_hom), med(v_het)
        delta = (m_hom - m_het) if not (np.isnan(m_hom) or np.isnan(m_het)) else float("nan")

        # per-microhap mean LD to others (projection)
        for w in ws:
            def mean_to_others(d):
                vals = [v for (a, b), v in d.items()
                        if (a == w or b == w) and not np.isnan(v)]
                return float(np.mean(vals)) if vals else float("nan")
            profile_rows.append({
                "lrr_id": L.lrr_id, "microhap_id": win_id[w], "chrom": L.chrom,
                "mid": win_mid[w],
                "mean_LD_to_others_HOM": mean_to_others(v_hom),
                "mean_LD_to_others_HET": mean_to_others(v_het),
            })

        n_pairs_hom = sum(1 for v in v_hom.values() if not np.isnan(v))
        n_pairs_het = sum(1 for v in v_het.values() if not np.isnan(v))
        thr, low = args.high_ld_threshold, args.low_ld_threshold
        if len(hom) < args.min_group_samples or n_pairs_hom == 0:
            interp = "low_confidence_HOM"
        else:
            hom_high = m_hom >= thr
            het_ok = not np.isnan(m_het) and n_pairs_het > 0
            if not het_ok:
                interp = "coherent_HOM_no_HET_data" if hom_high else "weak_HOM_no_HET_data"
            elif hom_high and m_het >= thr:
                interp = "coherent_in_both"
            elif hom_high and m_het >= low:
                interp = "coherent_HOM_reduced_HET"
            elif hom_high and m_het < low:
                interp = "coherent_HOM_disrupted_HET"
            elif (not hom_high) and m_het >= thr:
                interp = "coherent_HET_only"
            else:
                interp = "weak_both"

        block_rows.append({
            "lrr_id": L.lrr_id, "chrom": L.chrom, "start": L.bed_start, "end": L.bed_end,
            "n_windows": len(ws), "n_hom_samples": len(hom), "n_het_samples": len(het),
            "n_pairs_hom": n_pairs_hom, "n_pairs_het": n_pairs_het,
            "median_microhapLD_HOM": m_hom, "median_microhapLD_HET": m_het,
            "delta_HOM_minus_HET": delta, "interpretation": interp,
        })

    pd.DataFrame(block_rows).sort_values("delta_HOM_minus_HET", ascending=False,
                                         na_position="last").to_csv(
        args.outdir / "lrr_microhap_hom_het.tsv", sep="\t", index=False)
    pd.DataFrame(profile_rows).to_csv(
        args.outdir / "microhap_hom_het_profile.tsv", sep="\t", index=False)
    log(f"wrote lrr_microhap_hom_het.tsv ({len(block_rows)} LRRs) "
        f"and microhap_hom_het_profile.tsv ({len(profile_rows)} microhaps)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

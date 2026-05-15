#!/usr/bin/env python3
"""
STEP_TRC_02_traversal_scan.py — per-site sliding-window CO breakpoint refinement

Optional second pass on top of STEP_TRC_01. For each CO call from
tract_classifications.tsv, slides a window over per_site_haplotype_calls.tsv
and locates the maximum concordance gradient. Emits refined breakpoints
with confidence intervals.

Only runs on intervals classified as CO. Skips NCO/DCO/etc — gene
conversion and double crossover events have no breakpoint to refine.

Reads per_site_haplotype_calls.tsv from Stage 3 (must have been emitted
with --emit-per-site).
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd


def log(msg: str) -> None:
    sys.stderr.write(f"[{time.strftime('%F %T')}] [ngsTracts:trc02] {msg}\n")
    sys.stderr.flush()


def die(msg: str, code: int = 1) -> None:
    log(f"[ERROR] {msg}")
    sys.exit(code)


def find_max_gradient_breakpoint(positions: np.ndarray, states: np.ndarray,
                                 window: int = 21):
    """
    Slide a centered window across the per-site inferred_state sequence.
    Compute fraction of sites in window that are 'hapA' (states encoded
    as 0=hapA, 1=hapB). Breakpoint = position of maximum |gradient|.

    Returns (breakpoint_pos, ci_left, ci_right, gradient_magnitude) or
    None if too few sites.
    """
    n = len(positions)
    if n < window + 1:
        return None
    half = window // 2

    # Fraction-hapB in centered window, vectorized via cumulative sum
    csum = np.concatenate([[0], np.cumsum(states.astype(np.int64))])
    win_starts = np.arange(half, n - half)
    win_sums = csum[win_starts + half + 1] - csum[win_starts - half]
    frac_B = win_sums / window

    if len(frac_B) < 2:
        return None

    # Gradient between consecutive windows (centered on the right index)
    grad = np.diff(frac_B)
    if len(grad) == 0 or np.all(grad == 0):
        return None
    idx = int(np.argmax(np.abs(grad)))
    grad_mag = float(np.abs(grad[idx]))
    # idx is into grad which lives between windows; map back to a position
    # The breakpoint is between win_starts[idx] and win_starts[idx + 1]
    left_window_center = win_starts[idx]
    right_window_center = win_starts[idx + 1]
    bp_pos = int((positions[left_window_center] + positions[right_window_center]) / 2)

    # CI: positions where |gradient| >= grad_mag - 0.1
    thresh = max(grad_mag - 0.1, 0.0)
    big = np.where(np.abs(grad) >= thresh)[0]
    if len(big) == 0:
        ci_left = bp_pos
        ci_right = bp_pos
    else:
        ci_left = int(positions[win_starts[big[0]]])
        ci_right = int(positions[win_starts[big[-1]] + 1])

    return bp_pos, ci_left, ci_right, grad_mag


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="ngsTracts STEP_TRC_02 — refine CO breakpoints")
    p.add_argument("--tract-classifications", required=True, type=Path,
                   help="Output of STEP_TRC_01 (tract_classifications.tsv)")
    p.add_argument("--per-site-file", required=True, type=Path,
                   help="ngsPedigree Stage 3 per_site_haplotype_calls.tsv")
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--window-size", type=int, default=21,
                   help="sliding window size in informative sites (default 21)")
    args = p.parse_args(argv)

    args.outdir.mkdir(parents=True, exist_ok=True)

    log(f"loading tract classifications from {args.tract_classifications}")
    tc = pd.read_csv(args.tract_classifications, sep="\t", dtype=str)
    if "class" not in tc.columns:
        die("input tract_classifications.tsv missing 'class' column")
    co_calls = tc[tc["class"] == "CO"].copy()
    if len(co_calls) == 0:
        log("no CO calls to refine; writing empty output")
        out = args.outdir / "traversal_breakpoints.tsv"
        pd.DataFrame(columns=[
            "interval_id", "parent_id", "offspring_id", "chrom",
            "refined_breakpoint_bp", "ci_left_bp", "ci_right_bp",
            "gradient_magnitude", "n_window_sites",
        ]).to_csv(out, sep="\t", index=False)
        return 0
    log(f"{len(co_calls):,} CO calls to refine")

    # Load per-site once; index for fast filtering
    log(f"loading per-site calls from {args.per_site_file}")
    ps = pd.read_csv(args.per_site_file, sep="\t", dtype={
        "parent_id": str, "offspring_id": str, "chrom": str,
        "pos": np.int64, "inferred_state": str, "departure": str,
    }, usecols=["parent_id", "offspring_id", "chrom", "pos", "inferred_state"])
    log(f"loaded {len(ps):,} per-site rows")

    # Encode states: 0=hapA, 1=hapB. Drop unphased.
    ps = ps[ps["inferred_state"].isin(["hapA", "hapB"])].copy()
    ps["state_int"] = (ps["inferred_state"] == "hapB").astype(np.int8)
    ps = ps.sort_values(["parent_id", "offspring_id", "chrom", "pos"]).reset_index(drop=True)

    # Index for grouped lookups
    indexed = ps.set_index(["parent_id", "offspring_id", "chrom"], drop=False)

    co_calls["start"] = pd.to_numeric(co_calls["start"]).astype(np.int64)
    co_calls["end"]   = pd.to_numeric(co_calls["end"]).astype(np.int64)

    refined = []
    n_ok = 0
    n_skip = 0
    for r in co_calls.itertuples():
        key = (r.parent_id, r.offspring_id, r.chrom)
        try:
            sub = indexed.loc[key]
        except KeyError:
            n_skip += 1
            continue
        if isinstance(sub, pd.Series):  # only one row matched
            n_skip += 1
            continue
        # Pull sites in a generous window around the interval (interval +/- 100kb)
        pad = 100_000
        lo = max(1, int(r.start) - pad)
        hi = int(r.end) + pad
        local = sub[(sub["pos"] >= lo) & (sub["pos"] <= hi)]
        if len(local) < args.window_size + 1:
            n_skip += 1
            continue

        result = find_max_gradient_breakpoint(
            local["pos"].to_numpy(),
            local["state_int"].to_numpy(),
            window=args.window_size,
        )
        if result is None:
            n_skip += 1
            continue
        bp, ciL, ciR, gm = result
        refined.append({
            "interval_id": r.interval_id,
            "parent_id":   r.parent_id,
            "offspring_id": r.offspring_id,
            "chrom":       r.chrom,
            "refined_breakpoint_bp": bp,
            "ci_left_bp":  ciL,
            "ci_right_bp": ciR,
            "gradient_magnitude": gm,
            "n_window_sites": int(len(local)),
        })
        n_ok += 1

    out_df = pd.DataFrame(refined)
    out_path = args.outdir / "traversal_breakpoints.tsv"
    out_df.to_csv(out_path, sep="\t", index=False)
    log(f"wrote {out_path}: refined={n_ok}, skipped={n_skip}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

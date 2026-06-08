#!/usr/bin/env python3
"""
STEP_TRC_06_lrr_suppression.py — rank candidate low-recombination regions (LRRs)
by crossover suppression, with a POPULATION-level support score.

For each candidate LRR, split into edge zones vs interior, count CO-like
switches inside the interior, and ask across the cohort: of the dyads that
actually had the OPPORTUNITY to show a CO in this interior (enough informative
sites spanning it, or an observed event there), how many show no interior CO?

Family-level "no CO" is weak on its own (a region with no informative
transmissions trivially has zero COs). So the headline statistic is a
conservative population support score:

    population_support_score = Wilson lower 95% bound of
        (n opportunity-dyads with NO interior CO) / (n opportunity-dyads)

One consistent family scores low (wide CI); many consistent families score
high. Regions are ranked by this score.

This does NOT compute GHSL / FIS / divergence — those are external
(popstats / unified-ancestry). STEP_TRC_06 supplies the CO-suppression ranking
to overlay against them.

Inputs:
  --lrr-bed                 candidate LRRs, BED (chrom,start,end[,name]); 0-based
  --tract-classifications   STEP_TRC_01 output (CO/NCO/DCO calls)
Optional:
  --per-site-file           Stage 3 per_site_haplotype_calls.tsv -> real
                            informative-opportunity denominator (recommended)

Output (to --outdir): lrr_co_suppression.tsv (one row per LRR, ranked).
"""
from __future__ import annotations

import argparse
import sys
import time
from math import sqrt
from pathlib import Path

import numpy as np
import pandas as pd


def log(msg: str) -> None:
    sys.stderr.write(f"[{time.strftime('%F %T')}] [ngsTracts:trc06] {msg}\n")
    sys.stderr.flush()


def die(msg: str, code: int = 1) -> None:
    log(f"[ERROR] {msg}")
    sys.exit(code)


def wilson_lower(k: int, n: int, z: float = 1.96) -> float:
    """Wilson score lower confidence bound for a binomial proportion k/n."""
    if n == 0:
        return float("nan")
    phat = k / n
    denom = 1.0 + z * z / n
    center = phat + z * z / (2 * n)
    margin = z * sqrt(phat * (1 - phat) / n + z * z / (4 * n * n))
    return max(0.0, (center - margin) / denom)


def load_lrr(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, comment="#", dtype=str)
    if df.shape[1] < 3:
        die("--lrr-bed needs at least chrom,start,end")
    df = df.rename(columns={0: "chrom", 1: "bed_start", 2: "bed_end"})
    df["bed_start"] = pd.to_numeric(df["bed_start"]).astype(np.int64)
    df["bed_end"] = pd.to_numeric(df["bed_end"]).astype(np.int64)
    if df.shape[1] >= 4:
        df["lrr_id"] = df[3].astype(str)
    else:
        df["lrr_id"] = [f"LRR_{i+1:06d}" for i in range(len(df))]
    return df[["lrr_id", "chrom", "bed_start", "bed_end"]]


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="ngsTracts STEP_TRC_06 — LRR CO-suppression ranking")
    p.add_argument("--lrr-bed", required=True, type=Path)
    p.add_argument("--tract-classifications", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--per-site-file", type=Path, default=None,
                   help="Stage 3 per_site_haplotype_calls.tsv for opportunity denominator")
    p.add_argument("--edge-frac", type=float, default=0.10,
                   help="fraction of LRR length at each end treated as edge (default 0.10)")
    p.add_argument("--edge-bp", type=int, default=None,
                   help="fixed edge width in bp (overrides --edge-frac if larger)")
    p.add_argument("--min-opportunity-sites", type=int, default=10,
                   help="informative interior sites for a dyad to 'have a chance' (default 10)")
    p.add_argument("--min-opportunity-dyads", type=int, default=2,
                   help="min opportunity dyads or the LRR is low_confidence (default 2)")
    p.add_argument("--many-co-threshold", type=int, default=2,
                   help="interior CO count at/above which an LRR is 'many_CO_inside' (default 2)")
    p.add_argument("--z", type=float, default=1.96, help="Wilson z (default 1.96 = 95%)")
    args = p.parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)

    lrr = load_lrr(args.lrr_bed)

    tc = pd.read_csv(args.tract_classifications, sep="\t", dtype=str)
    for c in ("parent_id", "offspring_id", "chrom", "start", "end", "class"):
        if c not in tc.columns:
            die(f"tract_classifications.tsv missing column: {c}")
    tc["start"] = pd.to_numeric(tc["start"]).astype(np.int64)
    tc["end"] = pd.to_numeric(tc["end"]).astype(np.int64)
    tc["mid"] = (tc["start"] + tc["end"]) / 2.0
    tc["dyad"] = tc["parent_id"] + "||" + tc["offspring_id"]

    # Opportunity from per-site informative coverage (recommended).
    persite_idx: dict[tuple, np.ndarray] = {}
    opportunity_source = "proxy(any-event-on-chrom)"
    if args.per_site_file is not None and args.per_site_file.exists():
        ps = pd.read_csv(args.per_site_file, sep="\t",
                         usecols=["parent_id", "offspring_id", "chrom", "pos", "inferred_state"],
                         dtype={"parent_id": str, "offspring_id": str, "chrom": str,
                                "pos": np.int64, "inferred_state": str})
        ps = ps[ps["inferred_state"].isin(["hapA", "hapB"])]
        for (pa, of, ch), g in ps.groupby(["parent_id", "offspring_id", "chrom"]):
            persite_idx[(pa, of, ch)] = np.sort(g["pos"].to_numpy())
        opportunity_source = "per-site informative coverage"
        log(f"loaded per-site opportunity for {len(persite_idx):,} dyad-chroms")
    else:
        log("no --per-site-file; using weak proxy opportunity (any event on chrom)")

    rows = []
    for L in lrr.itertuples():
        s1, e1 = L.bed_start + 1, L.bed_end          # 1-based inclusive
        length = e1 - s1 + 1
        edge = int(round(args.edge_frac * length))
        if args.edge_bp is not None:
            edge = max(edge, args.edge_bp)
        ilo, ihi = s1 + edge, e1 - edge
        interior_ok = ilo <= ihi

        ev = tc[tc["chrom"] == L.chrom]
        co = ev[ev["class"] == "CO"]
        in_int = lambda d: d[(d["mid"] >= ilo) & (d["mid"] <= ihi)] if interior_ok else d.iloc[0:0]
        in_edge = lambda d: d[((d["mid"] >= s1) & (d["mid"] < ilo)) |
                              ((d["mid"] > ihi) & (d["mid"] <= e1))]
        co_int = in_int(co)
        n_interior_CO = int(len(co_int))
        n_edge_CO = int(len(in_edge(co)))
        n_interior_NCO = int(len(in_int(ev[ev["class"] == "NCO"])))
        n_interior_DCO = int(len(in_int(ev[ev["class"] == "DCO"])))
        dyads_with_interior_CO = set(co_int["dyad"].unique())

        # Opportunity dyads = enough interior informative sites, OR an observed
        # interior event (a dyad that produced any interior event clearly had a chance).
        opp_dyads = set(dyads_with_interior_CO)
        interior_site_counts = []
        if persite_idx and interior_ok:
            seen = set(k for k in persite_idx if k[2] == L.chrom)
            for (pa, of, ch) in seen:
                pos = persite_idx[(pa, of, ch)]
                cnt = int(np.searchsorted(pos, ihi, "right") - np.searchsorted(pos, ilo, "left"))
                if cnt >= args.min_opportunity_sites:
                    opp_dyads.add(f"{pa}||{of}")
                    interior_site_counts.append(cnt)
        else:
            # proxy: any dyad with any event on this chrom
            opp_dyads |= set(ev["dyad"].unique())
        # also fold in interior events of any class as opportunity evidence
        opp_dyads |= set(in_int(ev)["dyad"].unique())

        n_opp = len(opp_dyads)
        n_no_interior_CO = len(opp_dyads - dyads_with_interior_CO)
        consistency = (n_no_interior_CO / n_opp) if n_opp else float("nan")
        score = wilson_lower(n_no_interior_CO, n_opp, args.z) if n_opp else float("nan")

        # pattern label
        if not interior_ok or n_opp < args.min_opportunity_dyads:
            pattern = "low_confidence"
        elif n_interior_CO == 0 and n_edge_CO > 0:
            pattern = "CO_at_boundary_only"
        elif n_interior_CO == 0:
            pattern = "no_CO_inside"
        elif n_interior_CO >= args.many_co_threshold:
            pattern = "many_CO_inside"
        else:
            pattern = "few_CO_inside"

        rows.append({
            "lrr_id": L.lrr_id, "chrom": L.chrom, "start": s1, "end": e1,
            "length_bp": length, "edge_bp": edge,
            "pattern": pattern,
            "population_support_score": score,
            "suppression_consistency": consistency,
            "n_dyads_with_opportunity": n_opp,
            "n_dyads_no_interior_CO": n_no_interior_CO,
            "opportunity_source": opportunity_source,
            "mean_interior_sites_per_dyad": (float(np.mean(interior_site_counts))
                                             if interior_site_counts else float("nan")),
            "n_interior_CO": n_interior_CO,
            "n_interior_NCO": n_interior_NCO,
            "n_interior_DCO": n_interior_DCO,
            "n_edge_CO": n_edge_CO,
            "flag_low_opportunity": int(n_opp < args.min_opportunity_dyads),
        })

    out = pd.DataFrame(rows).sort_values(
        ["population_support_score", "n_dyads_with_opportunity"],
        ascending=[False, False], na_position="last").reset_index(drop=True)
    out_path = args.outdir / "lrr_co_suppression.tsv"
    out.to_csv(out_path, sep="\t", index=False)
    log(f"wrote {out_path}  ({len(out)} LRRs; opportunity from {opportunity_source})")
    return 0


if __name__ == "__main__":
    sys.exit(main())

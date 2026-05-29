#!/usr/bin/env python3
"""
STEP_TRC_03_exchange_mask.py — build an exchange/gene-flux segment mask

Downstream helper on top of STEP_TRC_01. Reads tract_classifications.tsv and
emits a BED of putative meiotic-exchange / gene-flux segments (the NCO / DCO /
MOSAIC_* classes), so they can be masked or separately reported when computing
arrangement-specific diversity (pi, DAF, dXY) on inversion candidates.

IMPORTANT scope:
  - These segments are per-(parent, offspring) OBSERVED meiotic events from the
    pedigree, not population-level gene-flux inferred from diversity. Aggregate
    across dyads (support count) to find exchange-PRONE regions; treat as a
    prior / cross-check for a diversity-derived mask, not a 1:1 truth set.
  - This step does NOT compute pi/DAF/dXY. It only produces the segment list.

Outputs (to --outdir):
  exchange_segments.bed         per-tract BED6 (one row per dyad-tract)
  exchange_segments.merged.bed  merged across dyads; name = support count
  exchange_mask_summary.tsv     per-class + (if --inversion-atlas) per-inversion
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd


# Classes that represent a sequence exchange / gene-flux footprint.
# CO is a clean crossover (handled by STEP_TRC_02 breakpoint refinement, not a
# tract to mask); AMBIG / LOW_CONFIDENCE are not interpretable exchanges.
EXCHANGE_CLASSES_DEFAULT = "NCO,DCO,MOSAIC_SHORT,MOSAIC_LONG"

CONF_SCORE = {"high": 1000, "medium": 500, "low": 100}


def log(msg: str) -> None:
    sys.stderr.write(f"[{time.strftime('%F %T')}] [ngsTracts:trc03] {msg}\n")
    sys.stderr.flush()


def die(msg: str, code: int = 1) -> None:
    log(f"[ERROR] {msg}")
    sys.exit(code)


def load_atlas(path: Path) -> pd.DataFrame:
    """Read inversion atlas (needs inversion_id, chrom, start, end)."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    need = {"inversion_id", "chrom", "start", "end"}
    if not need.issubset(df.columns):
        die(f"inversion atlas missing columns; need {sorted(need)}")
    df["start"] = pd.to_numeric(df["start"], errors="raise").astype(np.int64)
    df["end"] = pd.to_numeric(df["end"], errors="raise").astype(np.int64)
    return df[["inversion_id", "chrom", "start", "end"]]


def merge_intervals(seg: pd.DataFrame) -> pd.DataFrame:
    """
    Merge overlapping/abutting BED intervals per chrom (0-based, half-open).
    Returns merged intervals with a support count = number of input tracts
    that contributed to each merged block.
    """
    out_rows = []
    for chrom, g in seg.sort_values(["chrom", "bed_start", "bed_end"]).groupby("chrom"):
        cur_s = cur_e = None
        support = 0
        for s, e in zip(g["bed_start"].to_numpy(), g["bed_end"].to_numpy()):
            if cur_s is None:
                cur_s, cur_e, support = s, e, 1
            elif s <= cur_e:                      # overlap or abut
                cur_e = max(cur_e, e)
                support += 1
            else:
                out_rows.append((chrom, cur_s, cur_e, support))
                cur_s, cur_e, support = s, e, 1
        if cur_s is not None:
            out_rows.append((chrom, cur_s, cur_e, support))
    return pd.DataFrame(out_rows, columns=["chrom", "bed_start", "bed_end", "support"])


def assign_inversions(seg: pd.DataFrame, atlas: pd.DataFrame) -> pd.DataFrame:
    """For each tract, list overlapping inversion_id(s) (1-based-inclusive atlas)."""
    inv_of = []
    for r in seg.itertuples():
        hits = atlas[
            (atlas["chrom"] == r.chrom)
            & (atlas["start"] <= r.end_1based)      # 1-based inclusive overlap
            & (atlas["end"] >= r.start_1based)
        ]["inversion_id"].tolist()
        inv_of.append(",".join(hits) if hits else "-")
    seg = seg.copy()
    seg["inversion_id"] = inv_of
    return seg


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description="ngsTracts STEP_TRC_03 — build exchange/gene-flux segment mask")
    p.add_argument("--tract-classifications", required=True, type=Path,
                   help="Output of STEP_TRC_01 (tract_classifications.tsv)")
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--classes", default=EXCHANGE_CLASSES_DEFAULT,
                   help=f"comma list of classes to treat as exchange tracts "
                        f"(default {EXCHANGE_CLASSES_DEFAULT})")
    p.add_argument("--inside-only", action="store_true",
                   help="keep only tracts with inside_inversion in {yes, partial}")
    p.add_argument("--pad", type=int, default=0,
                   help="expand each tract by this many bp on both sides "
                        "before merging (breakpoint buffer; default 0)")
    p.add_argument("--inversion-atlas", type=Path, default=None,
                   help="optional atlas (inversion_id,chrom,start,end) for "
                        "per-inversion aggregation")
    p.add_argument("--per-sample", action="store_true",
                   help="also emit one BED per carrier sample (for GL-based "
                        "per-sample masking, where a global BED would drop a "
                        "site for all samples)")
    p.add_argument("--sample-key", choices=["offspring", "parent", "both"],
                   default="offspring",
                   help="which sample carries the tract for --per-sample. The "
                        "offspring inherits the exchanged sequence, so it is "
                        "the carrier (default offspring)")
    args = p.parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)

    keep = {c.strip() for c in args.classes.split(",") if c.strip()}
    log(f"exchange classes: {sorted(keep)}  inside_only={args.inside_only}  pad={args.pad}")

    tc = pd.read_csv(args.tract_classifications, sep="\t", dtype=str)
    for col in ("class", "chrom", "start", "end", "interval_id"):
        if col not in tc.columns:
            die(f"tract_classifications.tsv missing required column: {col}")
    tc["start"] = pd.to_numeric(tc["start"], errors="raise").astype(np.int64)
    tc["end"] = pd.to_numeric(tc["end"], errors="raise").astype(np.int64)

    seg = tc[tc["class"].isin(keep)].copy()
    if args.inside_only and "inside_inversion" in seg.columns:
        seg = seg[seg["inside_inversion"].isin(["yes", "partial"])].copy()
    log(f"selected {len(seg):,} exchange tracts (of {len(tc):,} classified intervals)")

    # 1-based-inclusive -> 0-based half-open BED, with optional pad
    seg["start_1based"] = seg["start"]
    seg["end_1based"] = seg["end"]
    seg["bed_start"] = np.maximum(seg["start"] - 1 - args.pad, 0)
    seg["bed_end"] = seg["end"] + args.pad

    if args.inversion_atlas is not None:
        seg = assign_inversions(seg, load_atlas(args.inversion_atlas))

    # Reusable BED6 name/score columns.
    if len(seg):
        seg["_score"] = (seg["confidence"].map(CONF_SCORE).fillna(0).astype(int)
                         if "confidence" in seg.columns else 0)
        seg["_name"] = (seg["class"] + "|" + seg["interval_id"].astype(str)
                        + "|" + seg.get("parent_id", "?").astype(str)
                        + ">" + seg.get("offspring_id", "?").astype(str))

    # ---- per-tract BED6 (global; pooled over all samples) ----
    bed_path = args.outdir / "exchange_segments.bed"
    if len(seg):
        pd.DataFrame({
            "chrom": seg["chrom"], "start": seg["bed_start"], "end": seg["bed_end"],
            "name": seg["_name"], "score": seg["_score"], "strand": ".",
        }).to_csv(bed_path, sep="\t", index=False, header=False)
    else:
        bed_path.write_text("")
    log(f"wrote {bed_path}")

    # ---- per-sample BEDs (for GL-based per-sample masking) ----
    # A global BED would drop a site for ALL samples; with genotype likelihoods
    # we want to exclude only the carrier sample's contribution at these sites.
    # The offspring inherits the exchanged sequence and is the carrier.
    if args.per_sample:
        roles = {"offspring": ["offspring_id"], "parent": ["parent_id"],
                 "both": ["offspring_id", "parent_id"]}[args.sample_key]
        for r in roles:
            if r not in seg.columns:
                die(f"--per-sample needs column {r} in tract_classifications.tsv")
        by_dir = args.outdir / "by_sample"
        by_dir.mkdir(parents=True, exist_ok=True)
        long_cols = ["sample_id", "role", "chrom", "bed_start", "bed_end",
                     "class", "interval_id", "parent_id", "offspring_id"]
        parts = []
        if len(seg):
            for r in roles:
                role_name = "offspring" if r == "offspring_id" else "parent"
                part = seg.copy()
                part["sample_id"] = part[r]
                part["role"] = role_name
                parts.append(part)
        long = pd.concat(parts, ignore_index=True) if parts else \
            pd.DataFrame(columns=long_cols)

        for sid, g in long.groupby("sample_id"):
            pd.DataFrame({
                "chrom": g["chrom"], "start": g["bed_start"], "end": g["bed_end"],
                "name": g["_name"], "score": g["_score"], "strand": ".",
            }).sort_values(["chrom", "start"]).to_csv(
                by_dir / f"{sid}.exchange.bed", sep="\t", index=False, header=False)

        long_path = args.outdir / "exchange_segments.by_sample.tsv"
        (long[long_cols] if len(long) else long).to_csv(long_path, sep="\t", index=False)
        n_samples = long["sample_id"].nunique() if len(long) else 0
        log(f"wrote {n_samples} per-sample BEDs to {by_dir}/ and {long_path}")

    # ---- merged across dyads (support = pileup depth) ----
    merged_path = args.outdir / "exchange_segments.merged.bed"
    if len(seg):
        merged = merge_intervals(seg)
        merged.to_csv(merged_path, sep="\t", index=False, header=False)
    else:
        merged = pd.DataFrame(columns=["chrom", "bed_start", "bed_end", "support"])
        merged_path.write_text("")
    log(f"wrote {merged_path}  ({len(merged):,} merged blocks)")

    # ---- summary ----
    sum_path = args.outdir / "exchange_mask_summary.tsv"
    rows = []
    for cls in sorted(keep):
        sub = seg[seg["class"] == cls]
        rows.append({
            "group": "class", "key": cls, "n_tracts": int(len(sub)),
            "total_bp": int((sub["bed_end"] - sub["bed_start"]).sum()) if len(sub) else 0,
        })
    rows.append({
        "group": "ALL", "key": "exchange_tracts", "n_tracts": int(len(seg)),
        "total_bp": int((seg["bed_end"] - seg["bed_start"]).sum()) if len(seg) else 0,
    })
    rows.append({
        "group": "ALL", "key": "merged_blocks", "n_tracts": int(len(merged)),
        "total_bp": int((merged["bed_end"] - merged["bed_start"]).sum()) if len(merged) else 0,
    })
    if args.inversion_atlas is not None and len(seg):
        per_inv = (seg.assign(bp=seg["bed_end"] - seg["bed_start"])
                      .groupby("inversion_id")
                      .agg(n_tracts=("interval_id", "count"), total_bp=("bp", "sum"))
                      .reset_index())
        for r in per_inv.itertuples():
            rows.append({"group": "inversion", "key": r.inversion_id,
                         "n_tracts": int(r.n_tracts), "total_bp": int(r.total_bp)})
    pd.DataFrame(rows).to_csv(sum_path, sep="\t", index=False)
    log(f"wrote {sum_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

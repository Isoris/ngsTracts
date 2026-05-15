#!/usr/bin/env python3
"""
STEP_TRC_01_classify_intervals.py — main ngsTracts interval classifier

Reads ngsPedigree Stage 3 outputs (departure_intervals.tsv, optional
chromosome_background_intervals.tsv, optional inversion_atlas.tsv) and
emits per-interval classifications + per-dyad and per-chrom summaries.

Vectorized over the intervals dataframe. No per-row Python loops in the
classification logic.

See docs/METHODOLOGY.md for the decision tree.
See docs/SCHEMA.md for the input contract.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd


SUPPORTED_SCHEMA_VERSIONS = {"0.1"}
NGSTRACTS_VERSION = "0.1.0"


# --------------------------------------------------------------------------
# Defaults — keep in sync with docs/METHODOLOGY.md §7
# --------------------------------------------------------------------------
DEFAULTS = dict(
    min_n_sites_per_interval=3,
    short_tract_max_bp=50_000,
    mosaic_short_max_bp=200_000,
    co_breakpoint_max_resolution=1_000_000,
    implausible_dco_min_bp=5_000_000,
    implausible_nco_min_bp=100_000,
    discordant_frac_strong=0.9,
    n_sites_high_conf_min=5,
    chrom_end_proximity_bp=2_000_000,
    prior_max_dco_rate=0.01,
    prior_gc_constant=0.001,
)


CLASS_VALUES = ("NCO", "CO", "DCO", "MOSAIC_SHORT", "MOSAIC_LONG", "AMBIG", "LOW_CONFIDENCE")


def log(msg: str) -> None:
    sys.stderr.write(f"[{time.strftime('%F %T')}] [ngsTracts] {msg}\n")
    sys.stderr.flush()


def die(msg: str, code: int = 1) -> None:
    log(f"[ERROR] {msg}")
    sys.exit(code)


# --------------------------------------------------------------------------
# Schema check
# --------------------------------------------------------------------------
def check_schema(args_path: Path) -> dict[str, str]:
    if not args_path.exists():
        die(f"Stage 3 args file not found: {args_path}")
    rows = pd.read_csv(args_path, sep="\t", dtype=str)
    if list(rows.columns[:2]) != ["key", "value"]:
        # Tolerate no-header form
        rows = pd.read_csv(args_path, sep="\t", header=None, names=["key", "value"], dtype=str)
    args_dict = dict(zip(rows["key"], rows["value"]))
    sv = args_dict.get("schema_version", "")
    sv_short = ".".join(sv.split(".")[:2]) if sv else ""
    if sv_short not in SUPPORTED_SCHEMA_VERSIONS:
        die(
            f"Schema version mismatch.\n"
            f"  Stage 3 emitted schema_version: {sv!r}\n"
            f"  ngsTracts {NGSTRACTS_VERSION} supports: {sorted(SUPPORTED_SCHEMA_VERSIONS)}\n"
            f"  See docs/SCHEMA_COMPATIBILITY.md"
        )
    log(f"schema check OK: stage3={sv} ngstracts={NGSTRACTS_VERSION}")
    return args_dict


# --------------------------------------------------------------------------
# IO
# --------------------------------------------------------------------------
DEPARTURE_REQUIRED = [
    "interval_id", "parent_id", "offspring_id", "chrom",
    "start", "end", "n_sites", "n_discordant", "span_bp",
    "flanking_left_state", "flanking_right_state", "departure_state",
    "distance_to_nearest_inv_bp", "inside_inversion", "confidence",
]


def load_departure_intervals(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    missing = [c for c in DEPARTURE_REQUIRED if c not in df.columns]
    if missing:
        die(f"departure_intervals.tsv missing required columns: {missing}")

    for c in ("start", "end", "n_sites", "n_discordant", "span_bp"):
        df[c] = pd.to_numeric(df[c], errors="raise").astype(np.int64)

    # distance can be "-" (no inversions on chrom)
    df["distance_to_nearest_inv_bp"] = pd.to_numeric(
        df["distance_to_nearest_inv_bp"].replace("-", np.nan),
        errors="coerce",
    )

    # Constraints
    bad_span = df["span_bp"] != (df["end"] - df["start"] + 1)
    if bad_span.any():
        die(f"span_bp != end-start+1 for {int(bad_span.sum())} intervals")
    if (df["n_discordant"] > df["n_sites"]).any():
        die("n_discordant > n_sites in some intervals")
    if (df["start"] > df["end"]).any():
        die("start > end in some intervals")

    for col, allowed in [
        ("flanking_left_state",  {"hapA", "hapB", "boundary"}),
        ("flanking_right_state", {"hapA", "hapB", "boundary"}),
        ("departure_state",      {"hapA", "hapB", "neither"}),
        ("inside_inversion",     {"yes", "no", "partial"}),
        ("confidence",           {"high", "medium", "low"}),
    ]:
        bad = ~df[col].isin(allowed)
        if bad.any():
            uniq = df.loc[bad, col].unique().tolist()
            die(f"{col} has invalid values: {uniq} (allowed: {sorted(allowed)})")

    if df["interval_id"].duplicated().any():
        die("interval_id has duplicates")

    log(f"loaded {len(df):,} departure intervals from {path}")
    return df


def load_chrom_lengths(fai_path: Path) -> dict[str, int]:
    if not fai_path.exists():
        die(f"FAI not found: {fai_path}")
    df = pd.read_csv(
        fai_path, sep="\t", header=None,
        names=["chrom", "length", "offset", "linebases", "linewidth"],
        usecols=["chrom", "length"], dtype={"chrom": str, "length": np.int64},
    )
    return dict(zip(df["chrom"], df["length"]))


def load_parent_karyotypes(path: Path | None):
    """Optional. Returns {(inversion_id, sample_id): karyotype 0/1/2} or None."""
    if path is None:
        return None
    if not path.exists():
        die(f"--parent-karyotypes file not found: {path}")
    df = pd.read_csv(path, sep="\t", dtype=str)
    need = {"inversion_id", "sample_id", "karyotype"}
    if not need.issubset(df.columns):
        die(f"karyotypes file missing columns; need {need}")
    df["karyotype"] = pd.to_numeric(df["karyotype"], errors="raise").astype(np.int64)
    return {(r.inversion_id, r.sample_id): int(r.karyotype) for r in df.itertuples()}


# --------------------------------------------------------------------------
# Vectorized classifier
# --------------------------------------------------------------------------
def classify(df: pd.DataFrame, params: dict, chrom_len: dict[str, int]) -> pd.DataFrame:
    """
    Apply the decision tree in METHODOLOGY §3 vectorized over the dataframe.
    First-matching-rule semantics: each rule only updates rows not yet
    finalized (`assigned` mask).
    """
    n = len(df)
    cls = np.full(n, "AMBIG", dtype=object)
    conf = np.full(n, "medium", dtype=object)
    review = np.zeros(n, dtype=np.int8)
    assigned = np.zeros(n, dtype=bool)

    span = df["span_bp"].to_numpy()
    nsites = df["n_sites"].to_numpy()
    ndisc = df["n_discordant"].to_numpy()
    inside = df["inside_inversion"].to_numpy()
    flank_L = df["flanking_left_state"].to_numpy()
    flank_R = df["flanking_right_state"].to_numpy()
    in_conf = df["confidence"].to_numpy()

    with np.errstate(divide="ignore", invalid="ignore"):
        disc_frac = np.where(nsites > 0, ndisc / nsites, 0.0)
    flank_same = (flank_L == flank_R)
    flanks_real = (flank_L != "boundary") & (flank_R != "boundary")

    # 3.1 confidence gate / min-sites
    mask = (in_conf == "low") | (nsites < params["min_n_sites_per_interval"])
    cls[mask] = "LOW_CONFIDENCE"
    conf[mask] = "low"
    assigned |= mask

    # 3.2 inside-inversion rule
    in_inv = (inside == "yes") & ~assigned

    short = in_inv & (span < params["short_tract_max_bp"])
    cls[short] = "NCO"
    conf[short] = np.where(nsites[short] >= params["n_sites_high_conf_min"], "high", "medium")
    assigned |= short

    mid = in_inv & (span >= params["short_tract_max_bp"]) & (span < params["mosaic_short_max_bp"]) & ~assigned
    cls[mid] = "MOSAIC_SHORT"
    conf[mid] = "medium"
    review[mid] = 1
    assigned |= mid

    long_in_inv = in_inv & (span >= params["mosaic_short_max_bp"]) & ~assigned
    cls[long_in_inv] = "MOSAIC_LONG"
    conf[long_in_inv] = "low"
    review[long_in_inv] = 1
    assigned |= long_in_inv

    # 3.3 short interval, outside inversion
    short_out = (inside == "no") & (span < params["short_tract_max_bp"]) & ~assigned
    short_out_ok = short_out & flank_same
    cls[short_out_ok] = "NCO"
    strong = (disc_frac >= params["discordant_frac_strong"]) & (nsites >= params["n_sites_high_conf_min"])
    conf[short_out_ok] = np.where(strong[short_out_ok], "high", "medium")
    assigned |= short_out_ok

    short_out_bad = short_out & ~flank_same
    cls[short_out_bad] = "AMBIG"
    conf[short_out_bad] = "low"
    review[short_out_bad] = 1
    assigned |= short_out_bad

    # 3.4 long interval with flank switch -> CO
    co = ~flank_same & ~assigned & flanks_real
    cls[co] = "CO"
    conf[co & (span <= 200_000)] = "high"
    conf[co & (span > 200_000) & (span <= params["co_breakpoint_max_resolution"])] = "medium"
    conf[co & (span > params["co_breakpoint_max_resolution"])] = "low"
    assigned |= co

    # 3.5 long, matching flanks
    long_match = flank_same & (span >= params["short_tract_max_bp"]) & ~assigned
    dco = long_match & (span <= params["mosaic_short_max_bp"]) & (disc_frac >= params["discordant_frac_strong"])
    cls[dco] = "DCO"
    conf[dco] = "medium"
    assigned |= dco

    ml = long_match & ~assigned
    cls[ml] = "MOSAIC_LONG"
    conf[ml] = "low"
    review[ml] = 1
    assigned |= ml

    # 3.6 everything not assigned -> AMBIG (default), low conf, review
    leftover = ~assigned
    review[leftover] = 1
    conf[leftover] = "low"

    # ---- Manual review flags (§6) ----
    review[(cls == "DCO") & (span > params["implausible_dco_min_bp"])] = 1
    review[(cls == "DCO") & (inside == "yes")] = 1
    review[(cls == "NCO") & (span > params["implausible_nco_min_bp"])] = 1
    review[conf == "low"] = 1
    review[nsites < params["n_sites_high_conf_min"]] = 1

    # ---- Position prior (informational) ----
    chrom_arr = df["chrom"].to_numpy()
    midpoint = (df["start"].to_numpy() + df["end"].to_numpy()) / 2
    clen = np.array([chrom_len.get(c, np.nan) for c in chrom_arr], dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = midpoint / clen
    prior_dco = params["prior_max_dco_rate"] * 4 * frac * (1 - frac)
    prior_log_ratio = np.log(np.maximum(prior_dco, 1e-12) / params["prior_gc_constant"])

    out = df.copy()
    out["class"] = cls
    out["confidence_out"] = conf  # avoid clobbering input column
    out["manual_review_flag"] = review
    out["prior_log_ratio_co_nco"] = prior_log_ratio
    out["refined_breakpoint_bp"] = pd.NA
    out["refined_ci_left"] = pd.NA
    out["refined_ci_right"] = pd.NA
    out["notes"] = ""

    return out


# --------------------------------------------------------------------------
# Summaries
# --------------------------------------------------------------------------
def dyad_summary(calls: pd.DataFrame, chrom_len: dict[str, int]) -> pd.DataFrame:
    counts = (
        calls
        .groupby(["parent_id", "offspring_id", "chrom", "class"])
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )
    for c in CLASS_VALUES:
        if c not in counts.columns:
            counts[c] = 0

    counts["chrom_len_bp"] = counts["chrom"].map(chrom_len).astype("Int64")
    counts["n_intervals_total"] = counts[list(CLASS_VALUES)].sum(axis=1)
    counts["fraction_classified"] = (
        (counts["n_intervals_total"] - counts["AMBIG"] - counts["LOW_CONFIDENCE"])
        / counts["n_intervals_total"].replace(0, np.nan)
    )

    mb = counts["chrom_len_bp"].astype(float) / 1e6
    counts["co_per_mb"] = (counts["CO"] + counts["DCO"]) / mb
    counts["nco_per_mb"] = counts["NCO"] / mb
    co_total = counts["CO"] + counts["DCO"]
    counts["co_nco_ratio"] = co_total / counts["NCO"].replace(0, np.nan)

    frac_inv = (
        calls.assign(in_inv=(calls["inside_inversion"] == "yes").astype(int))
             .groupby(["parent_id", "offspring_id", "chrom"])["in_inv"]
             .mean()
             .reset_index(name="fraction_in_inversions")
    )
    counts = counts.merge(frac_inv, on=["parent_id", "offspring_id", "chrom"], how="left")

    counts = counts.rename(columns={
        "CO": "n_co", "DCO": "n_dco", "NCO": "n_nco",
        "MOSAIC_SHORT": "n_mosaic_short", "MOSAIC_LONG": "n_mosaic_long",
        "AMBIG": "n_ambig", "LOW_CONFIDENCE": "n_low_conf",
    })

    col_order = [
        "parent_id", "offspring_id", "chrom", "chrom_len_bp",
        "n_co", "n_dco", "n_nco", "n_mosaic_short", "n_mosaic_long",
        "n_ambig", "n_low_conf",
        "co_per_mb", "nco_per_mb", "co_nco_ratio",
        "n_intervals_total", "fraction_classified",
        "fraction_in_inversions",
    ]
    return counts[col_order].sort_values(["parent_id", "offspring_id", "chrom"]).reset_index(drop=True)


def chrom_summary(dyad: pd.DataFrame) -> pd.DataFrame:
    out = dyad.groupby("chrom").agg(
        chrom_len_bp=("chrom_len_bp", "first"),
        n_dyads=("parent_id", "count"),
        total_co=("n_co", "sum"),
        total_dco=("n_dco", "sum"),
        total_nco=("n_nco", "sum"),
        mean_co_per_mb=("co_per_mb", "mean"),
        sd_co_per_mb=("co_per_mb", "std"),
        mean_nco_per_mb=("nco_per_mb", "mean"),
        sd_nco_per_mb=("nco_per_mb", "std"),
    ).reset_index()
    return out


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def parse_args(argv=None):
    p = argparse.ArgumentParser(description="ngsTracts STEP_TRC_01 — classify departure intervals")
    p.add_argument("--stage3-dir", required=True, type=Path,
                   help="ngsPedigree Stage 3 output dir (contains *.tsv and stage3.args.tsv)")
    p.add_argument("--fai", required=True, type=Path,
                   help="reference .fai for chromosome lengths")
    p.add_argument("--outdir", required=True, type=Path, help="output directory")
    p.add_argument("--parent-karyotypes", type=Path, default=None,
                   help="optional: parent karyotype table for inversion-aware NCO rule")
    p.add_argument("--departure-file", type=Path, default=None,
                   help="override default departure_intervals.tsv path")
    p.add_argument("--args-file", type=Path, default=None,
                   help="override default stage3.args.tsv path")
    for k, v in DEFAULTS.items():
        p.add_argument(f"--{k.replace('_', '-')}", type=type(v), default=v)
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)

    args_file = args.args_file or (args.stage3_dir / "stage3.args.tsv")
    dep_file = args.departure_file or (args.stage3_dir / "departure_intervals.tsv")

    stage3_args = check_schema(args_file)
    intervals = load_departure_intervals(dep_file)
    chrom_len = load_chrom_lengths(args.fai)
    _ = load_parent_karyotypes(args.parent_karyotypes)

    params = {k: getattr(args, k) for k in DEFAULTS}

    log("classifying...")
    calls = classify(intervals, params, chrom_len)

    out_cols = [
        "interval_id", "parent_id", "offspring_id", "chrom", "start", "end", "span_bp",
        "class", "confidence", "flanking_left_state", "flanking_right_state",
        "departure_state", "n_sites", "n_discordant", "inside_inversion",
        "distance_to_nearest_inv_bp",
        "prior_log_ratio_co_nco",
        "refined_breakpoint_bp", "refined_ci_left", "refined_ci_right",
        "manual_review_flag", "notes",
    ]
    # Build final output: rename confidence_out -> confidence
    final = calls.drop(columns=["confidence"]).rename(columns={"confidence_out": "confidence"})
    tract_path = args.outdir / "tract_classifications.tsv"
    final[out_cols].to_csv(tract_path, sep="\t", index=False)
    log(f"wrote {tract_path}  ({len(final):,} rows)")

    log("computing dyad summaries...")
    dyad = dyad_summary(final, chrom_len)
    dyad_path = args.outdir / "dyad_event_rates.tsv"
    dyad.to_csv(dyad_path, sep="\t", index=False)
    log(f"wrote {dyad_path}  ({len(dyad):,} rows)")

    log("computing chrom summaries...")
    chrom = chrom_summary(dyad)
    chrom_path = args.outdir / "chrom_summary.tsv"
    chrom.to_csv(chrom_path, sep="\t", index=False)
    log(f"wrote {chrom_path}  ({len(chrom):,} rows)")

    # Provenance
    args_out = args.outdir / "ngstracts.args.tsv"
    with args_out.open("w") as fh:
        fh.write("key\tvalue\n")
        fh.write(f"ngstracts_version\t{NGSTRACTS_VERSION}\n")
        fh.write(f"stage3_schema_version\t{stage3_args.get('schema_version', '?')}\n")
        fh.write(f"stage3_version\t{stage3_args.get('stage3_version', '?')}\n")
        fh.write(f"departure_file\t{dep_file}\n")
        fh.write(f"fai\t{args.fai}\n")
        fh.write(f"n_intervals\t{len(intervals)}\n")
        for k, v in params.items():
            fh.write(f"param_{k}\t{v}\n")
        fh.write(f"datetime\t{time.strftime('%FT%T%z')}\n")
    log(f"wrote {args_out}")

    log("class breakdown:")
    bd = final["class"].value_counts()
    for cls in CLASS_VALUES:
        n = int(bd.get(cls, 0))
        pct = 100 * n / max(len(final), 1)
        log(f"  {cls:<15} {n:>8}  ({pct:5.1f}%)")
    n_review = int(final["manual_review_flag"].sum())
    log(f"  manual_review_flag set on {n_review:,} ({100*n_review/max(len(final),1):.1f}%)")

    return 0


if __name__ == "__main__":
    sys.exit(main())

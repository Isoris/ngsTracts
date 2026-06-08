#!/usr/bin/env python3
"""
STEP_TRC_05_transmission_drive.py — Mendelian transmission / meiotic-drive test

Plugs in at the END OF POLARIZATION. The upstream caller (ngsPedigree) owns:
  Layer 1  karyotype calling           (HOM_REF / HET / HOM_INV per sample x inversion)
  Layer 2  dyad/triad compatibility    (impossible-genotype detection)
  Layer 3  transmission polarization    (choose REF/INV orientation; call which
                                         arrangement each HET parent transmitted)
This script is the ADAPTER from there onward: it consumes the polarized
transmission table and computes per-inversion arrangement-transmission and
meiotic-drive statistics, plus an optional join to the CO/NCO/DCO tract counts.

It does NOT re-call karyotypes, re-polarize, or re-derive haplotype phase; those
are ngsPedigree's (read via the adapter contract in
docs/SCHEMA_POLARIZATION_ADAPTER.md). Dependency-free binomial test (numpy/pandas
+ math.comb only).

Input (required): --transmissions  polarized_transmissions.tsv
  inversion_id, chrom, parent_id, offspring_id,
  parent_karyotype {HOM_REF,HET,HOM_INV}, parent_role {father,mother,unknown},
  transmitted_arrangement {REF,INV,ambiguous}, informative {0,1}
Optional: --polarity (per-inversion polarization summary, passed through),
          --tract-classifications + --inversion-atlas (adds CO/NCO/DCO counts).

Output (to --outdir): transmission_drive.tsv  (one row per inversion).
"""
from __future__ import annotations

import argparse
import sys
import time
from math import comb
from pathlib import Path

import numpy as np
import pandas as pd


def log(msg: str) -> None:
    sys.stderr.write(f"[{time.strftime('%F %T')}] [ngsTracts:trc05] {msg}\n")
    sys.stderr.flush()


def die(msg: str, code: int = 1) -> None:
    log(f"[ERROR] {msg}")
    sys.exit(code)


def binom_two_sided_p(k: int, n: int, p: float = 0.5) -> float:
    """Exact two-sided binomial test p-value (method of small p-values)."""
    if n == 0:
        return float("nan")
    pmf = [comb(n, i) * (p ** i) * ((1 - p) ** (n - i)) for i in range(n + 1)]
    obs = pmf[k]
    return float(min(1.0, sum(x for x in pmf if x <= obs * (1 + 1e-9))))


def drive_block(sub: pd.DataFrame) -> dict:
    """REF/INV transmission counts + two-sided binomial drive p for a HET-parent set."""
    n_ref = int((sub["transmitted_arrangement"] == "REF").sum())
    n_inv = int((sub["transmitted_arrangement"] == "INV").sum())
    tot = n_ref + n_inv
    rate = (n_inv / tot) if tot else float("nan")
    pval = binom_two_sided_p(n_inv, tot, 0.5) if tot else float("nan")
    return {"n_ref": n_ref, "n_inv": n_inv, "tot": tot, "rate": rate, "pval": pval}


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="ngsTracts STEP_TRC_05 — transmission / meiotic-drive test")
    p.add_argument("--transmissions", required=True, type=Path,
                   help="polarized_transmissions.tsv (end-of-polarization adapter input)")
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--polarity", type=Path, default=None,
                   help="optional per-inversion polarization summary to pass through")
    p.add_argument("--tract-classifications", type=Path, default=None,
                   help="optional STEP_TRC_01 output for CO/NCO/DCO candidate counts")
    p.add_argument("--inversion-atlas", type=Path, default=None,
                   help="required with --tract-classifications to map tracts to inversions")
    args = p.parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)

    tx = pd.read_csv(args.transmissions, sep="\t", dtype=str)
    need = {"inversion_id", "chrom", "parent_id", "offspring_id",
            "parent_karyotype", "parent_role", "transmitted_arrangement", "informative"}
    miss = need - set(tx.columns)
    if miss:
        die(f"polarized_transmissions.tsv missing columns: {sorted(miss)}")
    for col, allowed in [
        ("parent_karyotype", {"HOM_REF", "HET", "HOM_INV"}),
        ("parent_role", {"father", "mother", "unknown"}),
        ("transmitted_arrangement", {"REF", "INV", "ambiguous"}),
    ]:
        bad = ~tx[col].isin(allowed)
        if bad.any():
            die(f"{col} has invalid values: {tx.loc[bad, col].unique().tolist()}")
    tx["informative"] = pd.to_numeric(tx["informative"], errors="raise").astype(int)

    # Only HET parents with an informative, unambiguous transmission carry drive signal.
    het = tx[(tx["parent_karyotype"] == "HET")
             & (tx["informative"] == 1)
             & (tx["transmitted_arrangement"].isin(["REF", "INV"]))].copy()
    log(f"{len(tx):,} transmissions; {len(het):,} informative HET-parent transmissions")

    polarity = None
    if args.polarity is not None and args.polarity.exists():
        polarity = pd.read_csv(args.polarity, sep="\t", dtype=str).set_index("inversion_id")

    tract_counts = None
    if args.tract_classifications is not None:
        if args.inversion_atlas is None or not args.inversion_atlas.exists():
            die("--tract-classifications needs --inversion-atlas to map tracts to inversions")
        tc = pd.read_csv(args.tract_classifications, sep="\t", dtype=str)
        tc["start"] = pd.to_numeric(tc["start"]).astype(np.int64)
        tc["end"] = pd.to_numeric(tc["end"]).astype(np.int64)
        atlas = pd.read_csv(args.inversion_atlas, sep="\t", dtype=str)
        atlas["start"] = pd.to_numeric(atlas["start"]).astype(np.int64)
        atlas["end"] = pd.to_numeric(atlas["end"]).astype(np.int64)
        rows_tc = []
        for inv in atlas.itertuples():
            ov = tc[(tc["chrom"] == inv.chrom) & (tc["start"] <= inv.end) & (tc["end"] >= inv.start)]
            rows_tc.append({
                "inversion_id": inv.inversion_id,
                "n_CO_candidates": int((ov["class"] == "CO").sum()),
                "n_NCO_candidates": int((ov["class"] == "NCO").sum()),
                "n_DCO_candidates": int((ov["class"] == "DCO").sum()),
            })
        tract_counts = pd.DataFrame(rows_tc).set_index("inversion_id")

    out_rows = []
    for inv_id in sorted(tx["inversion_id"].unique()):
        sub_all = tx[tx["inversion_id"] == inv_id]
        sub = het[het["inversion_id"] == inv_id]
        comb_b = drive_block(sub)
        pat = drive_block(sub[sub["parent_role"] == "father"])
        mat = drive_block(sub[sub["parent_role"] == "mother"])

        row = {
            "inversion_id": inv_id,
            "chrom": sub_all["chrom"].mode().iat[0] if len(sub_all) else "NA",
            "n_dyads_tested": int(sub_all["offspring_id"].nunique()),
            "n_informative_het_parent_transmissions": comb_b["tot"],
            "n_REF_transmitted": comb_b["n_ref"],
            "n_INV_transmitted": comb_b["n_inv"],
            "INV_transmission_rate": comb_b["rate"],
            "mendelian_drive_pvalue": comb_b["pval"],
            "paternal_INV_transmission_rate": pat["rate"],
            "paternal_drive_pvalue": pat["pval"],
            "maternal_INV_transmission_rate": mat["rate"],
            "maternal_drive_pvalue": mat["pval"],
        }
        # pass-through polarization summary (computed upstream by ngsPedigree)
        for f in ("chosen_polarity", "polarization_score_REF_INV",
                  "polarization_score_INV_REF", "n_dyad_contradictions",
                  "n_triads_tested", "n_triad_contradictions"):
            row[f] = (polarity.loc[inv_id, f]
                      if polarity is not None and inv_id in polarity.index and f in polarity.columns
                      else pd.NA)
        for f in ("n_CO_candidates", "n_NCO_candidates", "n_DCO_candidates"):
            row[f] = (int(tract_counts.loc[inv_id, f])
                      if tract_counts is not None and inv_id in tract_counts.index else pd.NA)
        out_rows.append(row)

    cols = [
        "inversion_id", "chrom",
        "chosen_polarity", "polarization_score_REF_INV", "polarization_score_INV_REF",
        "n_dyads_tested", "n_dyad_contradictions", "n_triads_tested", "n_triad_contradictions",
        "n_informative_het_parent_transmissions", "n_REF_transmitted", "n_INV_transmitted",
        "INV_transmission_rate", "mendelian_drive_pvalue",
        "paternal_INV_transmission_rate", "paternal_drive_pvalue",
        "maternal_INV_transmission_rate", "maternal_drive_pvalue",
        "n_CO_candidates", "n_NCO_candidates", "n_DCO_candidates",
    ]
    out = pd.DataFrame(out_rows)
    for c in cols:
        if c not in out.columns:
            out[c] = pd.NA
    out_path = args.outdir / "transmission_drive.tsv"
    out[cols].to_csv(out_path, sep="\t", index=False)
    log(f"wrote {out_path}  ({len(out)} inversions)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

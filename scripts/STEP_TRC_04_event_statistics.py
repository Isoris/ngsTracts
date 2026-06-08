#!/usr/bin/env python3
"""
STEP_TRC_04_event_statistics.py — cohort event-statistics report

Reads STEP_TRC_01 tract_classifications.tsv (+ optional traversal breakpoints,
inversion atlas, parent karyotypes, sample-sex table, feature annotations) and
emits a tidy long-format table of recombination / gene-conversion statistics
modelled on the pedigree-recombination literature (Smeds et al. 2016 collared
flycatcher; Halldorsson et al. 2016/2019; Korunes & Noor 2019; Navarro et al.
1997; Betran et al. 1997).

Design discipline (same as the rest of ngsTracts): compute every statistic the
available artifacts actually support; for the rest emit an explicit
status="requires:<input>" or "model:<param>" row rather than fabricating a
value. The output is therefore the full intended statistic catalog, with each
row honestly flagged as computed / model-based / needs-external-data.

This step computes NO biological diversity statistic (pi/DAF/dXY/rho) and no
allele-level gBGC test itself — those are external (popstats / unified-ancestry,
LD engines). Such rows are emitted as requires:.

Output (to --outdir):
  event_statistics.tsv            tidy: statistic, stratum, value, n, unit,
                                  status, inspired_by
  nco_length_distribution.tsv     empirical NCO tract-length quantiles
  event_statistics.provenance.tsv which optional inputs were supplied
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd


EVENT_CLASSES = ("CO", "NCO", "DCO")
log_prefix = "[ngsTracts:trc04]"


def log(msg: str) -> None:
    sys.stderr.write(f"[{time.strftime('%F %T')}] {log_prefix} {msg}\n")
    sys.stderr.flush()


def die(msg: str, code: int = 1) -> None:
    log(f"[ERROR] {msg}")
    sys.exit(code)


def med(x) -> float:
    x = np.asarray(x, dtype=float)
    x = x[~np.isnan(x)]
    return float(np.median(x)) if len(x) else float("nan")


def load_chrom_lengths(fai_path: Path) -> dict[str, int]:
    if not fai_path.exists():
        die(f"FAI not found: {fai_path}")
    df = pd.read_csv(fai_path, sep="\t", header=None,
                     names=["chrom", "length", "o", "lb", "lw"],
                     usecols=["chrom", "length"],
                     dtype={"chrom": str, "length": np.int64})
    return dict(zip(df["chrom"], df["length"]))


def load_bed(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, comment="#",
                     usecols=[0, 1, 2], names=["chrom", "start", "end"],
                     dtype={"chrom": str})
    df["start"] = pd.to_numeric(df["start"]).astype(np.int64)  # 0-based
    df["end"] = pd.to_numeric(df["end"]).astype(np.int64)
    return df


def overlaps_any(mid_chrom, mid_pos, bed: pd.DataFrame) -> np.ndarray:
    """Boolean: is each (chrom,pos) inside any BED interval (0-based half-open)."""
    out = np.zeros(len(mid_pos), dtype=bool)
    for i, (c, p) in enumerate(zip(mid_chrom, mid_pos)):
        g = bed[bed["chrom"] == c]
        if len(g) and ((g["start"] <= p) & (p < g["end"])).any():
            out[i] = True
    return out


def nearest_distance(mid_chrom, mid_pos, bed: pd.DataFrame) -> np.ndarray:
    """Distance (bp) from each point to the nearest BED interval edge; NaN if none on chrom."""
    out = np.full(len(mid_pos), np.nan)
    for i, (c, p) in enumerate(zip(mid_chrom, mid_pos)):
        g = bed[bed["chrom"] == c]
        if not len(g):
            continue
        inside = (g["start"] <= p) & (p < g["end"])
        if inside.any():
            out[i] = 0.0
            continue
        d = np.minimum((g["start"] - p).abs(), (g["end"] - 1 - p).abs())
        out[i] = float(d.min())
    return out


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="ngsTracts STEP_TRC_04 — event statistics report")
    p.add_argument("--tract-classifications", required=True, type=Path)
    p.add_argument("--fai", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--traversal-breakpoints", type=Path, default=None,
                   help="STEP_TRC_02 output; refines CO resolution to CI width")
    p.add_argument("--inversion-atlas", type=Path, default=None,
                   help="inversion_id,chrom,start,end[,...]; enables inside-bp rates")
    p.add_argument("--parent-karyotypes", type=Path, default=None,
                   help="inversion_id,sample_id,karyotype(0/1/2); enables Het/Hom split")
    p.add_argument("--sample-sex", type=Path, default=None,
                   help="sample_id,sex(male/female); enables parent-sex effect")
    p.add_argument("--annotation", action="append", default=[],
                   metavar="NAME=BED",
                   help="feature BED for density (repeatable), e.g. genes=genes.bed")
    p.add_argument("--centromeres", type=Path, default=None,
                   help="BED of centromere/low-recombination regions")
    p.add_argument("--switch-error-rate", type=float, default=None,
                   help="per-informative-site switch-error rate for false-NCO bound")
    p.add_argument("--co-resolution-threshold-bp", type=int, default=10_000,
                   help="resolution benchmark threshold (default 10kb, as Smeds)")
    args = p.parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)

    tc = pd.read_csv(args.tract_classifications, sep="\t", dtype=str)
    for c in ("class", "chrom", "start", "end", "span_bp", "parent_id",
              "offspring_id", "inside_inversion"):
        if c not in tc.columns:
            die(f"tract_classifications.tsv missing column: {c}")
    for c in ("start", "end", "span_bp"):
        tc[c] = pd.to_numeric(tc[c], errors="raise").astype(np.int64)
    if "distance_to_nearest_inv_bp" in tc.columns:
        tc["dist_bp"] = pd.to_numeric(tc["distance_to_nearest_inv_bp"], errors="coerce")
    else:
        tc["dist_bp"] = np.nan
    tc["mid"] = (tc["start"] + tc["end"]) / 2.0

    chrom_len = load_chrom_lengths(args.fai)
    genome_bp = int(sum(chrom_len.values()))
    tc["telomere_dist"] = [
        min(s - 1, chrom_len.get(c, np.nan) - e) if c in chrom_len else np.nan
        for s, e, c in zip(tc["start"], tc["end"], tc["chrom"])
    ]

    ev = tc[tc["class"].isin(EVENT_CLASSES)].copy()
    co = ev[ev["class"] == "CO"]
    nco = ev[ev["class"] == "NCO"]
    dco = ev[ev["class"] == "DCO"]

    rows = []

    def add(stat, stratum, value, n, unit, status, inspired):
        rows.append({"statistic": stat, "stratum": stratum, "value": value,
                     "n": n, "unit": unit, "status": status, "inspired_by": inspired})

    # ---- counts (Smeds; Korunes; Navarro) ----
    add("n_CO", "ALL", int(len(co)), int(len(co)), "events", "computed", "Smeds2016;deCODE2019")
    add("n_NCO", "ALL", int(len(nco)), int(len(nco)), "events", "computed", "Smeds2016;Halldorsson2016;Korunes2019")
    add("n_DCO", "ALL", int(len(dco)), int(len(dco)), "events", "computed", "Navarro1997")
    add("n_events_total", "ALL", int(len(ev)), int(len(ev)), "events", "computed", "Smeds2016")
    for cls, sub in (("NCO", nco), ("CO", co), ("DCO", dco)):
        for strat in ("yes", "no", "partial"):
            k = int((sub["inside_inversion"] == strat).sum())
            add(f"n_{cls}", f"inside_inversion={strat}", k, k, "events", "computed", "Korunes2019")

    # ---- ratios (gene-flux / inversion diagnostics) ----
    nco_co = len(nco) / len(co) if len(co) else float("nan")
    dco_co = len(dco) / len(co) if len(co) else float("nan")
    add("NCO_CO_ratio", "ALL", nco_co, int(len(ev)), "ratio", "computed", "ngsTracts(gene-flux)")
    add("DCO_CO_ratio", "ALL", dco_co, int(len(ev)), "ratio", "computed", "inversion-heterokaryotype")
    add("gene_flux_index_count", "ALL", int(len(nco) + len(dco)), int(len(nco) + len(dco)),
        "events", "computed", "Navarro1997")
    add("gene_flux_index_bp", "ALL", int(nco["span_bp"].sum() + dco["span_bp"].sum()),
        int(len(nco) + len(dco)), "bp", "computed", "Navarro1997")

    # ---- CO resolution (Smeds: median, %<10kb) ----
    co_res = co["span_bp"].astype(float).to_numpy()
    res_source = "span_bp(flanking-marker interval)"
    if args.traversal_breakpoints is not None and args.traversal_breakpoints.exists():
        tb = pd.read_csv(args.traversal_breakpoints, sep="\t", dtype=str)
        if {"interval_id", "ci_left_bp", "ci_right_bp"}.issubset(tb.columns) and "interval_id" in co.columns:
            tb["w"] = pd.to_numeric(tb["ci_right_bp"], errors="coerce") - \
                      pd.to_numeric(tb["ci_left_bp"], errors="coerce")
            wmap = dict(zip(tb["interval_id"], tb["w"]))
            refined = co["interval_id"].map(wmap).astype(float)
            co_res = refined.fillna(co["span_bp"]).to_numpy()
            res_source = "refined CI width (traversal) w/ span_bp fallback"
    thr = args.co_resolution_threshold_bp
    add("CO_resolution_median_bp", "ALL", med(co_res), int(len(co)), f"bp [{res_source}]",
        "computed", "Smeds2016;deCODE2019")
    pct = 100 * float(np.mean(co_res < thr)) if len(co_res) else float("nan")
    add(f"CO_resolution_lt{thr}bp_pct", "ALL", pct, int(len(co)), "%", "computed", "Smeds2016")

    # ---- NCO tract length (Korunes; Betran; Halldorsson) ----
    for strat, sub in (("ALL", nco),
                       ("inside_inversion=yes", nco[nco["inside_inversion"] == "yes"]),
                       ("inside_inversion=no", nco[nco["inside_inversion"] == "no"])):
        L = sub["span_bp"].astype(float).to_numpy()
        add("NCO_length_median_bp", strat, med(L), int(len(sub)), "bp", "computed", "Korunes2019;Betran1997")
    Lall = nco["span_bp"].astype(float).to_numpy()
    if len(Lall):
        add("NCO_length_mean_bp", "ALL", float(np.mean(Lall)), int(len(Lall)), "bp", "computed", "Halldorsson2016")
        add("NCO_length_geometric_p_hat", "ALL", float(1.0 / np.mean(Lall)), int(len(Lall)),
            "1/bp", "computed", "Betran1997")
    else:
        add("NCO_length_mean_bp", "ALL", float("nan"), 0, "bp", "computed", "Halldorsson2016")
        add("NCO_length_geometric_p_hat", "ALL", float("nan"), 0, "1/bp", "computed", "Betran1997")

    # ---- distance to features (Smeds spatial; Korunes/Koury breakpoint) ----
    for cls, sub in (("CO", co), ("NCO", nco), ("DCO", dco)):
        add("event_distance_to_breakpoint_median_bp", cls, med(sub["dist_bp"]),
            int(sub["dist_bp"].notna().sum()), "bp", "computed", "Korunes2019;Koury2023")
        add("event_distance_to_telomere_median_bp", cls, med(sub["telomere_dist"]),
            int(sub["telomere_dist"].notna().sum()), "bp", "computed", "Smeds2016")

    # ---- inside vs outside rate (Korunes) ----
    if args.inversion_atlas is not None and args.inversion_atlas.exists():
        atlas = load_bed(args.inversion_atlas) if False else None  # atlas is 1-based; load directly
        a = pd.read_csv(args.inversion_atlas, sep="\t", dtype=str)
        a["start"] = pd.to_numeric(a["start"]).astype(np.int64)
        a["end"] = pd.to_numeric(a["end"]).astype(np.int64)
        inside_bp = int((a["end"] - a["start"] + 1).sum())
        outside_bp = max(genome_bp - inside_bp, 1)
        for cls, sub in (("NCO", nco), ("CO", co)):
            n_in = int((sub["inside_inversion"] == "yes").sum())
            n_out = int((sub["inside_inversion"] == "no").sum())
            r_in = n_in / (inside_bp / 1e6)
            r_out = n_out / (outside_bp / 1e6)
            add(f"{cls}_rate_inside_per_Mb", "inside_inversion=yes", r_in, n_in, "events/Mb",
                "computed", "Korunes2019")
            add(f"{cls}_rate_outside_per_Mb", "inside_inversion=no", r_out, n_out, "events/Mb",
                "computed", "Korunes2019")
            ratio = r_in / r_out if r_out else float("nan")
            add(f"{cls}_rate_ratio_inside_outside", "ALL", ratio, n_in + n_out, "ratio",
                "computed", "Korunes2019")
        # NCO preserved while CO suppressed (true-inversion vs cold-region)
        co_in = int((co["inside_inversion"] == "yes").sum()) / (inside_bp / 1e6)
        co_out = int((co["inside_inversion"] == "no").sum()) / (outside_bp / 1e6)
        nco_in = int((nco["inside_inversion"] == "yes").sum()) / (inside_bp / 1e6)
        nco_out = int((nco["inside_inversion"] == "no").sum()) / (outside_bp / 1e6)
        co_supp = 1 - (co_in / co_out) if co_out else float("nan")
        nco_ret = (nco_in / nco_out) if nco_out else float("nan")
        flag = int((co_supp > 0.5) and (nco_ret > 0.5)) if not (np.isnan(co_supp) or np.isnan(nco_ret)) else -1
        add("NCO_preserved_CO_suppressed_flag", "ALL", flag, int(len(ev)), "bool(-1=undef)",
            "computed", "inversion-vs-cold-region")
    else:
        for cls in ("NCO", "CO"):
            add(f"{cls}_rate_ratio_inside_outside", "ALL", float("nan"), 0, "ratio",
                "requires:inversion-atlas", "Korunes2019")
        add("NCO_preserved_CO_suppressed_flag", "ALL", float("nan"), 0, "bool",
            "requires:inversion-atlas", "inversion-vs-cold-region")

    # ---- by karyotype: Het vs Hom (the two-generation novelty) ----
    have_kary = (args.parent_karyotypes is not None and args.parent_karyotypes.exists()
                 and args.inversion_atlas is not None and args.inversion_atlas.exists())
    if have_kary:
        a = pd.read_csv(args.inversion_atlas, sep="\t", dtype=str)
        a["start"] = pd.to_numeric(a["start"]).astype(np.int64)
        a["end"] = pd.to_numeric(a["end"]).astype(np.int64)
        kt = pd.read_csv(args.parent_karyotypes, sep="\t", dtype=str)
        kt["karyotype"] = pd.to_numeric(kt["karyotype"], errors="raise").astype(int)
        kmap = {(r.inversion_id, r.sample_id): r.karyotype for r in kt.itertuples()}

        def parent_state(row):
            hits = a[(a["chrom"] == row.chrom) & (a["start"] <= row.end) & (a["end"] >= row.start)]
            for inv in hits["inversion_id"]:
                k = kmap.get((inv, row.parent_id))
                if k is not None:
                    return "Het" if k == 1 else "Hom"
            return "NA"

        ev_in = ev[ev["inside_inversion"].isin(["yes", "partial"])].copy()
        ev_in["kstate"] = [parent_state(r) for r in ev_in.itertuples()]
        for cls in ("CO", "NCO"):
            for ks in ("Het", "Hom"):
                k = int(((ev_in["class"] == cls) & (ev_in["kstate"] == ks)).sum())
                add(f"{cls}_count_inside_{ks}", f"karyotype={ks}", k, k, "events",
                    "computed", "classical-inversion-theory;Korunes2019")
    else:
        for cls in ("CO", "NCO"):
            add(f"{cls}_count_inside_Het", "karyotype=Het", float("nan"), 0, "events",
                "requires:parent-karyotypes+inversion-atlas", "classical-inversion-theory")
            add(f"{cls}_count_inside_Hom", "karyotype=Hom", float("nan"), 0, "events",
                "requires:parent-karyotypes+inversion-atlas", "inversion-control")

    # ---- parent-sex effect (Smeds; Halldorsson) ----
    if args.sample_sex is not None and args.sample_sex.exists():
        sx = pd.read_csv(args.sample_sex, sep="\t", dtype=str)
        sx.columns = [c.lower() for c in sx.columns]
        if not {"sample_id", "sex"}.issubset(sx.columns):
            die("--sample-sex needs columns sample_id, sex")
        sexmap = {r.sample_id: str(r.sex).lower()[0] for r in sx.itertuples()}  # m/f
        ev = ev.assign(psex=ev["parent_id"].map(sexmap))
        for s, lab in (("m", "male"), ("f", "female")):
            g = ev[ev["psex"] == s]
            n_par = g["parent_id"].nunique()
            rate = (len(g) / n_par) if n_par else float("nan")
            add("event_rate_by_parent_sex", lab, rate, int(len(g)), "events/parent",
                "computed", "Smeds2016;Halldorsson2016")
        m = ev[ev["psex"] == "m"]; f = ev[ev["psex"] == "f"]
        rm = len(m) / max(m["parent_id"].nunique(), 1)
        rf = len(f) / max(f["parent_id"].nunique(), 1)
        add("event_rate_male_over_female", "ALL", (rm / rf) if rf else float("nan"),
            int(len(ev)), "ratio", "computed", "Smeds2016(male+52%)")
    else:
        add("event_rate_by_parent_sex", "ALL", float("nan"), 0, "events/parent",
            "requires:sample-sex", "Smeds2016;Halldorsson2016")

    # ---- interference proxy (Smeds: spacing between COs) ----
    gaps = []
    for _, g in co.groupby(["parent_id", "offspring_id", "chrom"]):
        if len(g) >= 2:
            mids = np.sort(g["mid"].to_numpy())
            gaps.extend(np.diff(mids).tolist())
    if gaps:
        gaps = np.array(gaps, float)
        add("CO_interference_inter_event_gap_median_bp", "ALL", float(np.median(gaps)),
            len(gaps), "bp", "computed", "Smeds2016(interference<=14Mb)")
        cv = float(np.std(gaps) / np.mean(gaps)) if np.mean(gaps) else float("nan")
        add("CO_inter_event_gap_CV", "ALL", cv, len(gaps), "ratio(<1=>positive interference)",
            "computed", "interference-literature")
    else:
        add("CO_interference_inter_event_gap_median_bp", "ALL", float("nan"), 0, "bp",
            "computed", "Smeds2016")

    # ---- NCO-CO spatial association (Smeds) ----
    nn = []
    for key, g in nco.groupby(["parent_id", "offspring_id", "chrom"]):
        cg = co[(co["parent_id"] == key[0]) & (co["offspring_id"] == key[1]) & (co["chrom"] == key[2])]
        if len(cg):
            for mid in g["mid"]:
                nn.append(float(np.min(np.abs(cg["mid"].to_numpy() - mid))))
    add("NCO_CO_spatial_assoc_nearest_median_bp", "ALL", med(nn), len(nn), "bp",
        "computed" if nn else "computed(no co-located CO)", "Smeds2016")

    # ---- feature density (Smeds: promoters/genes/repeats) ----
    if args.annotation:
        for spec in args.annotation:
            if "=" not in spec:
                die(f"--annotation must be NAME=BED, got {spec!r}")
            name, path = spec.split("=", 1)
            bed = load_bed(Path(path))
            ov = overlaps_any(ev["chrom"].to_numpy(), ev["mid"].to_numpy(), bed)
            frac = float(np.mean(ov)) if len(ov) else float("nan")
            add(f"event_density_{name}_frac_overlap", "ALL", frac, int(len(ev)), "fraction",
                "computed", "Smeds2016")
    else:
        add("event_density_feature_frac_overlap", "ALL", float("nan"), 0, "fraction",
            "requires:annotation(NAME=BED)", "Smeds2016")

    # ---- distance to centromere / cold region ----
    if args.centromeres is not None and args.centromeres.exists():
        cb = load_bed(args.centromeres)
        for cls, sub in (("CO", co), ("NCO", nco), ("DCO", dco)):
            d = nearest_distance(sub["chrom"].to_numpy(), sub["mid"].to_numpy(), cb)
            add("event_distance_to_centromere_median_bp", cls, med(d),
                int(np.sum(~np.isnan(d))), "bp", "computed", "inversion/cold-region-control")
    else:
        add("event_distance_to_centromere_median_bp", "ALL", float("nan"), 0, "bp",
            "requires:centromeres-bed", "cold-region-control")

    # ---- low-coverage safeguards ----
    if args.switch_error_rate is not None:
        eps = args.switch_error_rate
        n_inf = int(pd.to_numeric(tc.get("n_sites", pd.Series(dtype=float)),
                                  errors="coerce").fillna(0).sum())
        # crude upper bound: a false short tract needs two switch errors (leave+return)
        exp_false = n_inf * eps * eps
        add("switch_error_rate", "ALL", eps, n_inf, "per-site", "model:param", "low-coverage-safeguard")
        add("expected_false_NCO_upper_bound", "ALL", float(exp_false), n_inf, "events",
            "model:n_sites*eps^2", "low-coverage-safeguard")
    else:
        add("switch_error_rate", "ALL", float("nan"), 0, "per-site",
            "requires:switch-error-rate", "low-coverage-safeguard")
        add("expected_false_NCO_upper_bound", "ALL", float("nan"), 0, "events",
            "requires:switch-error-rate", "low-coverage-safeguard")

    # ---- external-only statistics (documented, not computed here) ----
    add("GC_bias_at_NCO", "ALL", float("nan"), 0, "fraction_toward_GC",
        "requires:external(popstats allele/ancestral)", "Smeds2016;Halldorsson2016(gBGC)")
    add("rho_inside_outside_ratio", "ALL", float("nan"), 0, "ratio",
        "requires:external(LD map: pyrho/LDhat)", "LDhat;pyrho")
    add("IBD_NCO_support_frac", "ALL", float("nan"), 0, "fraction",
        "requires:external(IBD)", "Browning2024")

    out = pd.DataFrame(rows, columns=["statistic", "stratum", "value", "n",
                                      "unit", "status", "inspired_by"])
    out_path = args.outdir / "event_statistics.tsv"
    out.to_csv(out_path, sep="\t", index=False)
    log(f"wrote {out_path}  ({len(out)} statistic rows)")

    # NCO length distribution (empirical quantiles)
    if len(Lall):
        q = np.percentile(Lall, [5, 10, 25, 50, 75, 90, 95])
        ld = pd.DataFrame({"quantile": ["p05", "p10", "p25", "p50", "p75", "p90", "p95"],
                           "nco_length_bp": q.astype(int)})
    else:
        ld = pd.DataFrame(columns=["quantile", "nco_length_bp"])
    ld.to_csv(args.outdir / "nco_length_distribution.tsv", sep="\t", index=False)

    # provenance: which optional inputs were present
    prov = [
        ("traversal_breakpoints", bool(args.traversal_breakpoints and args.traversal_breakpoints.exists())),
        ("inversion_atlas", bool(args.inversion_atlas and args.inversion_atlas.exists())),
        ("parent_karyotypes", bool(args.parent_karyotypes and args.parent_karyotypes.exists())),
        ("sample_sex", bool(args.sample_sex and args.sample_sex.exists())),
        ("annotation", bool(args.annotation)),
        ("centromeres", bool(args.centromeres and args.centromeres.exists())),
        ("switch_error_rate", args.switch_error_rate is not None),
    ]
    pd.DataFrame(prov, columns=["optional_input", "supplied"]).to_csv(
        args.outdir / "event_statistics.provenance.tsv", sep="\t", index=False)

    n_comp = int((out["status"] == "computed").sum())
    n_req = int(out["status"].str.startswith("requires:").sum())
    n_mod = int(out["status"].str.startswith("model:").sum())
    log(f"computed={n_comp}  model={n_mod}  requires-external/optional={n_req}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

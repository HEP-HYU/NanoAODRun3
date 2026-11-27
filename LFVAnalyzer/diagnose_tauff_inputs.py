#!/usr/bin/env python3
"""
Read-only diagnostic tool for Run 3 Tau Fake Factor input `plots.root` files
(and optional sibling raw `hist_*.root` files).

Checks:
  1. THStack::BuildStack negative-bin audit across all bins and MC process groups,
     plus algebraic consistency between stack.GetStack().Last() and sum(GetHists()).
  2. MC statistical uncertainty (variance) spikes broken down by (era, region, DM, group)
     and optionally by individual dataset `hist_*.root` files.
  3. Per-era Data vs GenTau vs Total MC yields, Data-GenTau pulls, Data/GenTau, and
     Data/MC ratios in `tos` (per DM and summed) and across all 4 ABCD regions.
"""

import argparse
import glob
import math
import os
import re
import sys

REGIONS = ("lss", "los", "tss", "tos")
DM_BINS = [
    ("DM 0", 0, 1),
    ("DM 1", 1, 2),
    ("DM 10", 10, 11),
    ("DM 11", 11, 12),
]
DEFAULT_YEARS = ["2022", "2022EE", "2023", "2023BPix", "2024"]
DATA_PREFIXES = (
    "hist_muon",
    "hist_singlemuon",
    "hist_egamma",
    "hist_singleelectron",
)


def fail(msg, code=1):
    """Print an [ERROR] message detectable by scripts/checkerror.py and exit."""
    print(f"[ERROR] {msg}")
    sys.exit(code)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Diagnose Run 3 TauFF input plots.root for negative bins, MC stat spikes, and per-era Data/MC normalization."
    )
    parser.add_argument(
        "-V", "--version", default="process_0812_v7_fake",
        help="Input folder prefix (with or without _fake, e.g. process_1111_v1_fake)",
    )
    parser.add_argument(
        "-C", "--channel", required=True, choices=["muon", "electron"],
        help="Analysis channel (muon or electron)",
    )
    parser.add_argument(
        "-Y", "--years", nargs="+", default=list(DEFAULT_YEARS),
        help="Years/Eras list (default: 2022 2022EE 2023 2023BPix 2024)",
    )
    parser.add_argument(
        "--postfix", default="",
        help="Optional postfix appended to figure_<year><postfix>",
    )
    parser.add_argument(
        "--hist", default="h_tau1_decayMode_S3",
        help="Histogram canvas name for all MC/Data (default: h_tau1_decayMode_S3)",
    )
    parser.add_argument(
        "--hist-gen", default="h_tau1_decayMode_gen_S3",
        help="Histogram canvas name for genuine tau MC (default: h_tau1_decayMode_gen_S3)",
    )
    parser.add_argument(
        "--err-threshold", type=float, default=10.0,
        help="Single (era, group) bin error threshold to flag variance spikes (default: 10.0)",
    )
    parser.add_argument(
        "--rel-err-threshold", type=float, default=0.15,
        help="Combined (region, DM) relative MC error threshold to flag (default: 0.15)",
    )
    parser.add_argument(
        "--no-raw-scan", action="store_true", default=False,
        help="Disable scanning sibling raw hist_*.root files even when present",
    )
    return parser.parse_args(argv)


def resolve_region_plot_file(version, region, channel, year, postfix=""):
    """Resolve plots.root path matching fake_factor_calculator / plot_tauff_independence."""
    base = version.rstrip("/")
    candidates = []
    if base.endswith("_fake"):
        candidates.append(os.path.join(f"{base}_{region}", channel, f"figure_{year}{postfix}", "plots.root"))
    else:
        candidates.append(os.path.join(f"{base}_fake_{region}", channel, f"figure_{year}{postfix}", "plots.root"))
        candidates.append(os.path.join(f"{base}_{region}", channel, f"figure_{year}{postfix}", "plots.root"))

    for path in candidates:
        if os.path.isfile(path):
            return path

    fail(
        f"Missing plots.root for region='{region}', channel='{channel}', era='{year}'. "
        f"Searched: {', '.join(candidates)}"
    )


def find_plotit_config_dir(channel):
    """Locate plotIt/configs/Run3_<channel> relative to this script or CWD."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        os.path.abspath(os.path.join(script_dir, "..", "plotIt", "configs", f"Run3_{channel}")),
        os.path.abspath(os.path.join("NanoAODRun3", "plotIt", "configs", f"Run3_{channel}")),
        os.path.abspath(os.path.join("plotIt", "configs", f"Run3_{channel}")),
    ]
    for path in candidates:
        if os.path.isdir(path):
            return path
    return None


def parse_plotit_xsecs(files_yml_path):
    """Parse active dataset cross-section values from a plotIt files YAML without PyYAML."""
    xsecs = {}
    if not files_yml_path or not os.path.isfile(files_yml_path):
        return xsecs

    current_file = None
    entry_re = re.compile(r"^['\"]?(hist_[^'\":\s]+\.root)['\"]?\s*:\s*(?:#.*)?$")
    xsec_re = re.compile(r"^\s+cross-section\s*:\s*([0-9eE.+\-]+)")

    with open(files_yml_path, "r", encoding="utf-8") as fh:
        for raw_line in fh:
            stripped = raw_line.lstrip()
            if not stripped or stripped.startswith("#"):
                continue
            m_entry = entry_re.match(raw_line)
            if m_entry:
                current_file = m_entry.group(1)
                continue
            if not raw_line[0].isspace():
                current_file = None
                continue
            if current_file:
                m_xsec = xsec_re.match(raw_line)
                if m_xsec:
                    try:
                        xsecs[current_file] = float(m_xsec.group(1))
                    except ValueError:
                        pass
    return xsecs


def parse_plotit_luminosity(config_yml_path):
    """Parse luminosity from a plotIt config_<era>.yml without PyYAML."""
    if not config_yml_path or not os.path.isfile(config_yml_path):
        return None
    lumi_re = re.compile(r"^\s+luminosity\s*:\s*([0-9eE.+\-]+)")
    with open(config_yml_path, "r", encoding="utf-8") as fh:
        for raw_line in fh:
            if raw_line.lstrip().startswith("#"):
                continue
            m = lumi_re.match(raw_line)
            if m:
                try:
                    return float(m.group(1))
                except ValueError:
                    return None
    return None


def format_group_label(hist_obj):
    """Build a readable MC group label from a TH1 in mc_stack_0."""
    name = (hist_obj.GetName() or "").strip()
    title = (hist_obj.GetTitle() or "").strip()
    if not name:
        return title or "<unnamed>"
    if title and title != name and not title.startswith(";"):
        return f"{name} ({title})"
    return name


def extract_all_inputs(ROOT, years, version, channel, postfix, hist_name, hist_gen_name):
    """
    Read all (year, region) plots.root files and extract:
      - group_yields[(year, region, kind, dm)][group_label] = (val, err)
      - ndata[(year, region, dm)] = float
      - negative_bins: list of dicts for any MC group bin < 0
      - neg_canvases: set of (year, region, hname) with >=1 negative bin
      - stack_check: {'max_val_diff': float, 'max_err_diff': float, 'worst_cell': str}
      - plot_paths[(year, region)] = plots_root_path
    """
    group_yields = {}
    ndata = {}
    negative_bins = []
    neg_canvases = set()
    plot_paths = {}
    max_val_diff = 0.0
    max_err_diff = 0.0
    worst_cell = "N/A"

    for year in years:
        for region in REGIONS:
            path = resolve_region_plot_file(version, region, channel, year, postfix=postfix)
            plot_paths[(year, region)] = path
            tf = ROOT.TFile.Open(path, "READ")
            if not tf or tf.IsZombie():
                fail(f"Failed to open ROOT file: '{path}'")

            try:
                for kind, hname in (("all", hist_name), ("gen", hist_gen_name)):
                    canvas = tf.Get(hname)
                    if not canvas:
                        fail(f"Canvas '{hname}' not found in '{path}'")
                    pad_hi = canvas.GetPrimitive("pad_hi")
                    if not pad_hi:
                        fail(f"Primitive 'pad_hi' not found in canvas '{hname}' of '{path}'")
                    stack = pad_hi.GetPrimitive("mc_stack_0")
                    if not stack:
                        fail(f"MC stack 'mc_stack_0' not found in canvas '{hname}' of '{path}'")
                    hists = stack.GetHists()
                    if not hists or hists.GetSize() == 0:
                        fail(f"MC stack 'mc_stack_0' has no constituent histograms in '{hname}' of '{path}'")

                    for _, dm, _ in DM_BINS:
                        group_yields[(year, region, kind, dm)] = {}

                    for hist_obj in hists:
                        label = format_group_label(hist_obj)
                        nbins = hist_obj.GetNbinsX()
                        for b in range(1, nbins + 1):
                            bval = float(hist_obj.GetBinContent(b))
                            if bval < 0.0:
                                berr = float(hist_obj.GetBinError(b))
                                bcenter = float(hist_obj.GetBinCenter(b))
                                neg_canvases.add((year, region, hname))
                                negative_bins.append({
                                    "year": year,
                                    "region": region,
                                    "hist": hname,
                                    "kind": kind,
                                    "group": label,
                                    "bin_idx": b,
                                    "bin_center": bcenter,
                                    "val": bval,
                                    "err": berr,
                                })

                        for _, dm, bin_idx in DM_BINS:
                            val = float(hist_obj.GetBinContent(bin_idx))
                            err = float(hist_obj.GetBinError(bin_idx))
                            if label in group_yields[(year, region, kind, dm)]:
                                prev_v, prev_e = group_yields[(year, region, kind, dm)][label]
                                val += prev_v
                                err = math.hypot(prev_e, err)
                            group_yields[(year, region, kind, dm)][label] = (val, err)

                    built_stack = stack.GetStack()
                    if not built_stack or not built_stack.Last():
                        fail(f"stack.GetStack().Last() returned null in '{hname}' of '{path}'")
                    h_last = built_stack.Last()

                    for _, dm, bin_idx in DM_BINS:
                        s_val = float(h_last.GetBinContent(bin_idx))
                        s_err = float(h_last.GetBinError(bin_idx))
                        gmap = group_yields[(year, region, kind, dm)]
                        sum_val = sum(v for v, _ in gmap.values())
                        quad_err = math.sqrt(sum(e * e for _, e in gmap.values()))
                        vdiff = abs(s_val - sum_val)
                        ediff = abs(s_err - quad_err)
                        if vdiff > max_val_diff or ediff > max_err_diff:
                            max_val_diff = max(max_val_diff, vdiff)
                            max_err_diff = max(max_err_diff, ediff)
                            worst_cell = f"{year}/{region}/{hname}/DM{dm}"

                    if kind == "all":
                        h_data = None
                        for prim in pad_hi.GetListOfPrimitives():
                            if prim.InheritsFrom("TH1"):
                                h_data = prim
                        if h_data is None:
                            fail(f"Data TH1 primitive not found in '{hname}' of '{path}'")
                        for _, dm, bin_idx in DM_BINS:
                            ndata[(year, region, dm)] = float(h_data.GetBinContent(bin_idx))
            finally:
                tf.Close()

    stack_check = {
        "max_val_diff": max_val_diff,
        "max_err_diff": max_err_diff,
        "worst_cell": worst_cell,
    }
    return group_yields, ndata, negative_bins, neg_canvases, stack_check, plot_paths


def compute_region_dm_stats(years, group_yields, ndata, err_threshold, rel_err_threshold):
    """
    Aggregate yields and errors across eras/groups for each (region, DM), rank single
    (year, region, DM, group) variance spikes, and identify flagged (region, DM) cells.
    """
    cell_totals = {}
    group_entries = []

    for region in REGIONS:
        for dm_label, dm, _ in DM_BINS:
            mc_val = 0.0
            mc_var = 0.0
            gen_val = 0.0
            gen_var = 0.0
            d_val = 0.0
            for year in years:
                d_val += ndata[(year, region, dm)]
                for v, e in group_yields[(year, region, "all", dm)].values():
                    mc_val += v
                    mc_var += e * e
                for v, e in group_yields[(year, region, "gen", dm)].values():
                    gen_val += v
                    gen_var += e * e

            mc_err = math.sqrt(mc_var)
            gen_err = math.sqrt(gen_var)
            rel_err = (mc_err / abs(mc_val)) if mc_val != 0.0 else float("inf")
            cell_totals[(region, dm)] = {
                "dm_label": dm_label,
                "mc_val": mc_val,
                "mc_err": mc_err,
                "mc_var": mc_var,
                "rel_err": rel_err,
                "gen_val": gen_val,
                "gen_err": gen_err,
                "data": d_val,
            }

            for year in years:
                for group_label, (v, e) in group_yields[(year, region, "all", dm)].items():
                    var_share = (100.0 * (e * e) / mc_var) if mc_var > 0.0 else 0.0
                    group_entries.append({
                        "year": year,
                        "region": region,
                        "dm": dm,
                        "dm_label": dm_label,
                        "group": group_label,
                        "val": v,
                        "err": e,
                        "var_share": var_share,
                    })

    group_entries.sort(key=lambda r: (r["err"], abs(r["val"])), reverse=True)

    flagged_cells = []
    for region in REGIONS:
        for _, dm, _ in DM_BINS:
            info = cell_totals[(region, dm)]
            max_single_err = max(
                (e["err"] for e in group_entries if e["region"] == region and e["dm"] == dm),
                default=0.0,
            )
            if info["rel_err"] > rel_err_threshold or max_single_err > err_threshold:
                flagged_cells.append((region, dm))

    return cell_totals, group_entries, flagged_cells


def find_raw_era_dir(plots_root_path, year, postfix):
    """Look for a sibling raw era directory containing hist_*.root next to figure_<year>."""
    chan_dir = os.path.dirname(os.path.dirname(os.path.abspath(plots_root_path)))
    candidates = [
        os.path.join(chan_dir, f"{year}_postprocess{postfix}"),
        os.path.join(chan_dir, f"{year}_postprocess"),
        os.path.join(chan_dir, f"{year}{postfix}"),
        os.path.join(chan_dir, f"{year}"),
    ]
    seen = set()
    for cand in candidates:
        if cand in seen:
            continue
        seen.add(cand)
        if os.path.isdir(cand) and glob.glob(os.path.join(cand, "hist_*.root")):
            return cand
    return None


def is_data_root_file(filename):
    low = filename.lower()
    return any(low.startswith(prefix) for prefix in DATA_PREFIXES)


def scan_raw_datasets_for_spike(ROOT, era_dir, year, channel, hist_name, dm_bin_idx, config_dir):
    """Inspect individual MC hist_*.root files in era_dir for a specific DM bin."""
    files_yml = None
    lumi = None
    xsecs = {}
    if config_dir:
        files_yml = os.path.join(config_dir, "files24.yml" if year == "2024" else "files.yml")
        xsecs = parse_plotit_xsecs(files_yml)
        lumi = parse_plotit_luminosity(os.path.join(config_dir, f"config_{year}.yml"))

    records = []
    for path in sorted(glob.glob(os.path.join(era_dir, "hist_*.root"))):
        fname = os.path.basename(path)
        if is_data_root_file(fname):
            continue
        tf = ROOT.TFile.Open(path, "READ")
        if not tf or tf.IsZombie():
            records.append({
                "file": fname, "error": "failed to open ROOT file",
                "val": 0.0, "err": 0.0, "n_eff": 0.0,
                "val_scaled": None, "err_scaled": None, "nevents": None,
            })
            continue
        try:
            h = tf.Get(hist_name)
            if not h:
                continue
            val = float(h.GetBinContent(dm_bin_idx))
            err = float(h.GetBinError(dm_bin_idx))
            n_eff = (val / err) ** 2 if err > 0.0 else 0.0

            nevents = None
            hcounter = tf.Get("hcounter")
            if hcounter:
                nevents = float(hcounter.GetBinContent(2))
            else:
                h_nev = tf.Get("h_nevents")
                if h_nev:
                    nevents = float(h_nev.Integral())

            val_scaled = None
            err_scaled = None
            xsec = xsecs.get(fname)
            if xsec is not None and lumi is not None and nevents and nevents > 0.0:
                scale = xsec * lumi / nevents
                val_scaled = val * scale
                err_scaled = err * scale

            records.append({
                "file": fname,
                "val": val,
                "err": err,
                "n_eff": n_eff,
                "nevents": nevents,
                "xsec": xsec,
                "val_scaled": val_scaled,
                "err_scaled": err_scaled,
            })
        finally:
            tf.Close()

    has_scaled = any(r["err_scaled"] is not None for r in records)
    if has_scaled:
        records.sort(key=lambda r: (r["err_scaled"] if r["err_scaled"] is not None else -1.0, r["err"]), reverse=True)
    else:
        records.sort(key=lambda r: (r["err"], abs(r["val"])), reverse=True)
    neg_records = [r for r in records if r["val"] < 0.0]
    return records, neg_records, lumi is not None and bool(xsecs)


def build_negative_bin_section(negative_bins, neg_canvases):
    lines = [
        "=" * 100,
        "[Section 1] THStack::BuildStack Negative Bin Audit",
        "=" * 100,
        f"Canvases with >=1 MC group having bin content < 0: {len(neg_canvases)} "
        f"(total negative group-bin entries: {len(negative_bins)})",
    ]
    if not negative_bins:
        lines.append("No negative MC group bins detected in any histogram.")
        lines.append("")
        return "\n".join(lines)

    ordered = sorted(negative_bins, key=lambda r: r["val"])
    lines.append("")
    lines.append(
        f"{'Era':<8} | {'Reg':<4} | {'Hist':<24} | {'MC Group':<22} | {'Bin(center)':<11} | {'Yield +/- Err':<22}"
    )
    lines.append("-" * 100)
    for row in ordered:
        bin_str = f"{row['bin_idx']} ({row['bin_center']:.1f})"
        ye_str = f"{row['val']:+.3f} +/- {row['err']:.3f}"
        lines.append(
            f"{row['year']:<8} | {row['region']:<4} | {row['hist']:<24} | "
            f"{row['group'][:22]:<22} | {bin_str:<11} | {ye_str:<22}"
        )
    lines.append("")
    return "\n".join(lines)


def build_variance_spike_section(
    years, group_yields, cell_totals, group_entries, flagged_cells, err_threshold, rel_err_threshold
):
    lines = [
        "=" * 100,
        "[Section 2] MC Statistical Error Spike Breakdown",
        "=" * 100,
        f"Top 15 single (Era, Region, DM, MC Group) error contributors (flag thresholds: "
        f"single err > {err_threshold:g} or combined rel err > {100.0 * rel_err_threshold:.1f}%):",
        "",
        f"{'Rank':<4} | {'Era':<8} | {'Reg':<4} | {'DM':<6} | {'MC Group':<22} | {'Yield +/- Err':<22} | {'Var Share':>9}",
        "-" * 100,
    ]
    for idx, row in enumerate(group_entries[:15], start=1):
        ye_str = f"{row['val']:.2f} +/- {row['err']:.2f}"
        lines.append(
            f"{idx:<4} | {row['year']:<8} | {row['region']:<4} | {row['dm_label']:<6} | "
            f"{row['group'][:22]:<22} | {ye_str:<22} | {row['var_share']:8.2f}%"
        )

    if not flagged_cells:
        lines.append("")
        lines.append("No (region, DM) cells exceeded the MC error thresholds.")
        lines.append("")
        return "\n".join(lines)

    lines.append("")
    lines.append(f"Detailed Era x MC Group breakdown for {len(flagged_cells)} flagged (Region, DM) cell(s):")
    for region, dm in flagged_cells:
        info = cell_totals[(region, dm)]
        lines.append("")
        lines.append(
            f"--- Flagged Cell: region={region}, {info['dm_label']} | "
            f"Combined MC = {info['mc_val']:.2f} +/- {info['mc_err']:.2f} "
            f"({100.0 * info['rel_err']:.2f}% rel err), Data = {info['data']:.0f}, "
            f"GenTau = {info['gen_val']:.2f} +/- {info['gen_err']:.2f} ---"
        )
        lines.append(f"{'Era':<8} | {'MC Group':<24} | {'Yield +/- Err':<24} | {'Var Share':>9} | {'Note'}")
        lines.append("-" * 88)
        sub_rows = []
        for year in years:
            for grp, (v, e) in group_yields[(year, region, "all", dm)].items():
                share = (100.0 * (e * e) / info["mc_var"]) if info["mc_var"] > 0.0 else 0.0
                sub_rows.append((e, share, year, grp, v))
        sub_rows.sort(key=lambda item: (item[0], abs(item[4])), reverse=True)
        for e, share, year, grp, v in sub_rows:
            if e <= 0.0 and v == 0.0:
                continue
            note = "SPIKE" if (e > err_threshold or share >= 10.0) else ("NEG" if v < 0.0 else "")
            if share >= 0.5 or e >= 1.0 or v < 0.0:
                lines.append(
                    f"{year:<8} | {grp[:24]:<24} | {f'{v:.2f} +/- {e:.2f}':<24} | {share:8.2f}% | {note}"
                )
    lines.append("")
    return "\n".join(lines)


def build_raw_drilldown_section(ROOT, years, plot_paths, group_entries, flagged_cells, channel, hist_name):
    """Run raw hist_*.root drill-down for top unique (year, region, dm) spikes."""
    flagged_set = set(flagged_cells)
    top_triples = []
    seen = set()
    for row in group_entries:
        if flagged_set and (row["region"], row["dm"]) not in flagged_set:
            continue
        triple = (row["year"], row["region"], row["dm"])
        if triple not in seen:
            seen.add(triple)
            top_triples.append(triple)
        if len(top_triples) >= 5:
            break

    if not top_triples:
        return ""

    config_dir = find_plotit_config_dir(channel)
    dm_to_bin = {dm: (dm_lbl, bidx) for dm_lbl, dm, bidx in DM_BINS}
    lines = []
    found_any_raw = False

    for year, region, dm in top_triples:
        era_dir = find_raw_era_dir(plot_paths[(year, region)], year, "")
        if not era_dir:
            continue
        found_any_raw = True
        dm_lbl, bin_idx = dm_to_bin[dm]
        records, neg_records, has_norm = scan_raw_datasets_for_spike(
            ROOT, era_dir, year, channel, hist_name, bin_idx, config_dir
        )
        lines.append(f"Raw dataset drill-down: era={year}, region={region}, {dm_lbl} (dir: {era_dir})")
        if has_norm:
            lines.append("  [Note] Scaled yields use cross-section * lumi / hcounter[bin=2] from plotIt configs.")
            lines.append(
                f"  {'Dataset':<42} | {'Est. Scaled (Y +/- Err)':<24} | {'Raw (Y +/- Err)':<22} | {'N_eff':>8}"
            )
            lines.append("  " + "-" * 104)
            show = list(records[:10])
            shown_names = {r["file"] for r in show}
            for nr in neg_records:
                if nr["file"] not in shown_names:
                    show.append(nr)
                    shown_names.add(nr["file"])
            for r in show:
                sc_str = (
                    f"{r['val_scaled']:.2f} +/- {r['err_scaled']:.2f}"
                    if r["val_scaled"] is not None else "N/A (missing xsec/cnt)"
                )
                raw_str = f"{r['val']:.2f} +/- {r['err']:.2f}"
                lines.append(f"  {r['file'][:42]:<42} | {sc_str:<24} | {raw_str:<22} | {r['n_eff']:8.1f}")
        else:
            lines.append("  [Note] plotIt xsec/lumi configs not found; showing raw unscaled bin values and N_eff.")
            lines.append(f"  {'Dataset':<42} | {'Raw Yield +/- Err':<24} | {'N_eff':>8}")
            lines.append("  " + "-" * 80)
            for r in records[:10]:
                raw_str = f"{r['val']:.2f} +/- {r['err']:.2f}"
                lines.append(f"  {r['file'][:42]:<42} | {raw_str:<24} | {r['n_eff']:8.1f}")
        lines.append("")

    if not found_any_raw:
        return "Raw hist_*.root drill-down skipped: no sibling raw era directory found next to figure_<year>.\n"
    return "\n".join(lines)


def aggregate_era_region(years, group_yields, ndata, region, dm_list):
    """Compute per-era and combined Data, GenTau, TotalMC, pull, Data/GenTau, Data/MC."""
    rows = []
    tot_d = 0.0
    tot_g = 0.0
    tot_g_var = 0.0
    tot_m = 0.0
    tot_m_var = 0.0

    for year in years:
        d = sum(ndata[(year, region, dm)] for dm in dm_list)
        g = 0.0
        g_var = 0.0
        m = 0.0
        m_var = 0.0
        for dm in dm_list:
            for v, e in group_yields[(year, region, "gen", dm)].values():
                g += v
                g_var += e * e
            for v, e in group_yields[(year, region, "all", dm)].values():
                m += v
                m_var += e * e
        g_err = math.sqrt(g_var)
        m_err = math.sqrt(m_var)
        tot_d += d
        tot_g += g
        tot_g_var += g_var
        tot_m += m
        tot_m_var += m_var
        rows.append(_make_norm_row(year, d, g, g_err, m, m_err))

    rows.append(_make_norm_row("Combined", tot_d, tot_g, math.sqrt(tot_g_var), tot_m, math.sqrt(tot_m_var)))
    return rows


def _make_norm_row(era_label, d, g, g_err, m, m_err):
    diff = d - g
    pull = (diff / g_err) if g_err > 0.0 else float("nan")
    d_over_g = (d / g) if g != 0.0 else float("nan")
    d_over_m = (d / m) if m != 0.0 else float("nan")
    return {
        "era": era_label,
        "data": d,
        "gen": g,
        "gen_err": g_err,
        "mc": m,
        "mc_err": m_err,
        "diff": diff,
        "pull": pull,
        "d_over_g": d_over_g,
        "d_over_m": d_over_m,
    }


def build_per_era_norm_section(years, group_yields, ndata):
    lines = [
        "=" * 100,
        "[Section 3] Per-Era Data vs GenTau vs Total MC Breakdown",
        "=" * 100,
    ]
    dm_specs = [(lbl, [dm]) for lbl, dm, _ in DM_BINS] + [("Sum(DM 0,1,10,11)", [0, 1, 10, 11])]

    tos_sum_rows = None
    for title, dm_list in dm_specs:
        rows = aggregate_era_region(years, group_yields, ndata, "tos", dm_list)
        if title.startswith("Sum"):
            tos_sum_rows = rows
        lines.append(f"Region = tos | {title}")
        lines.append(
            f"{'Era':<9} | {'Data':>8} | {'GenTau +/- Err':<20} | {'TotalMC +/- Err':<20} | "
            f"{'Data-GenTau (pull)':<20} | {'Data/Gen':>8} | {'Data/MC':>8}"
        )
        lines.append("-" * 105)
        for r in rows:
            if r["era"] == "Combined":
                lines.append("-" * 105)
            gen_s = f"{r['gen']:.2f} +/- {r['gen_err']:.2f}"
            mc_s = f"{r['mc']:.2f} +/- {r['mc_err']:.2f}"
            pull_s = f"{r['diff']:+.2f} ({r['pull']:+.1f}s)" if not math.isnan(r["pull"]) else f"{r['diff']:+.2f} (N/A)"
            dg_s = f"{r['d_over_g']:.3f}" if not math.isnan(r["d_over_g"]) else "N/A"
            dm_s = f"{r['d_over_m']:.3f}" if not math.isnan(r["d_over_m"]) else "N/A"
            lines.append(
                f"{r['era']:<9} | {r['data']:8.0f} | {gen_s:<20} | {mc_s:<20} | "
                f"{pull_s:<20} | {dg_s:>8} | {dm_s:>8}"
            )
        lines.append("")

    lines.append("Compact Per-Era Sum(DM 0,1,10,11) Summary Across All 4 ABCD Regions:")
    lines.append(
        f"{'Region':<7} | {'Era':<9} | {'Data':>8} | {'GenTau +/- Err':<20} | {'TotalMC +/- Err':<20} | {'Data/Gen':>8} | {'Data/MC':>8}"
    )
    lines.append("-" * 92)
    for region in REGIONS:
        rows = aggregate_era_region(years, group_yields, ndata, region, [0, 1, 10, 11])
        for r in rows:
            gen_s = f"{r['gen']:.2f} +/- {r['gen_err']:.2f}"
            mc_s = f"{r['mc']:.2f} +/- {r['mc_err']:.2f}"
            dg_s = f"{r['d_over_g']:.3f}" if not math.isnan(r["d_over_g"]) else "N/A"
            dm_s = f"{r['d_over_m']:.3f}" if not math.isnan(r["d_over_m"]) else "N/A"
            lines.append(
                f"{region:<7} | {r['era']:<9} | {r['data']:8.0f} | {gen_s:<20} | {mc_s:<20} | {dg_s:>8} | {dm_s:>8}"
            )
        lines.append("-" * 92)
    lines.append("")
    return "\n".join(lines), tos_sum_rows


def build_takeaways_section(neg_canvases, negative_bins, stack_check, group_entries, tos_sum_rows):
    lines = [
        "=" * 100,
        "[Section 4] Automated Diagnostic Takeaways",
        "=" * 100,
    ]
    lines.append(
        f"1. Negative Bins & THStack Sum Integrity: {len(negative_bins)} negative group-bin entries "
        f"across {len(neg_canvases)} canvases. Max |GetStack().Last() - sum(GetHists())| = "
        f"{stack_check['max_val_diff']:.3e} (yield), {stack_check['max_err_diff']:.3e} (err) at {stack_check['worst_cell']}."
    )
    if stack_check["max_val_diff"] < 1e-6 and stack_check["max_err_diff"] < 1e-6:
        lines.append(
            "   -> Confirmed: THStack::BuildStack negative-bin warnings are purely visual and do NOT corrupt algebraic MC integrals."
        )
    else:
        lines.append(
            "   -> [WARNING] Non-negligible difference observed between GetStack().Last() and direct group sum!"
        )

    if group_entries:
        top = group_entries[0]
        lines.append(
            f"2. Top MC Variance Outlier: era={top['year']}, region={top['region']}, {top['dm_label']}, "
            f"group='{top['group']}' with yield = {top['val']:.2f} +/- {top['err']:.2f} "
            f"({top['var_share']:.2f}% of total combined {top['region']} {top['dm_label']} MC variance)."
        )

    era_rows = [r for r in (tos_sum_rows or []) if r["era"] != "Combined" and not math.isnan(r["d_over_m"])]
    if era_rows:
        min_dm_row = min(era_rows, key=lambda r: r["d_over_m"])
        max_dm_row = max(era_rows, key=lambda r: r["d_over_m"])
        min_dg_row = min(era_rows, key=lambda r: r["d_over_g"])
        max_dg_row = max(era_rows, key=lambda r: r["d_over_g"])
        lines.append(
            f"3. Per-Era tos Normalization Spread (Sum DM 0,1,10,11): "
            f"Data/MC spans [{min_dm_row['d_over_m']:.3f} ({min_dm_row['era']}) .. {max_dm_row['d_over_m']:.3f} ({max_dm_row['era']})]; "
            f"Data/GenTau spans [{min_dg_row['d_over_g']:.3f} ({min_dg_row['era']}) .. {max_dg_row['d_over_g']:.3f} ({max_dg_row['era']})]."
        )
        if max_dm_row["d_over_m"] - min_dm_row["d_over_m"] > 0.15:
            lines.append(
                f"   -> Large era-to-era spread (>0.15): suspect a missing/incomplete data file or lumi mismatch in '{min_dm_row['era']}'."
            )
        elif max_dm_row["d_over_m"] < 0.92:
            lines.append(
                "   -> Uniformly low Data/MC across all eras: suspect a channel-wide normalization/SF offset (e.g. Electron RECO/ID/Trig SF or Tau VSe Tight SF)."
            )
    lines.append("")
    return "\n".join(lines)


def main(argv=None):
    args = parse_args(argv)

    import ROOT
    ROOT.gROOT.SetBatch(True)

    group_yields, ndata, negative_bins, neg_canvases, stack_check, plot_paths = extract_all_inputs(
        ROOT, args.years, args.version, args.channel, args.postfix, args.hist, args.hist_gen
    )
    cell_totals, group_entries, flagged_cells = compute_region_dm_stats(
        args.years, group_yields, ndata, args.err_threshold, args.rel_err_threshold
    )

    print("=" * 100)
    print(f"Run 3 TauFF Input Diagnostics | channel={args.channel} | version={args.version}")
    print(f"Eras: {', '.join(args.years)} | postfix='{args.postfix}' | hist={args.hist} | hist_gen={args.hist_gen}")
    print(
        f"Stack Consistency Check (max |GetStack.Last - sum(GetHists)|): "
        f"yield_diff={stack_check['max_val_diff']:.3e}, err_diff={stack_check['max_err_diff']:.3e} "
        f"(at {stack_check['worst_cell']})"
    )
    print("")
    print(build_negative_bin_section(negative_bins, neg_canvases))
    print(
        build_variance_spike_section(
            args.years, group_yields, cell_totals, group_entries,
            flagged_cells, args.err_threshold, args.rel_err_threshold,
        )
    )
    if not args.no_raw_scan:
        raw_section = build_raw_drilldown_section(
            ROOT, args.years, plot_paths, group_entries, flagged_cells, args.channel, args.hist
        )
        if raw_section:
            print(raw_section)
    norm_section, tos_sum_rows = build_per_era_norm_section(args.years, group_yields, ndata)
    print(norm_section)
    print(build_takeaways_section(neg_canvases, negative_bins, stack_check, group_entries, tos_sum_rows))


if __name__ == "__main__":
    main()

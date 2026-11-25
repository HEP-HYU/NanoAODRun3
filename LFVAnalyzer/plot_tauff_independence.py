#!/usr/bin/env python3
"""
Verify ABCD variable independence (DeepTau VSjet WP vs lepton-tau charge sign)
for the Run 3 Tau Fake Factor measurement using existing lss/los/tss/tos plots.root.

Generates:
  1. MisID probability vs tau pT: OS MC (R_OS = TOS/LOS) vs SS MC (R_SS = TSS/LSS) + R_OS/R_SS ratio
  2. MisID probability vs tau DM (0, 1, 10, 11): OS MC vs SS MC + R_OS/R_SS ratio
  3. SS closure in tau pT: SS Data vs SS MC + Data/MC ratio
  4. SS closure in tau DM (0, 1, 10, 11): SS Data vs SS MC + Data/MC ratio
  5. Text summary table with Pearson 2x2 phi correlation, R_OS/R_SS, and non-closure Difference.
"""

import argparse
import ctypes
import math
import os
import sys

REGIONS = ('lss', 'los', 'tss', 'tos')
DM_BINS = [
    ('DM 0 (1p)', 1),
    ('DM 1 (1p+#pi^{0})', 2),
    ('DM 10 (3p)', 11),
    ('DM 11 (3p+#pi^{0})', 12),
]
DM_SHORT_NAMES = ['Inclusive', 'DM 0', 'DM 1', 'DM 10', 'DM 11']


def fail(msg, code=1):
    """Print an [ERROR] message detectable by scripts/checkerror.py and exit."""
    print(f"[ERROR] {msg}")
    sys.exit(code)


def resolve_region_plot_file(version, region, channel, year, postfix=''):
    """Resolve plots.root for a given fake region, channel, and era."""
    base = version.rstrip('/')
    candidates = []
    if base.endswith('_fake'):
        candidates.append(os.path.join(f'{base}_{region}', channel, f'figure_{year}{postfix}', 'plots.root'))
    else:
        candidates.append(os.path.join(f'{base}_fake_{region}', channel, f'figure_{year}{postfix}', 'plots.root'))
        candidates.append(os.path.join(f'{base}_{region}', channel, f'figure_{year}{postfix}', 'plots.root'))

    for path in candidates:
        if os.path.isfile(path):
            return path

    fail(
        f"Missing plots.root for region='{region}', channel='{channel}', era='{year}'. "
        f"Searched: {', '.join(candidates)}"
    )


def compute_abcd_metrics(A, dA, B, dB, C, dC, D, dD):
    """
    Compute ABCD independence and closure metrics from disjoint yields:
      A = LSS (Loose, SS), B = TSS (Tight, SS),
      C = LOS (Loose, OS), D = TOS (Tight, OS).
    Returns a dictionary with rates, transfer ratio K = R_OS/R_SS,
    closure Difference = 1 - R_SS/R_OS, and 2x2 Pearson phi correlation.
    """
    out = {
        'A': A, 'dA': dA, 'B': B, 'dB': dB,
        'C': C, 'dC': dC, 'D': D, 'dD': dD,
        'r_ss': None, 'dr_ss': None,
        'r_os': None, 'dr_os': None,
        'k_os_ss': None, 'dk_os_ss': None,
        'diff': None, 'ddiff': None,
        'phi': None, 'dphi': None,
    }

    if A > 0.0 and B >= 0.0:
        r_ss = B / A
        dr_ss = abs(r_ss) * math.hypot(dB / B if B > 0 else 0.0, dA / A) if (B > 0 or dB > 0) else 0.0
        if B == 0.0:
            dr_ss = dB / A
        out['r_ss'], out['dr_ss'] = r_ss, dr_ss

    if C > 0.0 and D >= 0.0:
        r_os = D / C
        dr_os = abs(r_os) * math.hypot(dD / D if D > 0 else 0.0, dC / C) if (D > 0 or dD > 0) else 0.0
        if D == 0.0:
            dr_os = dD / C
        out['r_os'], out['dr_os'] = r_os, dr_os

    if A > 0.0 and B > 0.0 and C > 0.0 and D > 0.0:
        k = (A * D) / (B * C)
        rel_k = math.sqrt((dA / A) ** 2 + (dB / B) ** 2 + (dC / C) ** 2 + (dD / D) ** 2)
        out['k_os_ss'] = k
        out['dk_os_ss'] = abs(k) * rel_k

        inv_k = 1.0 / k
        out['diff'] = 1.0 - inv_k
        out['ddiff'] = abs(inv_k) * rel_k

    n1 = A + B
    n2 = C + D
    m1 = A + C
    m2 = B + D
    if n1 > 0.0 and n2 > 0.0 and m1 > 0.0 and m2 > 0.0:
        den = math.sqrt(n1 * n2 * m1 * m2)
        phi = (A * D - B * C) / den
        dphi_dA = (D / den) - 0.5 * phi * (1.0 / n1 + 1.0 / m1)
        dphi_dB = (-C / den) - 0.5 * phi * (1.0 / n1 + 1.0 / m2)
        dphi_dC = (-B / den) - 0.5 * phi * (1.0 / n2 + 1.0 / m1)
        dphi_dD = (A / den) - 0.5 * phi * (1.0 / n2 + 1.0 / m2)
        dphi = math.sqrt(
            (dphi_dA * dA) ** 2
            + (dphi_dB * dB) ** 2
            + (dphi_dC * dC) ** 2
            + (dphi_dD * dD) ** 2
        )
        out['phi'] = phi
        out['dphi'] = dphi

    return out


def format_summary_table(label, rows):
    """Format ABCD independence metrics into a readable text report."""
    lines = [
        f"=== ABCD Independence Summary: {label} ===",
        f"{'Category':<12} | {'R_SS = TSS/LSS':<18} | {'R_OS = TOS/LOS':<18} | {'K = R_OS/R_SS':<18} | {'Closure Diff (%)':<18} | {'Pearson phi':<18}",
        "-" * 116,
    ]
    for name, m in rows:
        r_ss_s = f"{m['r_ss']:.4f} +/- {m['dr_ss']:.4f}" if m['r_ss'] is not None else "N/A"
        r_os_s = f"{m['r_os']:.4f} +/- {m['dr_os']:.4f}" if m['r_os'] is not None else "N/A"
        k_s = f"{m['k_os_ss']:.4f} +/- {m['dk_os_ss']:.4f}" if m['k_os_ss'] is not None else "N/A"
        diff_s = f"{100.0 * m['diff']:+.2f} +/- {100.0 * m['ddiff']:.2f}" if m['diff'] is not None else "N/A"
        phi_s = f"{m['phi']:+.4f} +/- {m['dphi']:.4f}" if m['phi'] is not None else "N/A"
        lines.append(f"{name:<12} | {r_ss_s:<18} | {r_os_s:<18} | {k_s:<18} | {diff_s:<18} | {phi_s:<18}")
    return "\n".join(lines) + "\n"


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Plot TauFF ABCD independence checks (OS vs SS MisID probability and Pearson phi) for Run 3"
    )
    parser.add_argument("-V", "--version", required=True, type=str, help="Version/input prefix (e.g. process_1111_v1 or process_1111_v1_fake)")
    parser.add_argument("-C", "--channel", required=True, choices=["muon", "electron"], help="Analysis channel")
    parser.add_argument("-Y", "--years", nargs="+", default=["2022", "2022EE", "2023", "2023BPix", "2024"], help="Run 3 eras to combine")
    parser.add_argument("-O", "--outdir", default="", help="Output directory (default: ff_independence_<channel>)")
    parser.add_argument("--postfix", default="", help="Optional postfix appended to figure_<year><postfix>")
    parser.add_argument("--pt-bins", nargs="+", type=float, default=[40.0, 60.0, 80.0, 100.0, 140.0, 200.0, 400.0], help="Bin edges for tau pT MisID plots")
    parser.add_argument("--per-era", action="store_true", default=False, help="Also produce per-era plots in addition to combined Run 3")
    args = parser.parse_args(argv)
    if len(args.pt_bins) < 2 or sorted(args.pt_bins) != list(args.pt_bins) or len(set(args.pt_bins)) != len(args.pt_bins):
        fail("--pt-bins must contain at least 2 strictly increasing bin edges.")
    return args


def _extract_histograms(ROOT, plots_root_path, hist_name):
    """Extract detached (mc_total_TH1, data_TH1) clones from a plotIt plots.root canvas."""
    tf = ROOT.TFile.Open(plots_root_path, "READ")
    if not tf or tf.IsZombie():
        fail(f"Failed to open ROOT file: '{plots_root_path}'")
    try:
        canvas = tf.Get(hist_name)
        if not canvas:
            fail(f"Canvas '{hist_name}' not found in '{plots_root_path}'")
        pad_hi = canvas.GetPrimitive("pad_hi")
        if not pad_hi:
            fail(f"Primitive 'pad_hi' not found in canvas '{hist_name}' of '{plots_root_path}'")
        stack = pad_hi.GetPrimitive("mc_stack_0")
        if not stack or not stack.GetStack() or not stack.GetStack().Last():
            fail(f"MC stack 'mc_stack_0' missing/empty in '{hist_name}' of '{plots_root_path}'")
        h_mc = stack.GetStack().Last().Clone(f"{hist_name}_mc_tmp")
        h_mc.SetDirectory(0)
        if not h_mc.GetSumw2N():
            h_mc.Sumw2()

        h_data = None
        for prim in pad_hi.GetListOfPrimitives():
            if prim.InheritsFrom("TH1"):
                h_data = prim.Clone(f"{hist_name}_data_tmp")
                h_data.SetDirectory(0)
                if not h_data.GetSumw2N():
                    h_data.Sumw2()
                break
        if h_data is None:
            fail(f"Data TH1 primitive not found in '{hist_name}' of '{plots_root_path}'")
        return h_mc, h_data
    finally:
        tf.Close()


def _load_region_hists(ROOT, version, channel, years, postfix):
    """Load and sum MC-total, genuine-tau MC, and Data histograms per region across given eras."""
    hist_specs = [
        ("pt", "h_tau1_pt_S3", "h_tau1_gen_pt_S3"),
        ("dm", "h_tau1_decayMode_S3", "h_tau1_decayMode_gen_S3"),
    ]
    data = {reg: {} for reg in REGIONS}

    for reg in REGIONS:
        for yr in years:
            path = resolve_region_plot_file(version, reg, channel, yr, postfix)
            for tag, all_name, gen_name in hist_specs:
                h_mc, h_data = _extract_histograms(ROOT, path, all_name)
                h_gen, _ = _extract_histograms(ROOT, path, gen_name)

                h_nsub_mc = h_mc.Clone(f"nsub_mc_{reg}_{tag}_{yr}")
                h_nsub_mc.SetDirectory(0)
                h_nsub_mc.Add(h_gen, -1.0)

                h_nsub_data = h_data.Clone(f"nsub_data_{reg}_{tag}_{yr}")
                h_nsub_data.SetDirectory(0)
                h_nsub_data.Add(h_gen, -1.0)

                if tag not in data[reg]:
                    data[reg][tag] = {
                        "nsub_mc": h_nsub_mc,
                        "nsub_data": h_nsub_data,
                    }
                else:
                    data[reg][tag]["nsub_mc"].Add(h_nsub_mc)
                    data[reg][tag]["nsub_data"].Add(h_nsub_data)

    return data


def _compact_dm_hist(ROOT, src_hist, name):
    """Project 15-bin decayMode histogram onto 4 labeled bins (DM 0, 1, 10, 11)."""
    h_out = ROOT.TH1D(name, "", len(DM_BINS), 0.0, float(len(DM_BINS)))
    h_out.SetDirectory(0)
    h_out.Sumw2()
    for out_idx, (label, raw_bin) in enumerate(DM_BINS, start=1):
        h_out.GetXaxis().SetBinLabel(out_idx, label)
        h_out.SetBinContent(out_idx, src_hist.GetBinContent(raw_bin))
        h_out.SetBinError(out_idx, src_hist.GetBinError(raw_bin))
    return h_out


def _hist_integral_and_err(hist, first_bin=1, last_bin=-1):
    err = ctypes.c_double(0.0)
    if last_bin < 0:
        last_bin = hist.GetNbinsX()
    val = hist.IntegralAndError(first_bin, last_bin, err)
    return float(val), float(err.value)


def _draw_two_pad_ratio(
    ROOT,
    h_num,
    h_den,
    num_label,
    den_label,
    x_title,
    y_title,
    ratio_title,
    channel,
    era_label,
    out_base,
    summary_text="",
):
    """Draw a publication-style upper comparison + lower ratio canvas and save PDF/PNG."""
    c = ROOT.TCanvas(f"c_{os.path.basename(out_base)}", "", 650, 650)
    pad_hi = ROOT.TPad("pad_hi", "pad_hi", 0.0, 0.30, 1.0, 1.0)
    pad_lo = ROOT.TPad("pad_lo", "pad_lo", 0.0, 0.0, 1.0, 0.30)

    pad_hi.SetBottomMargin(0.025)
    pad_hi.SetLeftMargin(0.14)
    pad_hi.SetRightMargin(0.05)
    pad_hi.SetTopMargin(0.08)

    pad_lo.SetTopMargin(0.035)
    pad_lo.SetBottomMargin(0.32)
    pad_lo.SetLeftMargin(0.14)
    pad_lo.SetRightMargin(0.05)
    pad_lo.SetGridy(1)

    pad_hi.Draw()
    pad_lo.Draw()

    # Upper pad
    pad_hi.cd()
    h1 = h_num.Clone(f"{h_num.GetName()}_draw1")
    h2 = h_den.Clone(f"{h_den.GetName()}_draw2")
    h1.SetDirectory(0)
    h2.SetDirectory(0)

    h1.SetLineColor(ROOT.kAzure + 1)
    h1.SetMarkerColor(ROOT.kAzure + 1)
    h1.SetMarkerStyle(20)
    h1.SetMarkerSize(1.0)
    h1.SetLineWidth(2)

    h2.SetLineColor(ROOT.kRed + 1)
    h2.SetMarkerColor(ROOT.kRed + 1)
    h2.SetMarkerStyle(24)
    h2.SetMarkerSize(1.0)
    h2.SetLineWidth(2)

    ymax = max(h1.GetMaximum(), h2.GetMaximum(), 0.05)
    h1.GetYaxis().SetTitle(y_title)
    h1.GetYaxis().SetTitleSize(0.052)
    h1.GetYaxis().SetLabelSize(0.045)
    h1.GetYaxis().SetTitleOffset(1.15)
    h1.GetXaxis().SetLabelSize(0.0)
    h1.GetYaxis().SetRangeUser(0.0, max(0.5, 1.45 * ymax))
    h1.Draw("E1")
    h2.Draw("E1 SAME")

    leg = ROOT.TLegend(0.55, 0.72, 0.92, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.04)
    leg.AddEntry(h1, num_label, "lep")
    leg.AddEntry(h2, den_label, "lep")
    leg.Draw()

    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextFont(62)
    lat.SetTextSize(0.045)
    lat.DrawLatex(0.14, 0.935, "CMS")
    lat.SetTextFont(52)
    lat.SetTextSize(0.038)
    lat.DrawLatex(0.225, 0.935, "Work in progress")
    lat.SetTextFont(42)
    lat.SetTextAlign(31)
    ch_tex = "#mu#tau_{h}" if channel == "muon" else "e#tau_{h}"
    lat.DrawLatex(0.95, 0.935, f"{ch_tex} | {era_label}")

    if summary_text:
        lat.SetTextAlign(11)
        lat.SetTextSize(0.036)
        lat.DrawLatex(0.18, 0.83, summary_text)

    # Lower ratio pad
    pad_lo.cd()
    h_ratio = h1.Clone(f"{h_num.GetName()}_ratio")
    h_ratio.SetDirectory(0)
    h_ratio.Divide(h1, h2, 1.0, 1.0)
    h_ratio.SetLineColor(ROOT.kBlack)
    h_ratio.SetMarkerColor(ROOT.kBlack)
    h_ratio.SetMarkerStyle(20)
    h_ratio.SetMarkerSize(0.9)
    h_ratio.GetYaxis().SetTitle(ratio_title)
    h_ratio.GetYaxis().SetTitleSize(0.105)
    h_ratio.GetYaxis().SetLabelSize(0.095)
    h_ratio.GetYaxis().SetTitleOffset(0.55)
    h_ratio.GetYaxis().SetNdivisions(505)
    h_ratio.GetYaxis().SetRangeUser(0.0, 2.2)
    h_ratio.GetXaxis().SetTitle(x_title)
    h_ratio.GetXaxis().SetTitleSize(0.115)
    h_ratio.GetXaxis().SetLabelSize(0.105)
    h_ratio.GetXaxis().SetTitleOffset(1.05)
    h_ratio.Draw("E1")

    xmin = h_ratio.GetXaxis().GetXmin()
    xmax = h_ratio.GetXaxis().GetXmax()
    unity = ROOT.TLine(xmin, 1.0, xmax, 1.0)
    unity.SetLineStyle(2)
    unity.SetLineWidth(1)
    unity.SetLineColor(ROOT.kGray + 2)
    unity.Draw("SAME")

    c.SaveAs(out_base + ".pdf")
    c.SaveAs(out_base + ".png")
    c.Close()


def _run_for_era_set(ROOT, args, years, tag, era_label, outdir):
    import array as py_array

    data = _load_region_hists(ROOT, args.version, args.channel, years, args.postfix)
    pt_edges = py_array.array("d", args.pt_bins)
    n_pt = len(pt_edges) - 1

    # Rebin pT and project DM
    rebinned = {}
    for reg in REGIONS:
        rebinned[reg] = {
            "pt_mc": data[reg]["pt"]["nsub_mc"].Rebin(n_pt, f"pt_mc_{reg}_{tag}", pt_edges),
            "pt_data": data[reg]["pt"]["nsub_data"].Rebin(n_pt, f"pt_data_{reg}_{tag}", pt_edges),
            "dm_mc": _compact_dm_hist(ROOT, data[reg]["dm"]["nsub_mc"], f"dm_mc_{reg}_{tag}"),
            "dm_data": _compact_dm_hist(ROOT, data[reg]["dm"]["nsub_data"], f"dm_data_{reg}_{tag}"),
        }
        for h in rebinned[reg].values():
            h.SetDirectory(0)

    # Compute ABCD summary metrics (Inclusive + per DM) on MC fake yields
    A_inc, dA_inc = _hist_integral_and_err(rebinned["lss"]["dm_mc"])
    B_inc, dB_inc = _hist_integral_and_err(rebinned["tss"]["dm_mc"])
    C_inc, dC_inc = _hist_integral_and_err(rebinned["los"]["dm_mc"])
    D_inc, dD_inc = _hist_integral_and_err(rebinned["tos"]["dm_mc"])
    rows = [("Inclusive", compute_abcd_metrics(A_inc, dA_inc, B_inc, dB_inc, C_inc, dC_inc, D_inc, dD_inc))]

    for idx, dm_name in enumerate(DM_SHORT_NAMES[1:], start=1):
        A, dA = rebinned["lss"]["dm_mc"].GetBinContent(idx), rebinned["lss"]["dm_mc"].GetBinError(idx)
        B, dB = rebinned["tss"]["dm_mc"].GetBinContent(idx), rebinned["tss"]["dm_mc"].GetBinError(idx)
        C, dC = rebinned["los"]["dm_mc"].GetBinContent(idx), rebinned["los"]["dm_mc"].GetBinError(idx)
        D, dD = rebinned["tos"]["dm_mc"].GetBinContent(idx), rebinned["tos"]["dm_mc"].GetBinError(idx)
        m_dm = compute_abcd_metrics(A, dA, B, dB, C, dC, D, dD)
        if m_dm["k_os_ss"] is None:
            print(f"[WARNING] Undefined ABCD ratio for {tag} {dm_name} (A={A:.2f}, B={B:.2f}, C={C:.2f}, D={D:.2f})")
        rows.append((dm_name, m_dm))

    report = format_summary_table(f"{args.channel} / {era_label}", rows)
    print(report)
    summary_path = os.path.join(outdir, f"independence_summary_{tag}.txt")
    with open(summary_path, "w") as f:
        f.write(report)

    inc_m = rows[0][1]
    ann = ""
    if inc_m["phi"] is not None and inc_m["k_os_ss"] is not None:
        ann = f"#phi = {inc_m['phi']:+.4f} #pm {inc_m['dphi']:.4f},  K = {inc_m['k_os_ss']:.3f} #pm {inc_m['dk_os_ss']:.3f}"

    # 1. MisID prob vs pT: OS MC vs SS MC
    r_os_pt = rebinned["tos"]["pt_mc"].Clone(f"r_os_pt_{tag}")
    r_ss_pt = rebinned["tss"]["pt_mc"].Clone(f"r_ss_pt_{tag}")
    r_os_pt.SetDirectory(0)
    r_ss_pt.SetDirectory(0)
    r_os_pt.Divide(rebinned["tos"]["pt_mc"], rebinned["los"]["pt_mc"], 1.0, 1.0)
    r_ss_pt.Divide(rebinned["tss"]["pt_mc"], rebinned["lss"]["pt_mc"], 1.0, 1.0)
    _draw_two_pad_ratio(
        ROOT, r_os_pt, r_ss_pt,
        "OS MC (TOS / LOS)", "SS MC (TSS / LSS)",
        "#tau_{h} p_{T} (GeV)", "MisID Probability (Tight / Loose)", "R_{OS} / R_{SS}",
        args.channel, era_label, os.path.join(outdir, f"misid_prob_os_vs_ss_pt_{tag}"), ann,
    )

    # 2. MisID prob vs DM: OS MC vs SS MC
    r_os_dm = rebinned["tos"]["dm_mc"].Clone(f"r_os_dm_{tag}")
    r_ss_dm = rebinned["tss"]["dm_mc"].Clone(f"r_ss_dm_{tag}")
    r_os_dm.SetDirectory(0)
    r_ss_dm.SetDirectory(0)
    r_os_dm.Divide(rebinned["tos"]["dm_mc"], rebinned["los"]["dm_mc"], 1.0, 1.0)
    r_ss_dm.Divide(rebinned["tss"]["dm_mc"], rebinned["lss"]["dm_mc"], 1.0, 1.0)
    _draw_two_pad_ratio(
        ROOT, r_os_dm, r_ss_dm,
        "OS MC (TOS / LOS)", "SS MC (TSS / LSS)",
        "#tau_{h} Decay Mode", "MisID Probability (Tight / Loose)", "R_{OS} / R_{SS}",
        args.channel, era_label, os.path.join(outdir, f"misid_prob_os_vs_ss_dm_{tag}"), ann,
    )

    # 3. SS closure vs pT: SS Data vs SS MC
    r_ss_data_pt = rebinned["tss"]["pt_data"].Clone(f"r_ss_data_pt_{tag}")
    r_ss_data_pt.SetDirectory(0)
    r_ss_data_pt.Divide(rebinned["tss"]["pt_data"], rebinned["lss"]["pt_data"], 1.0, 1.0)
    _draw_two_pad_ratio(
        ROOT, r_ss_data_pt, r_ss_pt,
        "SS Data (TSS / LSS)", "SS MC (TSS / LSS)",
        "#tau_{h} p_{T} (GeV)", "MisID Probability (Tight / Loose)", "Data / MC",
        args.channel, era_label, os.path.join(outdir, f"misid_prob_ss_data_vs_mc_pt_{tag}"),
    )

    # 4. SS closure vs DM: SS Data vs SS MC
    r_ss_data_dm = rebinned["tss"]["dm_data"].Clone(f"r_ss_data_dm_{tag}")
    r_ss_data_dm.SetDirectory(0)
    r_ss_data_dm.Divide(rebinned["tss"]["dm_data"], rebinned["lss"]["dm_data"], 1.0, 1.0)
    _draw_two_pad_ratio(
        ROOT, r_ss_data_dm, r_ss_dm,
        "SS Data (TSS / LSS)", "SS MC (TSS / LSS)",
        "#tau_{h} Decay Mode", "MisID Probability (Tight / Loose)", "Data / MC",
        args.channel, era_label, os.path.join(outdir, f"misid_prob_ss_data_vs_mc_dm_{tag}"),
    )

    root_out = ROOT.TFile.Open(os.path.join(outdir, f"tauff_independence_{tag}.root"), "RECREATE")
    for obj in (r_os_pt, r_ss_pt, r_os_dm, r_ss_dm, r_ss_data_pt, r_ss_data_dm):
        obj.Write()
    root_out.Close()


def main(argv=None):
    args = parse_args(argv)
    import ROOT

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)

    outdir = args.outdir or f"ff_independence_{args.channel}"
    os.makedirs(outdir, exist_ok=True)

    run3_tag = "Run3" if len(args.years) > 1 else args.years[0]
    run3_label = f"Run 3 ({'+'.join(args.years)})" if len(args.years) > 1 else args.years[0]
    _run_for_era_set(ROOT, args, args.years, run3_tag, run3_label, outdir)

    if args.per_era and len(args.years) > 1:
        for yr in args.years:
            _run_for_era_set(ROOT, args, [yr], yr, yr, outdir)

    print(f"[INFO] Saved TauFF independence plots and summary to '{outdir}'")


if __name__ == "__main__":
    main()

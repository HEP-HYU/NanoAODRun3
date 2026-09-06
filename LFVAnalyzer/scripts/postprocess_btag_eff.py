#!/usr/bin/env python3
"""
postprocess_btag_eff.py

Calculates 2D b-tag efficiency maps from numerator and denominator histograms
produced by BTagEfficiencyAnalyzer (or merged across dataset chunks).
Performs binomial error computation (Divide with "B" option) and empty bin smoothing.

Usage:
  python postprocess_btag_eff.py -I <raw_analyzer_output.root> -O <final_btag_eff.root>
"""

import sys
import os
import argparse
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

DEFAULT_FALLBACK = {
    "b": 0.65,
    "c": 0.15,
    "light": 0.01,
}

def find_hist(tfile, candidates):
    """Search for a histogram name among a list of candidates (e.g. with/without _S4 suffix)."""
    for name in candidates:
        h = tfile.Get(name)
        if h and not h.IsZombie():
            return h
    return None

def process_file(in_path, out_path, fallback_eff=None):
    if fallback_eff is None:
        fallback_eff = DEFAULT_FALLBACK

    if not os.path.isfile(in_path):
        print(f"[ERROR] Input file does not exist: {in_path}", file=sys.stderr)
        return False

    fin = ROOT.TFile.Open(in_path, "READ")
    if not fin or fin.IsZombie():
        print(f"[ERROR] Failed to open input ROOT file: {in_path}", file=sys.stderr)
        return False

    out_dir = os.path.dirname(os.path.abspath(out_path))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    fout = ROOT.TFile.Open(out_path, "RECREATE")
    if not fout or fout.IsZombie():
        print(f"[ERROR] Failed to create output ROOT file: {out_path}", file=sys.stderr)
        fin.Close()
        return False

    flavors = ["b", "c", "light"]
    for fl in flavors:
        denom_names = [f"h2_denom_{fl}_S4", f"h2_denom_{fl}"]
        num_names   = [f"h2_num_{fl}_S4", f"h2_num_{fl}"]

        h_denom = find_hist(fin, denom_names)
        h_num   = find_hist(fin, num_names)

        if not h_denom or not h_num:
            print(f"[WARNING] Histograms for flavor '{fl}' not found in {in_path}! Skipping...", file=sys.stderr)
            continue

        h_denom_cl = h_denom.Clone(f"h2_denom_{fl}")
        h_num_cl   = h_num.Clone(f"h2_num_{fl}")

        h_denom_cl.SetDirectory(fout)
        h_num_cl.SetDirectory(fout)
        h_denom_cl.Write()
        h_num_cl.Write()

        # Compute efficiency: num / denom with binomial errors
        h2_eff = h_num_cl.Clone(f"h2_eff_{fl}")
        h2_eff.SetTitle(f"B-tagging Efficiency ({fl});Jet p_{{T}} [GeV];Jet |#eta|")
        h2_eff.Divide(h_num_cl, h_denom_cl, 1.0, 1.0, "B")

        nbins_x = h2_eff.GetNbinsX()
        nbins_y = h2_eff.GetNbinsY()
        n_smoothed = 0

        for bx in range(1, nbins_x + 1):
            for by in range(1, nbins_y + 1):
                eff = h2_eff.GetBinContent(bx, by)
                denom_val = h_denom_cl.GetBinContent(bx, by)
                if denom_val <= 0 or eff <= 0.0 or eff > 1.0:
                    h2_eff.SetBinContent(bx, by, fallback_eff[fl])
                    h2_eff.SetBinError(bx, by, 0.0)
                    n_smoothed += 1

        h2_eff.SetDirectory(fout)
        h2_eff.Write()

        mean_eff = h2_eff.Integral() / (nbins_x * nbins_y) if (nbins_x * nbins_y) > 0 else 0.0
        print(f"  [FLAVOR {fl:>5}] Mean Eff: {mean_eff:.4f} (Smoothed bins: {n_smoothed}/{nbins_x * nbins_y})")

    # Also copy over control histograms if present
    for ctrl in ["h_nevents", "h_ncleanjetspass"]:
        h_ctrl = find_hist(fin, [f"{ctrl}_S4", ctrl])
        if h_ctrl:
            h_ctrl_cl = h_ctrl.Clone(ctrl)
            h_ctrl_cl.SetDirectory(fout)
            h_ctrl_cl.Write()

    fout.Close()
    fin.Close()
    print(f"[DONE] Final b-tag efficiency maps written to: {out_path}")
    return True

def main():
    parser = argparse.ArgumentParser(description="Post-process raw BTagEfficiencyAnalyzer output into finalized efficiency maps")
    parser.add_argument("-I", "--infile", required=True, help="Input raw ROOT file")
    parser.add_argument("-O", "--outfile", required=True, help="Output finalized ROOT file")
    args = parser.parse_args()

    success = process_file(args.infile, args.outfile)
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()

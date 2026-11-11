#!/usr/bin/env python3
"""
postprocess_btag_eff.py

Calculates 2D b-tag efficiency maps from numerator and denominator histograms
produced by BTagEfficiencyAnalyzer (or merged across dataset chunks).

Option 2-3 smoothing policy (CMS BTV Method 1a):
  1. Each (pT, |eta|) bin is accepted directly from the process sample iff:
       denom > 0,  0.0 < eff < 1.0,  and  N_eff = denom^2 / sigma_denom^2 >= min_neff (default: 20.0).
  2. For base maps without an external SM reference (e.g. btag_eff_ttbar.root or
     btag_eff_signal.root), any tail/low-statistics bin failing (1) is smoothed
     using the nearest valid pT bin in the same (flavor, |eta|) slice.
  3. For SM process maps (singletop, ttx, dyjets, wjets, diboson, qcd, inclusive),
     any bin failing (1) is replaced by the corresponding (pT, |eta|) bin from the
     smoothed SM reference map (default: <out_dir>/btag_eff_ttbar.root).
  4. Never falls back to hardcoded constants; prints explicit [ERROR] and exits
     non-zero if neither a valid bin nor a valid reference map is available.

Usage:
  python postprocess_btag_eff.py -I <raw_analyzer_output.root> -O <final_btag_eff.root> [--ref-file <ref.root>] [--min-neff 20.0]
"""

import sys
import os
import math
import argparse

FLAVORS = ["b", "c", "light"]
SELF_SMOOTH_BASENAMES = {
    "btag_eff_ttbar.root",
    "btag_eff_signal.root",
    "btag_eff_signal_tc.root",
    "btag_eff_signal_tu.root",
}


def resolve_ref_path(out_path, ref_path):
    """Auto-discover <out_dir>/btag_eff_ttbar.root for SM non-ttbar maps unless overridden."""
    if ref_path:
        if ref_path.lower() == "none":
            return None
        return ref_path

    abs_out = os.path.abspath(out_path)
    out_base = os.path.basename(abs_out)
    if out_base in SELF_SMOOTH_BASENAMES:
        return None

    cand = os.path.join(os.path.dirname(abs_out), "btag_eff_ttbar.root")
    if os.path.abspath(cand) == abs_out:
        return None
    return cand


def compute_smoothed_grid(denom_vals, denom_vars, num_vals, num_vars, min_neff, ref_vals=None, ref_vars=None, fl="", in_path="", ref_path=None):
    """
    Computes (eff_vals, eff_vars, n_smoothed) over 2D shape (nx, ny).
    Raises RuntimeError with [ERROR] prefix on unrecoverable invalid inputs.
    """
    nx, ny = denom_vals.shape
    eff_vals = [[0.0] * ny for _ in range(nx)]
    eff_vars = [[0.0] * ny for _ in range(nx)]
    valid_mask = [[False] * ny for _ in range(nx)]

    for bx in range(nx):
        for by in range(ny):
            d_val = float(denom_vals[bx][by])
            d_var = float(denom_vars[bx][by])
            n_val = float(num_vals[bx][by])
            n_var = float(num_vars[bx][by])

            if d_val > 0.0 and d_var > 0.0:
                raw_eff = n_val / d_val
                neff = (d_val * d_val) / d_var
                if math.isfinite(raw_eff) and 0.0 < raw_eff < 1.0 and neff >= min_neff:
                    # Exact ROOT TH1::Divide(..., "B") binomial variance for weighted histograms
                    var_b = abs(((1.0 - 2.0 * raw_eff) * n_var + (raw_eff * raw_eff) * d_var) / (d_val * d_val))
                    eff_vals[bx][by] = raw_eff
                    eff_vars[bx][by] = var_b
                    valid_mask[bx][by] = True

    n_smoothed = 0
    for by in range(ny):
        valid_bxs = [bx for bx in range(nx) if valid_mask[bx][by]]
        for bx in range(nx):
            if valid_mask[bx][by]:
                continue
            n_smoothed += 1
            if ref_vals is not None:
                r_eff = float(ref_vals[bx][by])
                r_var = float(ref_vars[bx][by]) if ref_vars is not None else 0.0
                if not math.isfinite(r_eff) or r_eff <= 0.0 or r_eff >= 1.0:
                    raise RuntimeError(
                        f"[ERROR] postprocess_btag_eff: Reference map '{ref_path}' has invalid efficiency "
                        f"({r_eff}) at flavor='{fl}', bin=({bx + 1},{by + 1}) while smoothing '{in_path}'."
                    )
                eff_vals[bx][by] = r_eff
                eff_vars[bx][by] = max(0.0, r_var)
            else:
                if not valid_bxs:
                    raise RuntimeError(
                        f"[ERROR] postprocess_btag_eff: No valid bin with N_eff >= {min_neff} and 0 < eff < 1 "
                        f"for flavor='{fl}', eta_bin={by + 1} in '{in_path}', and no reference map provided."
                    )
                # Nearest valid pT bin (prefer lower-pT neighbor on tie)
                nearest_bx = min(valid_bxs, key=lambda k: (abs(k - bx), k > bx))
                eff_vals[bx][by] = eff_vals[nearest_bx][by]
                eff_vars[bx][by] = eff_vars[nearest_bx][by]

    return eff_vals, eff_vars, n_smoothed


def _find_hist_pyroot(tfile, candidates):
    for name in candidates:
        h = tfile.Get(name)
        if h and not h.IsZombie():
            return h
    return None


def _process_with_pyroot(ROOT, in_path, out_path, ref_path, min_neff):
    fin = ROOT.TFile.Open(in_path, "READ")
    if not fin or fin.IsZombie():
        print(f"[ERROR] postprocess_btag_eff: Failed to open input ROOT file: {in_path}", file=sys.stderr)
        return False

    fref = None
    if ref_path is not None:
        if not os.path.isfile(ref_path):
            print(
                f"[ERROR] postprocess_btag_eff: Required SM reference efficiency map does not exist: {ref_path}",
                file=sys.stderr,
            )
            fin.Close()
            return False
        fref = ROOT.TFile.Open(ref_path, "READ")
        if not fref or fref.IsZombie():
            print(f"[ERROR] postprocess_btag_eff: Failed to open reference ROOT file: {ref_path}", file=sys.stderr)
            fin.Close()
            return False

    tmp_out = out_path + ".tmp.root"
    fout = ROOT.TFile.Open(tmp_out, "RECREATE")
    if not fout or fout.IsZombie():
        print(f"[ERROR] postprocess_btag_eff: Failed to create output ROOT file: {tmp_out}", file=sys.stderr)
        if fref:
            fref.Close()
        fin.Close()
        return False

    try:
        for fl in FLAVORS:
            denom_names = [f"h2_denom_{fl}_S4", f"h2_denom_{fl}", f"h2_denom_{fl}_S0000"]
            num_names   = [f"h2_num_{fl}_S4", f"h2_num_{fl}", f"h2_num_{fl}_S0000"]

            h_denom = _find_hist_pyroot(fin, denom_names)
            h_num   = _find_hist_pyroot(fin, num_names)
            if not h_denom or not h_num:
                raise RuntimeError(
                    f"[ERROR] postprocess_btag_eff: Histograms for flavor '{fl}' not found in {in_path}!"
                )

            nx = h_denom.GetNbinsX()
            ny = h_denom.GetNbinsY()

            d_vals = [[h_denom.GetBinContent(bx, by) for by in range(1, ny + 1)] for bx in range(1, nx + 1)]
            d_vars = [[h_denom.GetBinError(bx, by) ** 2 for by in range(1, ny + 1)] for bx in range(1, nx + 1)]
            n_vals = [[h_num.GetBinContent(bx, by) for by in range(1, ny + 1)] for bx in range(1, nx + 1)]
            n_vars = [[h_num.GetBinError(bx, by) ** 2 for by in range(1, ny + 1)] for bx in range(1, nx + 1)]

            r_vals = None
            r_vars = None
            if fref is not None:
                h_ref = _find_hist_pyroot(fref, [f"h2_eff_{fl}"])
                if not h_ref or h_ref.GetNbinsX() != nx or h_ref.GetNbinsY() != ny:
                    raise RuntimeError(
                        f"[ERROR] postprocess_btag_eff: Missing or incompatible 'h2_eff_{fl}' in reference '{ref_path}'."
                    )
                r_vals = [[h_ref.GetBinContent(bx, by) for by in range(1, ny + 1)] for bx in range(1, nx + 1)]
                r_vars = [[h_ref.GetBinError(bx, by) ** 2 for by in range(1, ny + 1)] for bx in range(1, nx + 1)]

            eff_vals, eff_vars, n_smoothed = compute_smoothed_grid(
                d_vals, d_vars, n_vals, n_vars, min_neff,
                ref_vals=r_vals, ref_vars=r_vars, fl=fl, in_path=in_path, ref_path=ref_path,
            )

            h_denom_cl = h_denom.Clone(f"h2_denom_{fl}")
            h_num_cl   = h_num.Clone(f"h2_num_{fl}")
            h_denom_cl.SetDirectory(fout)
            h_num_cl.SetDirectory(fout)
            h_denom_cl.Write()
            h_num_cl.Write()

            h2_eff = h_num_cl.Clone(f"h2_eff_{fl}")
            h2_eff.Reset("ICES")
            h2_eff.SetTitle(f"B-tagging Efficiency ({fl});Jet p_{{T}} [GeV];Jet |#eta|")
            for bx in range(1, nx + 1):
                for by in range(1, ny + 1):
                    h2_eff.SetBinContent(bx, by, eff_vals[bx - 1][by - 1])
                    h2_eff.SetBinError(bx, by, math.sqrt(eff_vars[bx - 1][by - 1]))

            h2_eff.SetDirectory(fout)
            h2_eff.Write()

            mean_eff = h2_eff.Integral() / (nx * ny) if (nx * ny) > 0 else 0.0
            src_note = f"ref={os.path.basename(ref_path)}" if ref_path else "self-neighbor"
            print(f"  [FLAVOR {fl:>5}] Mean Eff: {mean_eff:.4f} (Smoothed bins: {n_smoothed}/{nx * ny}, {src_note})")

        for ctrl in ["h_nevents", "h_ncleanjetspass"]:
            h_ctrl = _find_hist_pyroot(fin, [f"{ctrl}_S4", ctrl, f"{ctrl}_S0000"])
            if h_ctrl:
                h_ctrl_cl = h_ctrl.Clone(ctrl)
                h_ctrl_cl.SetDirectory(fout)
                h_ctrl_cl.Write()

    except Exception as exc:
        print(str(exc), file=sys.stderr)
        fout.Close()
        if fref:
            fref.Close()
        fin.Close()
        if os.path.isfile(tmp_out):
            os.remove(tmp_out)
        return False

    fout.Close()
    if fref:
        fref.Close()
    fin.Close()
    os.replace(tmp_out, out_path)
    print(f"[DONE] Final b-tag efficiency maps written to: {out_path}")
    return True


def _find_key_uproot(uf, candidates):
    present = {k.split(";")[0]: k for k in uf.keys()}
    for name in candidates:
        if name in present:
            return uf[present[name]]
    return None


def _process_with_uproot(in_path, out_path, ref_path, min_neff):
    import numpy as np
    import uproot
    import uproot.writing.identify as uid

    if ref_path is not None and not os.path.isfile(ref_path):
        print(
            f"[ERROR] postprocess_btag_eff: Required SM reference efficiency map does not exist: {ref_path}",
            file=sys.stderr,
        )
        return False

    tmp_out = out_path + ".tmp.root"
    try:
        with uproot.open(in_path) as fin:
            ref_hists = {}
            if ref_path is not None:
                with uproot.open(ref_path) as fref:
                    for fl in FLAVORS:
                        h_ref = _find_key_uproot(fref, [f"h2_eff_{fl}"])
                        if h_ref is None:
                            raise RuntimeError(
                                f"[ERROR] postprocess_btag_eff: Missing 'h2_eff_{fl}' in reference '{ref_path}'."
                            )
                        ref_hists[fl] = (h_ref.values().copy(), h_ref.variances().copy())

            out_objs = {}
            for fl in FLAVORS:
                h_denom = _find_key_uproot(fin, [f"h2_denom_{fl}_S4", f"h2_denom_{fl}", f"h2_denom_{fl}_S0000"])
                h_num   = _find_key_uproot(fin, [f"h2_num_{fl}_S4", f"h2_num_{fl}", f"h2_num_{fl}_S0000"])
                if h_denom is None or h_num is None:
                    raise RuntimeError(
                        f"[ERROR] postprocess_btag_eff: Histograms for flavor '{fl}' not found in {in_path}!"
                    )

                d_vals = h_denom.values()
                d_vars = h_denom.variances()
                n_vals = h_num.values()
                n_vars = h_num.variances()
                nx, ny = d_vals.shape

                r_vals, r_vars = (None, None)
                if ref_path is not None:
                    r_vals, r_vars = ref_hists[fl]
                    if r_vals.shape != (nx, ny):
                        raise RuntimeError(
                            f"[ERROR] postprocess_btag_eff: Shape mismatch for 'h2_eff_{fl}' in reference '{ref_path}'."
                        )

                eff_vals, eff_vars, n_smoothed = compute_smoothed_grid(
                    d_vals, d_vars, n_vals, n_vars, min_neff,
                    ref_vals=r_vals, ref_vars=r_vars, fl=fl, in_path=in_path, ref_path=ref_path,
                )

                out_objs[f"h2_denom_{fl}"] = h_denom
                out_objs[f"h2_num_{fl}"] = h_num

                eff_flow = np.zeros((nx + 2, ny + 2), dtype=np.float64)
                var_flow = np.zeros((nx + 2, ny + 2), dtype=np.float64)
                eff_arr = np.asarray(eff_vals, dtype=np.float64)
                var_arr = np.asarray(eff_vars, dtype=np.float64)
                eff_flow[1:-1, 1:-1] = eff_arr
                var_flow[1:-1, 1:-1] = var_arr

                # ROOT TH2 internal 1D cell array is C-contiguous over (ny + 2, nx + 2)
                data_1d  = np.ascontiguousarray(eff_flow.T).ravel()
                sumw2_1d = np.ascontiguousarray(var_flow.T).ravel()
                tsumw = float(eff_arr.sum())
                tsumw2 = float(var_arr.sum())

                out_objs[f"h2_eff_{fl}"] = uid.to_TH2x(
                    f"h2_eff_{fl}",
                    f"B-tagging Efficiency ({fl})",
                    data_1d,
                    float(nx * ny),
                    tsumw,
                    tsumw2,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    sumw2_1d,
                    h_denom.member("fXaxis"),
                    h_denom.member("fYaxis"),
                    h_denom.member("fZaxis"),
                )

                mean_eff = tsumw / (nx * ny) if (nx * ny) > 0 else 0.0
                src_note = f"ref={os.path.basename(ref_path)}" if ref_path else "self-neighbor"
                print(f"  [FLAVOR {fl:>5}] Mean Eff: {mean_eff:.4f} (Smoothed bins: {n_smoothed}/{nx * ny}, {src_note})")

            for ctrl in ["h_nevents", "h_ncleanjetspass"]:
                h_ctrl = _find_key_uproot(fin, [f"{ctrl}_S4", ctrl, f"{ctrl}_S0000"])
                if h_ctrl is not None:
                    out_objs[ctrl] = h_ctrl

            with uproot.recreate(tmp_out) as fout:
                for k, obj in out_objs.items():
                    fout[k] = obj

    except Exception as exc:
        print(str(exc), file=sys.stderr)
        if os.path.isfile(tmp_out):
            os.remove(tmp_out)
        return False

    os.replace(tmp_out, out_path)
    print(f"[DONE] Final b-tag efficiency maps written to: {out_path}")
    return True


def process_file(in_path, out_path, ref_path=None, min_neff=20.0):
    if not os.path.isfile(in_path):
        print(f"[ERROR] postprocess_btag_eff: Input file does not exist: {in_path}", file=sys.stderr)
        return False

    out_dir = os.path.dirname(os.path.abspath(out_path))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    resolved_ref = resolve_ref_path(out_path, ref_path)

    try:
        import ROOT
        ROOT.gROOT.SetBatch(True)
        ROOT.PyConfig.IgnoreCommandLineOptions = True
        return _process_with_pyroot(ROOT, in_path, out_path, resolved_ref, min_neff)
    except ImportError:
        return _process_with_uproot(in_path, out_path, resolved_ref, min_neff)


def main():
    parser = argparse.ArgumentParser(description="Post-process raw BTagEfficiencyAnalyzer output into finalized efficiency maps")
    parser.add_argument("-I", "--infile", required=True, help="Input raw ROOT file")
    parser.add_argument("-O", "--outfile", required=True, help="Output finalized ROOT file")
    parser.add_argument("--ref-file", default="", help="SM reference efficiency ROOT file for low-stat bin replacement (default: auto <out_dir>/btag_eff_ttbar.root, or 'none' for self-smoothing)")
    parser.add_argument("--min-neff", type=float, default=20.0, help="Minimum effective unweighted denominator entries N_eff = (sum w)^2 / (sum w^2) per bin [default: 20.0]")
    args = parser.parse_args()

    success = process_file(args.infile, args.outfile, ref_path=args.ref_file, min_neff=args.min_neff)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()

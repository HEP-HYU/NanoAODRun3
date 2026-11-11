#!/usr/bin/env python3
"""
compute_btag_efficiency_rdataframe.py

Executes the C++ BTagEfficiencyAnalyzer class via ROOT RDataFrame to extract
per-flavor 2D b-tagging efficiency maps at native C++ speed.

Implements CMS BTV Method 1a sample combination:
  1. Runs BTagEfficiencyAnalyzer per MC sample directory with a sample-specific
     output filename (ensuring TopPtWeight is properly activated for TTto2L2Nu,
     TTtoLNu2Q, TTto4Q in CommonTools/src/impl_Reweight.cpp).
  2. Extracts the pre-skim generator weight sum (sum_gen_weight) from bin 2 of
     'hcounter' across all ROOT files of each sample.
  3. Scales each sample's raw denominator and numerator histograms by:
       - SM processes (ttbar, singletop, wjets, dyjets, qcd, diboson, ttx, inclusive):
           scale_i = xsec_i / sum_gen_weight_i
       - LFV signals (signal_tc, signal_tu):
           scale_i = 1.0 / sum_gen_weight_i  (equal generator normalization)
  4. Sums the scaled histograms (preserving Sumw2 variances for Kish's N_eff)
     into raw_<outfile>.root and runs postprocess_btag_eff.py (Option 2-3 smoothing).

Usage:
  python compute_btag_efficiency_rdataframe.py -Y 2022 -C muon -P ttbar \
      -D /path/to/mc/TTto2L2Nu /path/to/mc/TTtoLNu2Q /path/to/mc/TTto4Q \
      -O data/BTV/muon/2022/btag_eff_ttbar.root
"""

import sys
import os
import re
import glob
import argparse
from collections import OrderedDict

from run_all_btag_efficiency import (
    SM_PROCESS_GROUPS,
    SIGNAL_PROCESS_GROUPS,
    should_exclude_sample,
    classify_sample_dir,
)

HIST_NAMES = [
    "h2_denom_b",
    "h2_num_b",
    "h2_denom_c",
    "h2_num_c",
    "h2_denom_light",
    "h2_num_light",
    "h_nevents",
    "h_ncleanjetspass",
]

# ==============================================================================
# Authoritative SM Cross-Section Tables (in pb)
# Synchronized with plotIt/configs/Run3_{muon,electron}/files.yml (2022-2023BPix)
# and files24.yml (2024), including low-pT QCD samples.
# ==============================================================================
XSEC_2022_2023 = {
    # ttbar
    "TTtoLNu2Q": 405.75,
    "TTto2L2Nu": 98.04,
    "TTto4Q": 419.82,
    # Drell-Yan
    "DYto2L-2Jets_MLL-10to50": 20950.0,
    "DYto2L-2Jets_MLL-50": 6688.0,
    # Single Top
    "TBbar_s-channel": 2.278,
    "TbarB_s-channel": 1.43,
    "TQbarto2Q-t-channel": 87.75,
    "TQbartoLNu-t-channel": 42.25,
    "TbarQto2Q-t-channel": 53.028,
    "TbarQtoLNu-t-channel": 25.532,
    "TWminusto2L2Nu": 4.0,
    "TWminusto4Q": 16.0,
    "TWminustoLNu2Q": 16.0,
    "TbarWplusto2L2Nu": 4.01,
    "TbarWplusto4Q": 16.02,
    "TbarWplustoLNu2Q": 16.02,
    # tt+X
    "TTWtoQQ": 0.4678,
    "TTZtoQQ": 0.6603,
    "TTHto2B": 0.5742,
    "TTHtoNon2B": 0.5742,
    # W+Jets
    "WtoLNu-2Jets_0J": 55760.0,
    "WtoLNu-2Jets_1J": 9529.0,
    "WtoLNu-2Jets_2J": 3532.0,
    # Diboson
    "WW": 80.23,
    "WZ": 29.1,
    "ZZ": 13.15,
    # QCD MuEnriched (including Pt15To20)
    "QCD_Pt15To20_MuEnriched": 2982000.0,
    "QCD_Pt20To30_MuEnriched": 2679000.0,
    "QCD_Pt30To50_MuEnriched": 1497000.0,
    "QCD_Pt50To80_MuEnriched": 409500.0,
    "QCD_Pt80To120_MuEnriched": 95130.0,
    "QCD_Pt120To170_MuEnriched": 22980.0,
    "QCD_Pt170To300_MuEnriched": 7763.0,
    "QCD_Pt300To470_MuEnriched": 699.1,
    "QCD_Pt470To600_MuEnriched": 68.24,
    "QCD_Pt600To800_MuEnriched": 21.37,
    "QCD_Pt800To1000_MuEnriched": 3.913,
    "QCD_Pt1000_MuEnriched": 1.323,
    # QCD EMEnriched (including Pt10To30)
    "QCD_Pt10To30_EMEnriched": 6854000.0,
    "QCD_Pt30To50_EMEnriched": 6686000.0,
    "QCD_Pt50To80_EMEnriched": 2118000.0,
    "QCD_Pt80To120_EMEnriched": 394100.0,
    "QCD_Pt120To170_EMEnriched": 72170.0,
    "QCD_Pt170To300_EMEnriched": 17920.0,
    "QCD_Pt300_EMEnriched": 1224.0,
}

XSEC_2024 = {
    # ttbar
    "TTtoLNu2Q": 407.4,
    "TTto2L2Nu": 98.88,
    "TTto4Q": 419.82,
    # Drell-Yan (flavor-split in 2024)
    "DYto2E-2Jets_MLL-10-50": 7046.7,
    "DYto2E-2Jets_MLL-50": 2124.08,
    "DYto2Mu-2Jets_MLL-10-50": 7063.3,
    "DYto2Mu-2Jets_MLL-50": 2124.08,
    "DYto2Tau-2Jets_MLL-10-50": 7053.3,
    "DYto2Tau-2Jets_MLL-50": 2124.08,
    # Single Top
    "TBbarto2Q_s-channel": 4.065,
    "TBbartoLNu_s-channel": 2.033,
    "TbarBto2Q_s-channel": 2.55,
    "TbarBtoLNu_s-channel": 1.275,
    "TBbarQto2Q_t-channel": 77.26,
    "TBbarQtoLNu_t-channel": 38.6,
    "TbarBQto2Q_t-channel": 46.73,
    "TbarBQtoLNu_t-channel": 23.34,
    "TWminusto2L2Nu": 4.0,
    "TWminusto4Q": 16.0,
    "TWminustoLNu2Q": 16.0,
    "TbarWplusto2L2Nu": 4.01,
    "TbarWplusto4Q": 16.02,
    "TbarWplustoLNu2Q": 16.02,
    # tt+X
    "TTWtoQQ": 0.4678,
    "TTZtoQQ": 0.6337,
    "TTHto2B": 0.5742,
    "TTHtoNon2B": 0.5742,
    # W+Jets (flavor-split in 2024)
    "WtoENu-4Jets": 21132.0,
    "WtoMuNu-4Jets": 21132.0,
    "WtoTauNu-4Jets": 21449.0,
    # Diboson
    "WW": 80.23,
    "WZ": 29.1,
    "ZZ": 13.15,
    # QCD MuEnriched
    "QCD_Pt15To20_MuEnriched": 3018000.0,
    "QCD_Pt20To30_MuEnriched": 2701000.0,
    "QCD_Pt30To50_MuEnriched": 1461000.0,
    "QCD_Pt50To80_MuEnriched": 407600.0,
    "QCD_Pt80To120_MuEnriched": 26070.0,
    "QCD_Pt120To170_MuEnriched": 23140.0,
    "QCD_Pt170To300_MuEnriched": 7754.0,
    "QCD_Pt300To470_MuEnriched": 699.6,
    "QCD_Pt470To600_MuEnriched": 67.67,
    "QCD_Pt600To800_MuEnriched": 21.27,
    "QCD_Pt800To1000_MuEnriched": 3.89,
    "QCD_Pt1000_MuEnriched": 1.323,
    # QCD EMEnriched
    "QCD_Pt15To20_EMEnriched": 1444000.0,
    "QCD_Pt20To30_EMEnriched": 5309000.0,
    "QCD_Pt30To50_EMEnriched": 6849000.0,
    "QCD_Pt50To80_EMEnriched": 2130000.0,
    "QCD_Pt80To120_EMEnriched": 391400.0,
    "QCD_Pt120To170_EMEnriched": 71630.0,
    "QCD_Pt170To300_EMEnriched": 18010.0,
    "QCD_Pt300To470_EMEnriched": 1116.0,
    "QCD_Pt470To600_EMEnriched": 82.76,
    "QCD_Pt600To800_EMEnriched": 21.62,
    "QCD_Pt800To1000_EMEnriched": 3.361,
    "QCD_Pt1000_EMEnriched": 1.02,
}


def _parse_plotit_xsecs(year, channel):
    """Fallback parser reading plotIt/configs/Run3_<channel>/files*.yml if present."""
    this_dir = os.path.dirname(os.path.abspath(__file__))
    yml_name = "files24.yml" if "2024" in str(year) else "files.yml"
    yml_path = os.path.abspath(
        os.path.join(this_dir, "..", "..", "plotIt", "configs", f"Run3_{channel}", yml_name)
    )
    out = {}
    if not os.path.isfile(yml_path):
        return out

    cur_sample = None
    with open(yml_path, "r") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            m_file = re.match(r"^'hist_([^']+)\.root'\s*:", line)
            if m_file:
                cur_sample = m_file.group(1)
                continue
            if cur_sample:
                m_xs = re.match(r"^cross-section\s*:\s*([0-9.eE+-]+)", line)
                if m_xs:
                    out[cur_sample] = float(m_xs.group(1))
    return out


def get_sample_xsec(sample_name, year, channel, process_group):
    """
    Returns the sample normalization numerator:
      - 1.0 for signal_tc / signal_tu samples (equal generator-weight normalization).
      - Physical cross section (pb) for SM samples.
    Raises RuntimeError with [ERROR] prefix if missing or invalid.
    """
    sample_grp = classify_sample_dir(sample_name)
    if sample_grp in SIGNAL_PROCESS_GROUPS or process_group in ("signal", "signal_tc", "signal_tu"):
        if process_group in SM_PROCESS_GROUPS or process_group == "inclusive":
            raise RuntimeError(
                f"[ERROR] Signal sample '{sample_name}' must not be included in SM process group '{process_group}'."
            )
        return 1.0

    table = XSEC_2024 if "2024" in str(year) else XSEC_2022_2023
    xsec = table.get(sample_name)
    if xsec is None:
        dynamic_table = _parse_plotit_xsecs(year, channel)
        xsec = dynamic_table.get(sample_name)

    if xsec is None or xsec <= 0.0:
        raise RuntimeError(
            f"[ERROR] Missing or non-positive cross section for sample '{sample_name}' "
            f"(era='{year}', channel='{channel}', process='{process_group}')."
        )
    return float(xsec)


def _collect_root_files_in_dir(d):
    """Recursively collects all .root files under directory d."""
    out = []
    for root, _, filenames in os.walk(d):
        for fn in sorted(filenames):
            if fn.endswith(".root"):
                out.append(os.path.join(root, fn))
    return out


def collect_sample_groups(indirs, infiles, year, channel):
    """
    Groups input ROOT files by MC sample name: OrderedDict[sample_name -> list[filepath]].
    Applies should_exclude_sample() to filter out systematic variations ('__'),
    superseded inclusive/sliced samples, and channel-mismatched samples.
    """
    sample_map = OrderedDict()

    def _add_file(sample_name, fpath):
        if should_exclude_sample(sample_name, year, channel):
            return
        if sample_name not in sample_map:
            sample_map[sample_name] = []
        if fpath not in sample_map[sample_name]:
            sample_map[sample_name].append(fpath)

    if indirs:
        for d in indirs:
            if not os.path.exists(d):
                continue
            base = os.path.basename(d.rstrip("/"))
            if should_exclude_sample(base, year, channel):
                print(f"[INFO] Skipping excluded sample directory: {base}")
                continue
            if classify_sample_dir(base) != "unclassified":
                for f in _collect_root_files_in_dir(d):
                    _add_file(base, f)
            else:
                # Parent directory (e.g. era directory): inspect immediate subdirectories
                subdirs = sorted(os.listdir(d))
                found_subdir = False
                for sub in subdirs:
                    sub_path = os.path.join(d, sub)
                    if not os.path.isdir(sub_path):
                        continue
                    if should_exclude_sample(sub, year, channel):
                        continue
                    if classify_sample_dir(sub) != "unclassified":
                        found_subdir = True
                        for f in _collect_root_files_in_dir(sub_path):
                            _add_file(sub, f)
                if not found_subdir:
                    # Fallback for custom test directory name
                    for f in _collect_root_files_in_dir(d):
                        _add_file(base, f)

    if infiles:
        for pat in infiles:
            matched = [pat] if (os.path.isfile(pat) and pat.endswith(".root")) else [
                m for m in sorted(glob.glob(pat)) if m.endswith(".root") and os.path.isfile(m)
            ]
            for f in matched:
                parent_name = os.path.basename(os.path.dirname(os.path.abspath(f)))
                if classify_sample_dir(parent_name) != "unclassified" or should_exclude_sample(parent_name, year, channel):
                    sname = parent_name
                else:
                    fname = os.path.splitext(os.path.basename(f))[0]
                    sname = fname[5:] if fname.startswith("hist_") else fname
                _add_file(sname, f)

    # Remove any sample group that ended up with 0 files
    return OrderedDict((k, v) for k, v in sample_map.items() if v)


def main():
    parser = argparse.ArgumentParser(description="C++ RDataFrame B-tag Efficiency Calculator")
    parser.add_argument("-Y", "--year", dest="year", required=True, help="Year / Era (e.g. 2022, 2022EE, 2023, 2023BPix, 2024)")
    parser.add_argument("-C", "--channel", dest="channel", default="muon", choices=["muon", "electron"], help="Analysis channel (muon or electron)")
    parser.add_argument("-P", "--process", dest="process", default="ttbar", help="Process group name (e.g. ttbar, singletop, wjets, qcd, signal_tc, signal_tu, inclusive)")
    parser.add_argument("-D", "--indir", dest="indirs", nargs="+", default=[], help="Input directory or directories (searched recursively)")
    parser.add_argument("-I", "--infile", dest="infiles", nargs="+", default=[], help="Input ROOT file(s) or glob pattern(s)")
    parser.add_argument("-O", "--outfile", dest="outfile", default="", help="Final output ROOT file path")
    parser.add_argument("-j", "--threads", dest="threads", type=int, default=1, help="Number of implicit MT threads (default: 1)")
    parser.add_argument("--keep-raw", dest="keep_raw", action="store_true", help="Keep raw analyzer output ROOT file")
    parser.add_argument("--ref-file", dest="ref_file", default="", help="SM reference efficiency ROOT file for low-stat bin replacement")
    parser.add_argument("--min-neff", dest="min_neff", type=float, default=20.0, help="Minimum effective denominator entries per bin [default: 20.0]")
    args = parser.parse_args()

    # Determine default output file path if not given
    if not args.outfile:
        if args.process == "inclusive":
            args.outfile = f"data/BTV/{args.channel}/{args.year}/btag_eff.root"
        else:
            args.outfile = f"data/BTV/{args.channel}/{args.year}/btag_eff_{args.process}.root"

    out_dir = os.path.dirname(os.path.abspath(args.outfile))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    raw_outfile = os.path.join(out_dir, f"raw_{os.path.basename(args.outfile)}")

    # Collect files grouped by sample
    sample_groups = collect_sample_groups(args.indirs, args.infiles, args.year, args.channel)
    total_files = sum(len(flist) for flist in sample_groups.values())
    if not sample_groups or total_files == 0:
        print(f"[ERROR] No valid ROOT files found in specified inputs: dirs={args.indirs}, files={args.infiles}", file=sys.stderr)
        sys.exit(1)

    print("=" * 70)
    print(f"[BTagEfficiencyAnalyzer] Launching C++ RDataFrame job")
    print(f"  Era       : {args.year}")
    print(f"  Channel   : {args.channel}")
    print(f"  Process   : {args.process}")
    print(f"  Samples   : {len(sample_groups)} sample(s) ({total_files} ROOT files collected)")
    print(f"  Raw Out   : {raw_outfile}")
    print(f"  Final Out : {args.outfile}")
    print("=" * 70)

    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.PyConfig.IgnoreCommandLineOptions = True

    # Enable Implicit Multi-threading if threads > 1
    if args.threads > 1:
        ROOT.EnableImplicitMT(args.threads)

    # Load compiled shared library
    this_dir = os.path.dirname(os.path.abspath(__file__))
    lib_path = os.path.abspath(os.path.join(this_dir, "..", "libnanoadrdframe.so"))
    if not os.path.isfile(lib_path):
        alt_lib = os.path.abspath(os.path.join(this_dir, "..", "..", "LFVAnalyzer", "libnanoadrdframe.so"))
        if os.path.isfile(alt_lib):
            lib_path = alt_lib

    if os.path.isfile(lib_path):
        load_res = ROOT.gSystem.Load(lib_path)
        if load_res < 0:
            print(f"[ERROR] Failed to load shared library {lib_path}", file=sys.stderr)
            sys.exit(1)
    else:
        print(f"[WARNING] libnanoadrdframe.so not found at {lib_path}. Attempting ROOT autoloading...")

    # Pre-validate cross sections and hcounter sums for all samples before running
    sample_plan = []
    total_entries_all = 0
    for sample_name, sfiles in sample_groups.items():
        try:
            xsec = get_sample_xsec(sample_name, args.year, args.channel, args.process)
        except RuntimeError as e:
            print(str(e), file=sys.stderr)
            sys.exit(1)

        sum_gen_weight = 0.0
        valid_files = []
        sample_entries = 0

        for f in sfiles:
            tf = ROOT.TFile.Open(f)
            if not tf or tf.IsZombie():
                print(f"[ERROR] Corrupt or unreadable ROOT file '{f}' in sample '{sample_name}'.", file=sys.stderr)
                sys.exit(1)

            hc = tf.Get("hcounter")
            if not hc:
                tf.Close()
                print(f"[ERROR] Missing 'hcounter' histogram in file '{f}' (sample '{sample_name}').", file=sys.stderr)
                sys.exit(1)
            sum_gen_weight += float(hc.GetBinContent(2))

            tt = tf.Get("Events")
            if tt and tt.GetEntries() > 0:
                valid_files.append(f)
                sample_entries += int(tt.GetEntries())
            tf.Close()

        if sum_gen_weight <= 0.0:
            print(
                f"[ERROR] Non-positive generator weight sum (hcounter bin 2 = {sum_gen_weight}) "
                f"for sample '{sample_name}'.",
                file=sys.stderr,
            )
            sys.exit(1)

        scale_factor = xsec / sum_gen_weight
        total_entries_all += sample_entries
        sample_plan.append({
            "name": sample_name,
            "files": valid_files,
            "entries": sample_entries,
            "xsec": xsec,
            "sum_gen_weight": sum_gen_weight,
            "scale": scale_factor,
        })
        print(
            f"  [Sample] {sample_name:32s} | files={len(valid_files):3d}/{len(sfiles):3d} | "
            f"entries={sample_entries:9d} | xsec={xsec:12.6g} | sumW={sum_gen_weight:12.6g} | scale={scale_factor:.6e}"
        )

    if total_entries_all == 0:
        print(f"[ERROR] Total entries across all samples in process '{args.process}' is 0. Exiting.", file=sys.stderr)
        sys.exit(1)

    # Run BTagEfficiencyAnalyzer per sample and accumulate scaled histograms
    merged_hists = OrderedDict()
    raw_stem = os.path.splitext(os.path.basename(args.outfile))[0]

    for item in sample_plan:
        sname = item["name"]
        sfiles = item["files"]
        sentries = item["entries"]
        sscale = item["scale"]

        if sentries == 0 or not sfiles:
            print(f"[INFO] Sample '{sname}' has 0 skimmed Events entries; skipping RDataFrame execution.")
            continue

        sample_chain = ROOT.TChain("Events")
        for f in sfiles:
            sample_chain.Add(f)

        # Include sample name in temporary output path so CommonTools/src/impl_Reweight.cpp
        # accurately identifies TTto2L2Nu / TTtoLNu2Q / TTto4Q for TopPtWeight application.
        sample_raw_out = os.path.join(out_dir, f"raw_{raw_stem}_{sname}.root")

        try:
            analyzer = ROOT.BTagEfficiencyAnalyzer(
                sample_chain,
                sample_raw_out,
                args.year,
                args.channel,
                "nosyst",  # syst
                "",        # jsonfname
                False,     # applytauFF
                "",        # globaltag
                args.threads
            )
        except AttributeError as e:
            print(f"[ERROR] BTagEfficiencyAnalyzer not found in ROOT dictionary! Did you recompile with 'make'? Error: {e}", file=sys.stderr)
            sys.exit(1)

        analyzer.setupAnalysis()
        analyzer.run(False, "Events")

        tf_out = ROOT.TFile.Open(sample_raw_out, "READ")
        if not tf_out or tf_out.IsZombie():
            print(f"[ERROR] Failed to open analyzer output '{sample_raw_out}' for sample '{sname}'.", file=sys.stderr)
            sys.exit(1)

        for hname in HIST_NAMES:
            h = tf_out.Get(hname)
            if not h:
                tf_out.Close()
                print(f"[ERROR] Missing histogram '{hname}' in analyzer output '{sample_raw_out}'.", file=sys.stderr)
                sys.exit(1)
            if h.GetSumw2N() == 0:
                h.Sumw2()
            h.Scale(sscale)
            if hname not in merged_hists:
                h_clone = h.Clone(hname)
                h_clone.SetDirectory(0)
                merged_hists[hname] = h_clone
            else:
                merged_hists[hname].Add(h)

        tf_out.Close()
        if not args.keep_raw and os.path.isfile(sample_raw_out):
            try:
                os.remove(sample_raw_out)
            except OSError:
                pass

    if not merged_hists:
        print(f"[ERROR] No histograms were produced for process '{args.process}'.", file=sys.stderr)
        sys.exit(1)

    # Write combined scaled histograms to raw_outfile
    tf_merged = ROOT.TFile.Open(raw_outfile, "RECREATE")
    if not tf_merged or tf_merged.IsZombie():
        print(f"[ERROR] Cannot create merged raw output file '{raw_outfile}'.", file=sys.stderr)
        sys.exit(1)
    tf_merged.cd()
    for hname, hobj in merged_hists.items():
        hobj.Write(hname)
    tf_merged.Close()
    print(f">> C++ RDataFrame per-sample processing and cross-section scaling finished. Merged raw file: {raw_outfile}")

    # Post-process raw output into final efficiency maps
    from postprocess_btag_eff import process_file
    success = process_file(raw_outfile, args.outfile, ref_path=args.ref_file, min_neff=args.min_neff)
    if not success:
        print(f"[ERROR] Post-processing failed for {raw_outfile}", file=sys.stderr)
        sys.exit(1)

    if not args.keep_raw and os.path.isfile(raw_outfile):
        try:
            os.remove(raw_outfile)
        except OSError:
            pass

    print(f">> [SUCCESS] B-tagging efficiency map finalized at: {args.outfile}")

if __name__ == "__main__":
    main()

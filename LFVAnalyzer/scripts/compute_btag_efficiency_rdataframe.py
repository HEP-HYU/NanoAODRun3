#!/usr/bin/env python3
"""
compute_btag_efficiency_rdataframe.py

Executes the C++ BTagEfficiencyAnalyzer class via ROOT RDataFrame to extract
per-flavor 2D b-tagging efficiency maps at native C++ speed.
Supports recursive directory walking, multi-threading, and automated smoothing post-processing.

Usage:
  python compute_btag_efficiency_rdataframe.py -Y 2022 -C muon -P ttbar \
      -D /path/to/mc/dir1 /path/to/mc/dir2 -O data/BTV/muon/2022/btag_eff_ttbar.root
"""

import sys
import os
import glob
import argparse
import subprocess

def collect_root_files(indirs, infiles):
    """Recursively collects all .root files from directory trees and file lists."""
    files = []
    if infiles:
        for f in infiles:
            if os.path.isfile(f) and f.endswith(".root"):
                files.append(f)
            else:
                matched = glob.glob(f)
                files.extend([m for m in matched if m.endswith(".root") and os.path.isfile(m)])

    if indirs:
        for d in indirs:
            if not os.path.exists(d):
                continue
            for root, _, filenames in os.walk(d):
                for fn in filenames:
                    if fn.endswith(".root"):
                        files.append(os.path.join(root, fn))

    # Remove duplicates preserving order
    seen = set()
    unique_files = []
    for f in files:
        if f not in seen:
            seen.add(f)
            unique_files.append(f)
    return unique_files

def main():
    parser = argparse.ArgumentParser(description="C++ RDataFrame B-tag Efficiency Calculator")
    parser.add_argument("-Y", "--year", dest="year", required=True, help="Year / Era (e.g. 2022, 2022EE, 2023, 2023BPix, 2024)")
    parser.add_argument("-C", "--channel", dest="channel", default="muon", choices=["muon", "electron"], help="Analysis channel (muon or electron)")
    parser.add_argument("-P", "--process", dest="process", default="ttbar", help="Process group name (e.g. ttbar, singletop, wjets, qcd, inclusive)")
    parser.add_argument("-D", "--indir", dest="indirs", nargs="+", default=[], help="Input directory or directories (searched recursively)")
    parser.add_argument("-I", "--infile", dest="infiles", nargs="+", default=[], help="Input ROOT file(s) or glob pattern(s)")
    parser.add_argument("-O", "--outfile", dest="outfile", default="", help="Final output ROOT file path")
    parser.add_argument("-j", "--threads", dest="threads", type=int, default=1, help="Number of implicit MT threads (default: 1)")
    parser.add_argument("--keep-raw", dest="keep_raw", action="store_true", help="Keep raw analyzer output ROOT file")
    args = parser.parse_args()

    # Determine default output file path if not given
    if not args.outfile:
        args.outfile = f"data/BTV/{args.channel}/{args.year}/btag_eff_{args.process}.root"

    out_dir = os.path.dirname(os.path.abspath(args.outfile))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    raw_outfile = os.path.join(out_dir, f"raw_{os.path.basename(args.outfile)}")

    # Collect files
    files = collect_root_files(args.indirs, args.infiles)
    if not files:
        print(f"[ERROR] No ROOT files found in specified inputs: dirs={args.indirs}, files={args.infiles}", file=sys.stderr)
        sys.exit(1)

    print("=" * 70)
    print(f"[BTagEfficiencyAnalyzer] Launching C++ RDataFrame job")
    print(f"  Era       : {args.year}")
    print(f"  Channel   : {args.channel}")
    print(f"  Process   : {args.process}")
    print(f"  Files     : {len(files)} ROOT files collected")
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
        # Also try parent or LD_LIBRARY_PATH
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

    # Build TChain
    chain = ROOT.TChain("Events")
    valid_count = 0
    for f in files:
        tf = ROOT.TFile.Open(f)
        if tf and not tf.IsZombie():
            tt = tf.Get("Events")
            if tt and tt.GetEntries() > 0:
                chain.Add(f)
                valid_count += 1
            tf.Close()

    total_entries = chain.GetEntries()
    print(f">> Total valid files added: {valid_count}/{len(files)}")
    print(f">> Total entries to process: {total_entries}")

    if total_entries == 0:
        print("[ERROR] Total entries is 0. Exiting.", file=sys.stderr)
        sys.exit(1)

    # Instantiate BTagEfficiencyAnalyzer
    try:
        analyzer = ROOT.BTagEfficiencyAnalyzer(
            chain,
            raw_outfile,
            args.year,
            args.channel,
            "nosyst",  # syst
            "",        # jsonfname
            False,     # applytauFF
            "",        # globaltag
            args.threads
        )
    except AttributeError as e:
        print(f"[FATAL] BTagEfficiencyAnalyzer not found in ROOT dictionary! Did you recompile with 'make'? Error: {e}", file=sys.stderr)
        sys.exit(1)

    # Run analysis
    analyzer.setupAnalysis()
    analyzer.run(False, "Events")
    print(f">> C++ RDataFrame processing finished successfully. Raw file: {raw_outfile}")

    # Post-process raw output into final efficiency maps
    from postprocess_btag_eff import process_file
    success = process_file(raw_outfile, args.outfile)
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

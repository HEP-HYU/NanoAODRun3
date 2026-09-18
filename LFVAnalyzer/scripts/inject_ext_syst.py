#!/usr/bin/env python3
"""
inject_ext_syst.py
------------------
Injects systematic variation histograms from external ROOT files
(e.g., hist_TTto2L2Nu__hdampup.root, hist_TTto2L2Nu__tunedown.root)
into the corresponding nominal ROOT file (e.g., hist_TTto2L2Nu.root)
with the standard CMS variation naming convention: <histname>__<syst><var>.

This allows downstream tools such as plotIt and prepareShapesAndCards.py
to access external shape variations directly from nominal files.
"""

import os
import sys
import re
import argparse
import ROOT

ROOT.gROOT.SetBatch(True)


def parse_args():
    parser = argparse.ArgumentParser(description="Inject external systematic histograms into nominal ROOT files.")
    parser.add_argument("-d", "--dir", dest="directories", action="append", required=True,
                        help="Directory containing ROOT files (can be specified multiple times)")
    parser.add_argument("--systs", nargs="+", default=["hdamp", "tune"],
                        help="List of external systematics to inject (default: hdamp tune)")
    parser.add_argument("--filter", dest="filter_pattern", default="",
                        help="Optional regex pattern to filter files (e.g., 'TTto')")
    parser.add_argument("--dry-run", dest="dry_run", action="store_true", default=False,
                        help="Only show what would be done without modifying files")
    return parser.parse_args()


def inject_file(nom_filepath, syst_filepath, syst_tag, dry_run=False):
    """
    Reads histograms from syst_filepath and writes them into nom_filepath
    with the suffix '__' + syst_tag.
    """
    if not os.path.isfile(syst_filepath):
        print(f"  [SKIP] External file not found: {syst_filepath}")
        return 0

    if not os.path.isfile(nom_filepath):
        print(f"  [ERROR] Nominal file not found: {nom_filepath}")
        return 0

    f_syst = ROOT.TFile.Open(syst_filepath, "READ")
    if not f_syst or f_syst.IsZombie():
        print(f"  [ERROR] Failed to open systematic file: {syst_filepath}")
        return 0

    keys = [k.GetName() for k in f_syst.GetListOfKeys()]
    # Remove duplicate keys if any
    keys = list(dict.fromkeys(keys))

    # Identify histograms to inject (exclude counters and existing variation histograms)
    skip_keywords = ["counter", "Sum"]
    hists_to_inject = []
    for k in keys:
        if any(skip in k for skip in skip_keywords):
            continue
        if "__" in k:
            continue
        obj = f_syst.Get(k)
        if obj and obj.InheritsFrom("TH1"):
            hists_to_inject.append(k)

    if not hists_to_inject:
        print(f"  [INFO] No histograms to inject in {os.path.basename(syst_filepath)}")
        f_syst.Close()
        return 0

    if dry_run:
        print(f"  [DRY-RUN] Would inject {len(hists_to_inject)} histograms from "
              f"{os.path.basename(syst_filepath)} into {os.path.basename(nom_filepath)}")
        f_syst.Close()
        return len(hists_to_inject)

    f_nom = ROOT.TFile.Open(nom_filepath, "UPDATE")
    if not f_nom or f_nom.IsZombie():
        print(f"  [ERROR] Failed to open nominal file for update: {nom_filepath}")
        f_syst.Close()
        return 0

    count = 0
    f_nom.cd()
    for hname in hists_to_inject:
        h_src = f_syst.Get(hname)
        if not h_src:
            continue
        target_name = f"{hname}__{syst_tag}"
        h_dest = h_src.Clone(target_name)
        h_dest.SetDirectory(f_nom)
        h_dest.Write(target_name, ROOT.TObject.kOverwrite)
        count += 1

    f_nom.Close()
    f_syst.Close()

    print(f"  [DONE] Injected {count} histograms into {os.path.basename(nom_filepath)} (tag: __{syst_tag})")
    return count


def process_directory(directory, systs, filter_pattern, dry_run=False):
    if not os.path.isdir(directory):
        print(f"[WARN] Directory does not exist: {directory}")
        return

    print(f"\n=======================================================")
    print(f"Scanning directory: {directory}")
    print(f"=======================================================")

    all_files = os.listdir(directory)
    total_injected = 0

    # Build regex to match: hist_<process>__<syst><var>.root
    syst_regex = re.compile(r"^(hist_.*)__(" + "|".join(systs) + r")(up|down)\.root$")

    # Group external files by nominal file
    matched_pairs = []
    for fname in sorted(all_files):
        m = syst_regex.match(fname)
        if m:
            nom_base = m.group(1)
            syst_name = m.group(2)
            var_name = m.group(3)
            syst_tag = f"{syst_name}{var_name}"

            if filter_pattern and not re.search(filter_pattern, nom_base):
                continue

            nom_fname = f"{nom_base}.root"
            nom_path = os.path.join(directory, nom_fname)
            syst_path = os.path.join(directory, fname)
            matched_pairs.append((nom_path, syst_path, syst_tag))

    if not matched_pairs:
        print(f"No matching external systematic files found in {directory}.")
        return

    print(f"Found {len(matched_pairs)} systematic files to inject.")
    for nom_path, syst_path, syst_tag in matched_pairs:
        n = inject_file(nom_path, syst_path, syst_tag, dry_run=dry_run)
        total_injected += n

    print(f"Summary for {directory}: Total {total_injected} histograms injected.\n")


def main():
    args = parse_args()
    for d in args.directories:
        process_directory(d, args.systs, args.filter_pattern, dry_run=args.dry_run)


if __name__ == "__main__":
    main()

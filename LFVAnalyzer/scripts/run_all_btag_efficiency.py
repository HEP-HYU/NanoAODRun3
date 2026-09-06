#!/usr/bin/env python3
"""
run_all_btag_efficiency.py

Batch runner to discover all MC sample directories under the skimmed tree,
group them by process (ttbar, singletop, wjets, dyjets, qcd, diboson, ttx, signal, inclusive),
and compute b-tagging efficiency maps for muon and/or electron channels across eras.

Supports two engines:
  - 'rdataframe' (default): Uses high-performance C++ BTagEfficiencyAnalyzer via ROOT RDataFrame
  - 'python': Uses standalone compute_btag_efficiency.py

Usage examples:
  # Run both muon and electron channels across all eras using C++ RDataFrame (Recommended)
  python run_all_btag_efficiency.py --channel both --engine rdataframe -j 4

  # Run only muon channel for 2022 ttbar with dry-run
  python run_all_btag_efficiency.py --channel muon -Y 2022 -P ttbar --dry-run

  # Run electron channel with python engine
  python run_all_btag_efficiency.py --channel electron --engine python -j 4
"""

import sys
import os
import re
import argparse
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

DEFAULT_BASE_DIR_TEMPLATE = "/data2/common/skimmed_NanoAOD/skim_0812_LFV/{channel}/mc/"

# Classification of process groups based on sample folder naming conventions
# Aligned with files.yml / files24.yml definitions
PROCESS_GROUPS = [
    "ttbar",
    "singletop",
    "wjets",
    "dyjets",
    "qcd",
    "diboson",
    "ttx",
    "signal",
    "inclusive",
]

# Standard BTV campaign name mappings aligned with analysis_config.json
ERA_TO_BTAG_YEAR = {
    "2022": "2022_Summer22",
    "2022EE": "2022_Summer22EE",
    "2023": "2023_Summer23",
    "2023BPix": "2023_Summer23BPix",
    "2024": "2024_Summer24",
}

def parse_era_from_dirname(dirname):
    """
    Parses era string from directory name.
    Matches:
      v12_2022     -> 2022
      v12_2022EE   -> 2022EE
      v12_2023     -> 2023
      v12_2023BPix -> 2023BPix
      v13_2024     -> 2024
    """
    dirname = os.path.basename(dirname.rstrip("/"))
    m = re.search(r'(?:v\d+_)?(202[2-4](?:EE|BPix|postEE|preEE)?)', dirname, re.IGNORECASE)
    if m:
        raw = m.group(1)
        # Normalize
        if raw.lower() == "2022postee":
            return "2022EE"
        if raw.lower() == "2022pree":
            return "2022"
        return raw
    return None

def classify_sample_dir(dirname):
    """Classifies a sample directory into one of the standard process groups."""
    name = os.path.basename(dirname.rstrip("/"))
    lname = name.lower()

    # 1. Signal (LFV / TCMuTau)
    if "tcmu" in lname or "tcmutau" in lname or "tce" in lname or "tcetau" in lname or "lfv" in lname:
        return "signal"

    # 2. TTX (ttW, ttZ, ttH, ttTT, etc.)
    if (lname.startswith("ttw") or lname.startswith("ttz") or
        lname.startswith("tth") or lname.startswith("tttt")):
        return "ttx"

    # 3. Pure ttbar (TTto2L2Nu, TTto4Q, TTtoLNu2Q, TT_TuneCP5)
    if (lname.startswith("ttto") or lname.startswith("tt_") or
        "ttto2l2nu" in lname or "tttolnu2q" in lname or "ttto4q" in lname):
        return "ttbar"

    # 4. Single Top (t-channel, s-channel, tW, tbarW, TBbarQ, etc.)
    if (lname.startswith("tw") or lname.startswith("tbarw") or
        "t-channel" in lname or "s-channel" in lname or
        lname.startswith("tbbar") or lname.startswith("tq") or
        lname.startswith("st_")):
        return "singletop"

    # 5. Drell-Yan (DYto2L, DYJetsToLL)
    if lname.startswith("dy"):
        return "dyjets"

    # 6. W+Jets (WtoLNu, WJetsToLNu)
    if lname.startswith("wto") or lname.startswith("wjets") or lname.startswith("w_"):
        return "wjets"

    # 7. Diboson (WW, WZ, ZZ)
    if (lname.startswith("ww") or lname.startswith("wz") or lname.startswith("zz") or
        "_ww" in lname or "_wz" in lname or "_zz" in lname):
        return "diboson"

    # 8. QCD Multijet
    if lname.startswith("qcd"):
        return "qcd"

    return "unclassified"

def discover_jobs(base_dir, channel, target_years=None, target_procs=None, out_base="data/BTV"):
    """
    Scans base_dir for era subdirectories and sample directories,
    groups them by process, and returns a list of job specifications.
    """
    jobs = []
    if not os.path.isdir(base_dir):
        return []

    subdirs = sorted(os.listdir(base_dir))
    era_dirs = {}
    for sub in subdirs:
        full_sub = os.path.join(base_dir, sub)
        if not os.path.isdir(full_sub):
            continue
        era = parse_era_from_dirname(sub)
        if era:
            if target_years and era not in target_years:
                continue
            era_dirs[era] = full_sub

    for era, era_path in sorted(era_dirs.items()):
        samples = sorted(os.listdir(era_path))
        proc_samples = {p: [] for p in PROCESS_GROUPS if p != "inclusive"}
        all_samples = []

        for s in samples:
            sample_path = os.path.join(era_path, s)
            if not os.path.isdir(sample_path):
                continue
            grp = classify_sample_dir(s)
            if grp in proc_samples:
                proc_samples[grp].append(sample_path)
                all_samples.append(sample_path)
            else:
                print(f"[WARN] Unclassified sample directory in {channel}/{era}: {s}")

        # Build job specifications
        for proc in PROCESS_GROUPS:
            if target_procs and proc not in target_procs:
                continue

            if proc == "inclusive":
                dirs = all_samples
                outfile = os.path.join(out_base, channel, era, "btag_eff.root")
            else:
                dirs = proc_samples.get(proc, [])
                outfile = os.path.join(out_base, channel, era, f"btag_eff_{proc}.root")

            if not dirs:
                continue

            jobs.append({
                "channel": channel,
                "era": era,
                "process": proc,
                "input_dirs": dirs,
                "outfile": outfile,
            })

    return jobs

def run_single_job(job, script_path, engine, extra_args):
    """Executes the calculation script for a single job."""
    cmd = [
        sys.executable,
        script_path,
        "-Y", job["era"],
        "-C", job["channel"],
        "-P", job["process"],
        "-O", job["outfile"],
        "-D",
    ] + job["input_dirs"]

    if extra_args:
        cmd.extend(extra_args)

    tag = f"[{job['channel']} | {job['era']} | {job['process']}]"
    print(f"[START] {tag} -> {job['outfile']} ({len(job['input_dirs'])} dataset dirs) [Engine: {engine}]")
    try:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=True)
        print(f"[DONE ] {tag} successfully produced: {job['outfile']}")
        return True, job, proc.stdout
    except subprocess.CalledProcessError as e:
        print(f"[FAIL ] {tag} failed with code {e.returncode}!")
        if e.stdout:
            print(f"\n{'='*25} ERROR LOG: {tag} {'='*25}")
            print(e.stdout.strip())
            print(f"{'='*65}\n")
        return False, job, e.stdout

def create_era_symlinks(out_base, channels):
    """
    Creates bidirectional convenience symlinks between short era names (e.g. '2022')
    and BTV campaign names (e.g. '2022_Summer22') so that any analyzer can seamlessly
    find efficiency maps regardless of naming convention.
    """
    for ch in channels:
        ch_dir = os.path.join(out_base, ch)
        if not os.path.isdir(ch_dir):
            continue
        for short_era, btv_era in ERA_TO_BTAG_YEAR.items():
            short_path = os.path.join(ch_dir, short_era)
            btv_path = os.path.join(ch_dir, btv_era)
            if os.path.isdir(short_path) and not os.path.exists(btv_path):
                try:
                    os.symlink(short_era, btv_path)
                    print(f"[INFO] Created BTV era symlink: {ch}/{btv_era} -> {short_era}")
                except OSError:
                    pass
            elif os.path.isdir(btv_path) and not os.path.exists(short_path):
                try:
                    os.symlink(btv_era, short_path)
                    print(f"[INFO] Created short era symlink: {ch}/{short_era} -> {btv_era}")
                except OSError:
                    pass

def main():
    parser = argparse.ArgumentParser(description="Batch orchestrator to compute b-tag efficiencies across channels, eras, and processes")
    parser.add_argument("-C", "--channel", dest="channel", type=str, default="both",
                        choices=["muon", "electron", "both", "all"],
                        help="Analysis channel: 'muon', 'electron', or 'both'/'all' [default: both]")
    parser.add_argument("-B", "--base-dir", dest="base_dir", type=str, default="",
                        help="Override base directory for MC skims (supports '{channel}' placeholder). Default: /data2/common/skimmed_NanoAOD/skim_0812_LFV/{channel}/mc/")
    parser.add_argument("-E", "--engine", dest="engine", type=str, default="rdataframe",
                        choices=["rdataframe", "python"],
                        help="Calculation engine: 'rdataframe' (C++ BTagEfficiencyAnalyzer, recommended) or 'python' [default: rdataframe]")
    parser.add_argument("-Y", "--years", dest="years", nargs="*", default=[],
                        help="Specific eras to run (e.g. 2022 2022EE 2023 2023BPix 2024). Default: all found")
    parser.add_argument("-P", "--processes", dest="processes", nargs="*", default=[],
                        help=f"Specific processes to run (choices: {', '.join(PROCESS_GROUPS)}). Default: all")
    parser.add_argument("-O", "--out-base", dest="out_base", type=str, default="data/BTV",
                        help="Base output directory [default: data/BTV]")
    parser.add_argument("-j", "--jobs", dest="jobs", type=int, default=4,
                        help="Number of concurrent worker processes [default: 4]")
    parser.add_argument("--dry-run", dest="dry_run", action="store_true",
                        help="Display discovered jobs and command plan without executing")
    parser.add_argument("--skip-existing", dest="skip_existing", action="store_true",
                        help="Skip jobs whose output ROOT file already exists")
    parser.add_argument("--eta-bins", dest="eta_bins", type=str, default="single",
                        help="[Python engine only] Eta binning: 'single' (default), 'barrel-endcap', or comma-separated edges")
    parser.add_argument("--pt-bins", dest="pt_bins", type=str, default="",
                        help="[Python engine only] Custom comma-separated pT bin edges")
    parser.add_argument("--max-events", dest="max_events", type=int, default=-1,
                        help="Max events to process per job (-1 for all)")
    args = parser.parse_args()

    this_dir = os.path.dirname(os.path.abspath(__file__))

    # Determine execution script based on engine
    if args.engine == "rdataframe":
        compute_script = os.path.join(this_dir, "compute_btag_efficiency_rdataframe.py")
    else:
        compute_script = os.path.join(this_dir, "compute_btag_efficiency.py")

    if not os.path.isfile(compute_script):
        print(f"[ERROR] Engine script not found at {compute_script}!", file=sys.stderr)
        sys.exit(1)

    # Determine active channels
    if args.channel in ["both", "all"]:
        active_channels = ["muon", "electron"]
    else:
        active_channels = [args.channel]

    print("=" * 80)
    print("BTV Method 1a Automated Efficiency Map Batch Runner")
    print(f"Engine               : {args.engine.upper()} ({os.path.basename(compute_script)})")
    print(f"Channels             : {', '.join(active_channels)}")
    print(f"Output Base Directory: {args.out_base}")
    print(f"Target Eras          : {args.years if args.years else 'All discovered'}")
    print(f"Target Processes     : {args.processes if args.processes else 'All (ttbar, singletop, wjets, dyjets, qcd, diboson, ttx, signal, inclusive)'}")
    print(f"Concurrency (Workers): {args.jobs}")
    print("=" * 80)

    # Discover jobs across active channels
    all_jobs = []
    for ch in active_channels:
        if args.base_dir:
            if "{channel}" in args.base_dir:
                ch_base_dir = args.base_dir.format(channel=ch)
            else:
                ch_base_dir = args.base_dir
        else:
            ch_base_dir = DEFAULT_BASE_DIR_TEMPLATE.format(channel=ch)

        if not os.path.isdir(ch_base_dir):
            print(f"[WARN] Base directory for channel '{ch}' does not exist: {ch_base_dir}")
            continue

        ch_jobs = discover_jobs(
            base_dir=ch_base_dir,
            channel=ch,
            target_years=args.years,
            target_procs=args.processes,
            out_base=args.out_base
        )
        print(f"[INFO] Discovered {len(ch_jobs)} jobs for channel '{ch}' in {ch_base_dir}")
        all_jobs.extend(ch_jobs)

    if not all_jobs:
        print("[WARNING] No matching jobs found! Please check input directory paths or filters.")
        sys.exit(0)

    # Handle --skip-existing
    if args.skip_existing:
        filtered_jobs = []
        for j in all_jobs:
            if os.path.isfile(j["outfile"]):
                print(f"[SKIP] Output already exists: {j['outfile']}")
            else:
                filtered_jobs.append(j)
        all_jobs = filtered_jobs

    print(f"\n>> Total active jobs to run: {len(all_jobs)}")

    extra_args = []
    if args.engine == "python":
        extra_args.extend(["--eta-bins", args.eta_bins])
        if args.pt_bins:
            extra_args.extend(["--pt-bins", args.pt_bins])
    if args.max_events > 0:
        extra_args.extend(["--max-events", str(args.max_events)])

    if args.dry_run:
        print("\n" + "=" * 30 + " DRY RUN PLAN " + "=" * 30)
        for i, j in enumerate(all_jobs, 1):
            print(f"{i:>3}. [{j['channel']:^8} | {j['era']:^10} | {j['process']:^10}]")
            print(f"     Output: {j['outfile']}")
            print(f"     Inputs: {len(j['input_dirs'])} directories:")
            for d in j['input_dirs'][:3]:
                print(f"       - {os.path.basename(d)}")
            if len(j['input_dirs']) > 3:
                print(f"       ... and {len(j['input_dirs']) - 3} more")
        print("=" * 74)
        print("Dry run completed. No commands executed.")
        return

    # Execute jobs using ProcessPoolExecutor
    n_workers = min(args.jobs, len(all_jobs)) if all_jobs else 1
    print(f"\nStarting execution with {n_workers} concurrent workers...\n")

    results = []
    failed_jobs = []

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        future_to_job = {
            executor.submit(run_single_job, job, compute_script, args.engine, extra_args): job
            for job in all_jobs
        }

        for future in as_completed(future_to_job):
            success, job, _ = future.result()
            results.append((success, job))
            if not success:
                failed_jobs.append(job)
    print("EXECUTION SUMMARY")
    print(f"Total jobs scheduled : {len(all_jobs)}")
    print(f"Successfully finished: {len(all_jobs) - len(failed_jobs)}")
    print(f"Failed jobs          : {len(failed_jobs)}")

    # Create era symlinks for seamless analyzer loading
    create_era_symlinks(args.out_base, active_channels)

    if failed_jobs:
        print("\nFailed Jobs:")
        for fj in failed_jobs:
            print(f"  - [{fj['channel']} | {fj['era']} | {fj['process']}] -> {fj['outfile']}")
        print("=" * 80)
        sys.exit(1)
    else:
        print("\nAll b-tag efficiency maps were computed and saved successfully!")
        print("=" * 80)

if __name__ == "__main__":
    main()

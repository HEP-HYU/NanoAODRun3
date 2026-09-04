#!/usr/bin/env python3
"""
===============================================================================
run_all_btag_efficiency.py
-------------------------------------------------------------------------------
Orchestrator script to automatically discover skimmed MC sample directories
(under e.g. /data2/common/skimmed_NanoAOD/skim_0812_LFV/muon/mc/), classify them
into plotIt-aligned process groups (ttbar, singletop, wjets, dyjets, qcd,
diboson, ttx, signal, inclusive), and run compute_btag_efficiency.py for each
process across all eras (2022, 2022EE, 2023, 2023BPix, 2024).

Features:
  - Recursive folder scanning per era
  - 1:1 classification aligned with plotIt files.yml / files24.yml
  - Multi-process parallel execution (--jobs / -j)
  - Dry-run mode (--dry-run) to inspect discovery & execution plan
  - Skip existing outputs (--skip-existing)
  - Selective year/process filtering (--years, --processes)
===============================================================================
"""

import os
import sys
import argparse
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

DEFAULT_BASE_DIR = "/data2/common/skimmed_NanoAOD/skim_0812_LFV/muon/mc/"

ERA_FOLDERS = {
    "v12_2022": "2022",
    "v12_2022EE": "2022EE",
    "v12_2023": "2023",
    "v12_2023BPix": "2023BPix",
    "v15_2024": "2024",
}

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


def parse_era_from_dirname(dirname):
    """Map directory name (e.g. v12_2022, v12_2022EE, v15_2024) to era name."""
    clean = os.path.basename(os.path.normpath(dirname))
    if clean in ERA_FOLDERS:
        return ERA_FOLDERS[clean]
    for key, era in ERA_FOLDERS.items():
        if key in clean:
            return era
    # Fallback search by year string
    for era in ["2022EE", "2022", "2023BPix", "2023", "2024"]:
        if era in clean:
            return era
    return None


def classify_sample_dir(dirname):
    """
    Classify a sample directory name into one of the plotIt process groups:
      - signal    (LFV signals)
      - ttx       (TTW, TTZ, TTH)
      - ttbar     (TTto2L2Nu, TTtoLNu2Q, TTto4Q)
      - singletop (TBbar, TbarB, TQbar, TbarQ, TWminus, TbarWplus, etc.)
      - wjets     (WtoLNu, WtoENu, WtoMuNu, WtoTauNu)
      - dyjets    (DYto2L, DYto2E, DYto2Mu, DYto2Tau)
      - diboson   (WW, WZ, ZZ)
      - qcd       (QCD_Pt*)
    """
    low = os.path.basename(os.path.normpath(dirname)).lower()

    # 1. Signal (LFV) — check first to prevent TTto*LFV being caught as ttbar
    if "lfv" in low:
        return "signal"

    # 2. tt+X / ttV — check before ttbar
    if any(k in low for k in ["ttw", "ttz", "tth"]):
        return "ttx"

    # 3. ttbar
    if any(k in low for k in ["ttto", "ttbar", "tt_"]):
        return "ttbar"

    # 4. Single Top
    if any(k in low for k in ["tbbar", "tbarb", "tqbar", "tbarq", "twminus", "tbarwplus", "singletop", "st_"]):
        return "singletop"

    # 5. W+Jets (including 2024 flavor-split samples)
    if any(k in low for k in ["wtolnu", "wtoenu", "wtomunu", "wtotaunu", "wto", "wjets"]):
        return "wjets"

    # 6. Z+Jets / Drell-Yan (including 2024 flavor-split samples)
    if any(k in low for k in ["dyto2l", "dyto2e", "dyto2mu", "dyto2tau", "dyto", "dyjets", "zjets"]):
        return "dyjets"

    # 7. Diboson (WW, WZ, ZZ)
    if low in ["ww", "wz", "zz"] or any(k in low for k in ["ww_", "wz_", "zz_", "_ww", "_wz", "_zz", "diboson"]):
        return "diboson"

    # 8. QCD
    if "qcd" in low:
        return "qcd"

    return "unclassified"


def discover_jobs(base_dir, target_years=None, target_procs=None, out_base="data/BTV"):
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
                print(f"[WARN] Unclassified sample directory in {era}: {s}")

        # Build job specifications
        for proc in PROCESS_GROUPS:
            if target_procs and proc not in target_procs:
                continue

            if proc == "inclusive":
                dirs = all_samples
                outfile = os.path.join(out_base, era, "btag_eff.root")
            else:
                dirs = proc_samples.get(proc, [])
                outfile = os.path.join(out_base, era, f"btag_eff_{proc}.root")

            if not dirs:
                continue

            jobs.append({
                "era": era,
                "process": proc,
                "input_dirs": dirs,
                "outfile": outfile,
            })

    return jobs


def run_single_job(job, script_path, extra_args):
    """Executes compute_btag_efficiency.py for a single job."""
    cmd = [
        sys.executable,
        script_path,
        "-Y", job["era"],
        "-P", job["process"],
        "-O", job["outfile"],
        "-D",
    ] + job["input_dirs"]

    if extra_args:
        cmd.extend(extra_args)

    era_proc = f"[{job['era']} | {job['process']}]"
    print(f"[START] {era_proc} -> {job['outfile']} ({len(job['input_dirs'])} dataset dirs)")
    try:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=True)
        print(f"[DONE ] {era_proc} successfully produced: {job['outfile']}")
        return True, job, proc.stdout
    except subprocess.CalledProcessError as e:
        print(f"[FAIL ] {era_proc} failed with code {e.returncode}!")
        return False, job, e.stdout


def main():
    parser = argparse.ArgumentParser(description="Batch orchestrator to compute b-tag efficiencies across eras and processes")
    parser.add_argument("-B", "--base-dir", dest="base_dir", type=str, default=DEFAULT_BASE_DIR,
                        help=f"Base directory containing skimmed MC era folders [default: {DEFAULT_BASE_DIR}]")
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
    parser.add_argument("--channel", dest="channel", type=str, default="muon",
                        choices=["muon", "electron", "both"], help="Lepton selection channel [default: muon]")
    parser.add_argument("--eta-bins", dest="eta_bins", type=str, default="single",
                        help="Eta binning: 'single' (default), 'barrel-endcap', or comma-separated edges")
    parser.add_argument("--pt-bins", dest="pt_bins", type=str, default="",
                        help="Custom comma-separated pT bin edges")
    parser.add_argument("--max-events", dest="max_events", type=int, default=-1,
                        help="Max events to process per job (-1 for all)")
    args = parser.parse_args()

    this_dir = os.path.dirname(os.path.abspath(__file__))
    compute_script = os.path.join(this_dir, "compute_btag_efficiency.py")

    if not os.path.isfile(compute_script):
        print(f"ERROR: compute_btag_efficiency.py not found at {compute_script}!")
        sys.exit(1)

    print("=" * 80)
    print("BTV Method 1a Automated Efficiency Map Batch Runner")
    print(f"Base Skim Directory : {args.base_dir}")
    print(f"Output Base Directory: {args.out_base}")
    print(f"Target Eras          : {args.years if args.years else 'All discovered'}")
    print(f"Target Processes     : {args.processes if args.processes else 'All (ttbar, singletop, wjets, dyjets, qcd, diboson, ttx, signal, inclusive)'}")
    print(f"Concurrency (Workers): {args.jobs}")
    print(f"Eta Binning          : {args.eta_bins}")
    print("=" * 80)

    extra_args = [
        "--channel", args.channel,
        "--eta-bins", args.eta_bins,
    ]
    if args.pt_bins:
        extra_args.extend(["--pt-bins", args.pt_bins])
    if args.max_events > 0:
        extra_args.extend(["--max-events", str(args.max_events)])

    # Discover jobs from filesystem
    jobs = discover_jobs(
        base_dir=args.base_dir,
        target_years=args.years,
        target_procs=args.processes,
        out_base=args.out_base,
    )

    if not jobs:
        # Fallback to predefined topology if running on a machine without /data2 access
        print("\n[INFO] /data2 not locally accessible. Generating execution plan based on canonical dataset topology...")
        mock_samples_22_23 = [
            "DYto2L-2Jets_MLL-10to50", "DYto2L-2Jets_MLL-50",
            "QCD_Pt1000_MuEnriched", "QCD_Pt120To170_MuEnriched", "QCD_Pt15To20_MuEnriched",
            "QCD_Pt170To300_MuEnriched", "QCD_Pt20To30_MuEnriched", "QCD_Pt300To470_MuEnriched",
            "QCD_Pt30To50_MuEnriched", "QCD_Pt470To600_MuEnriched", "QCD_Pt50To80_MuEnriched",
            "QCD_Pt600To800_MuEnriched", "QCD_Pt800To1000_MuEnriched", "QCD_Pt80To120_MuEnriched",
            "TBbarQ_t-channel", "TBbar_s-channel", "TCMuTau-LFV-Scalar", "TCMuTau-LFV-Tensor", "TCMuTau-LFV-Vector",
            "TQbarto2Q-t-channel", "TQbartoLNu-t-channel", "TTHto2B", "TTHtoNon2B", "TTWtoQQ", "TTZtoQQ",
            "TTto2L2Nu", "TTto4Q", "TTtoCMuTau-LFV-Scalar", "TTtoCMuTau-LFV-Tensor", "TTtoCMuTau-LFV-Vector",
            "TTtoLNu2Q", "TTtoUMuTau-LFV-Scalar", "TTtoUMuTau-LFV-Tensor", "TTtoUMuTau-LFV-Vector",
            "TUMuTau-LFV-Scalar", "TUMuTau-LFV-Tensor", "TUMuTau-LFV-Vector", "TWminusto2L2Nu", "TWminusto4Q", "TWminustoLNu2Q",
            "TbarBQ_t-channel", "TbarB_s-channel", "TbarQto2Q-t-channel", "TbarQtoLNu-t-channel",
            "TbarWplusto2L2Nu", "TbarWplusto4Q", "TbarWplustoLNu2Q", "WW", "WZ", "ZZ",
            "WtoLNu-2Jets_0J", "WtoLNu-2Jets_1J", "WtoLNu-2Jets_2J"
        ]
        mock_samples_24 = [
            "DYto2E-2Jets_MLL-10-50", "DYto2E-2Jets_MLL-50", "DYto2Mu-2Jets_MLL-10-50", "DYto2Mu-2Jets_MLL-50",
            "DYto2Tau-2Jets_MLL-10-50", "DYto2Tau-2Jets_MLL-50",
            "QCD_Pt1000_MuEnriched", "QCD_Pt120To170_MuEnriched", "QCD_Pt15To20_MuEnriched",
            "QCD_Pt170To300_MuEnriched", "QCD_Pt20To30_MuEnriched", "QCD_Pt300To470_MuEnriched",
            "QCD_Pt30To50_MuEnriched", "QCD_Pt470To600_MuEnriched", "QCD_Pt50To80_MuEnriched",
            "QCD_Pt600To800_MuEnriched", "QCD_Pt800To1000_MuEnriched", "QCD_Pt80To120_MuEnriched",
            "TBbarQto2Q_t-channel", "TBbarQtoLNu_t-channel", "TBbarto2Q_s-channel", "TBbartoLNu_s-channel",
            "TCMuTau-LFV-Scalar", "TCMuTau-LFV-Tensor", "TCMuTau-LFV-Vector", "TTHto2B", "TTHtoNon2B", "TTWtoQQ", "TTZtoQQ",
            "TTto2L2Nu", "TTto4Q", "TTtoCMuTau-LFV-Scalar", "TTtoCMuTau-LFV-Tensor", "TTtoCMuTau-LFV-Vector",
            "TTtoLNu2Q", "TTtoUMuTau-LFV-Scalar", "TTtoUMuTau-LFV-Tensor", "TTtoUMuTau-LFV-Vector",
            "TUMuTau-LFV-Scalar", "TUMuTau-LFV-Tensor", "TUMuTau-LFV-Vector", "TWminusto2L2Nu", "TWminusto4Q", "TWminustoLNu2Q",
            "TbarBQto2Q_t-channel", "TbarBQtoLNu_t-channel", "TbarBto2Q_s-channel", "TbarBtoLNu_s-channel",
            "TbarWplusto2L2Nu", "TbarWplusto4Q", "TbarWplustoLNu2Q", "WW", "WZ", "ZZ",
            "WtoENu-4Jets", "WtoLNu-4Jets_1J", "WtoLNu-4Jets_2J", "WtoLNu-4Jets_3J", "WtoLNu-4Jets_4J",
            "WtoMuNu-4Jets", "WtoTauNu-4Jets"
        ]
        eras_spec = [
            ("2022", "v12_2022", mock_samples_22_23),
            ("2022EE", "v12_2022EE", mock_samples_22_23),
            ("2023", "v12_2023", mock_samples_22_23),
            ("2023BPix", "v12_2023BPix", mock_samples_22_23),
            ("2024", "v15_2024", mock_samples_24),
        ]
        for era, era_folder, smp_list in eras_spec:
            if args.years and era not in args.years:
                continue
            proc_map = {p: [] for p in PROCESS_GROUPS if p != "inclusive"}
            all_list = []
            for s in smp_list:
                s_path = os.path.join(args.base_dir, era_folder, s)
                grp = classify_sample_dir(s)
                if grp in proc_map:
                    proc_map[grp].append(s_path)
                    all_list.append(s_path)
            for proc in PROCESS_GROUPS:
                if args.processes and proc not in args.processes:
                    continue
                if proc == "inclusive":
                    d_list = all_list
                    out_f = os.path.join(args.out_base, era, "btag_eff.root")
                else:
                    d_list = proc_map.get(proc, [])
                    out_f = os.path.join(args.out_base, era, f"btag_eff_{proc}.root")
                if d_list:
                    jobs.append({
                        "era": era,
                        "process": proc,
                        "input_dirs": d_list,
                        "outfile": out_f,
                    })

    # Filter out existing outputs if requested
    pending_jobs = []
    for job in jobs:
        if args.skip_existing and os.path.isfile(job["outfile"]):
            print(f"[SKIP] Output already exists: {job['outfile']}")
            continue
        pending_jobs.append(job)

    print(f"\nTotal planned jobs: {len(pending_jobs)} (across {len(set(j['era'] for j in pending_jobs))} eras)")
    print("-" * 80)
    print(f"{'Era':<10} | {'Process':<12} | {'# Dirs':<8} | {'Target Output File'}")
    print("-" * 80)
    for j in pending_jobs:
        print(f"{j['era']:<10} | {j['process']:<12} | {len(j['input_dirs']):<8} | {j['outfile']}")
    print("-" * 80)

    if args.dry_run:
        print("\n[INFO] Dry-run completed. No commands executed.")
        return

    if not os.path.isdir(args.base_dir):
        print(f"\n[INFO] Base directory {args.base_dir} is not present on this machine.")
        print("       Run this command on the analysis server where /data2 is mounted:")
        print(f"       python3 {compute_script} ...")
        return

    # Execute jobs using ProcessPoolExecutor
    print(f"\nLaunching {len(pending_jobs)} jobs with {args.jobs} worker processes...\n")
    success_count = 0
    fail_count = 0

    with ProcessPoolExecutor(max_workers=args.jobs) as executor:
        future_to_job = {
            executor.submit(run_single_job, job, compute_script, extra_args): job
            for job in pending_jobs
        }
        for future in as_completed(future_to_job):
            job = future_to_job[future]
            try:
                success, _, stdout = future.result()
                if success:
                    success_count += 1
                else:
                    fail_count += 1
            except Exception as exc:
                print(f"[ERROR] Exception occurred for job {job['era']} | {job['process']}: {exc}")
                fail_count += 1

    print("\n" + "=" * 80)
    print(f"Batch Execution Summary: {success_count} succeeded, {fail_count} failed out of {len(pending_jobs)} jobs.")
    print("=" * 80)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compute_btag_efficiency.py

Extracts 2D b-tagging efficiency maps (h2_eff_b, h2_eff_c, h2_eff_light) as a function
of (pT, |eta|) for CMS BTV Method 1a (Event Reweighting) from Skimmed NanoAOD MC files.

Usage:
  python compute_btag_efficiency.py -Y 2022 -I /path/to/skimmed_mc.root -O data/BTV/2022_Summer22/btag_eff.root
  python compute_btag_efficiency.py -Y 2024 -D /path/to/skimmed_TTto2L2Nu_dir/ -O data/BTV/2024_Summer24/btag_eff.root
"""

import sys, os, glob, argparse, math
import array

# Working point thresholds (Medium WP) per era
# Matches s_btagcut_map in CommonTools/src/impl_Selection.cpp & analysis_config.json
BTAG_WP_M = {
    "2022": 0.2450,
    "2022_Summer22": 0.2450,
    "2022EE": 0.2605,
    "2022_Summer22EE": 0.2605,
    "2023": 0.1917,
    "2023_Summer23": 0.1917,
    "2023BPix": 0.1919,
    "2023_Summer23BPix": 0.1919,
    "2024": 0.1272,
    "2024_Summer24": 0.1272,
}

ETA_BINS = [0.0, 0.6, 1.2, 2.1, 2.5]

def get_tagger_col(year):
    if "2024" in year:
        return "Jet_btagUParTAK4B"
    return "Jet_btagPNetB"

def get_wp_threshold(year):
    for key, val in BTAG_WP_M.items():
        if key.lower() == year.lower() or key.lower() in year.lower():
            return val
    return 0.2450

def to_int(val, default=0, signed=False):
    """
    Safely convert PyROOT TTree values (int, UChar_t, Char_t, Bool_t, bytes, string) to int.
    Handles ROOT's Python 3 string mapping for char/uchar branches (e.g. '\x1e' -> 30, '\x05' -> 5).
    """
    if val is None:
        return default
    if isinstance(val, int):
        return val
    if isinstance(val, bool):
        return int(val)
    if isinstance(val, (bytes, bytearray)):
        b = val[0] if len(val) > 0 else default
        if signed and b >= 128:
            return b - 256
        return b
    if isinstance(val, str):
        if not val:
            return default
        if val.isdigit() or (val.startswith('-') and val[1:].isdigit()):
            return int(val)
        if len(val) == 1:
            code = ord(val)
            if signed and code >= 128:
                return code - 256
            return code
        try:
            return int(val)
        except Exception:
            code = ord(val[0])
            if signed and code >= 128:
                return code - 256
            return code
    try:
        return int(val)
    except Exception:
        return default

def to_float(val, default=0.0):
    """Safely convert value to float, handling string/char buffer artifacts."""
    if val is None:
        return default
    if isinstance(val, (float, int)):
        return float(val)
    if isinstance(val, str):
        if len(val) == 1:
            return float(ord(val))
        try:
            return float(val)
        except ValueError:
            return default
    try:
        return float(val)
    except Exception:
        return default

def safe_len(obj, fallback_count=None):
    """Safely determine collection or buffer length without throwing TypeError on C-buffers."""
    if fallback_count is not None:
        return to_int(fallback_count, 0)
    if hasattr(obj, "__len__"):
        try:
            return len(obj)
        except TypeError:
            pass
    if hasattr(obj, "size"):
        try:
            return int(obj.size())
        except Exception:
            pass
    return 0

def delta_r(eta1, phi1, eta2, phi2):
    deta = eta1 - eta2
    dphi = math.remainder(phi1 - phi2, 2.0 * math.pi)
    return math.hypot(deta, dphi)

def find_root_files(inputs):
    """
    Recursively find all .root files from a list of files, directories, or glob patterns.
    Handles nested subdirectories (e.g. QCD pt-binned samples in subfolders).
    """
    found = []
    if isinstance(inputs, str):
        inputs = [inputs]
    for inp in inputs:
        # Split comma-separated paths
        for token in inp.split(","):
            token = token.strip()
            if not token:
                continue
            matched = glob.glob(token)
            if not matched and os.path.exists(token):
                matched = [token]
            for item in matched:
                if os.path.isfile(item) and item.endswith(".root"):
                    found.append(item)
                elif os.path.isdir(item):
                    for root_dir, _, files in os.walk(item, followlinks=True):
                        for f in files:
                            if f.endswith(".root"):
                                found.append(os.path.join(root_dir, f))
    return sorted(list(set(found)))

def main():
    parser = argparse.ArgumentParser(description="Compute 2D MC b-tag efficiency maps for Method 1a aligned with analysis selection")
    parser.add_argument("-I", "--infile", dest="infile", nargs="*", default=[], help="Single/multiple input ROOT files, comma-separated files, or glob patterns")
    parser.add_argument("-D", "--indir", dest="indir", nargs="*", default=[], help="Directory or directories containing skimmed MC files (searched recursively)")
    parser.add_argument("-O", "--outfile", dest="outfile", type=str, default="", help="Output ROOT file path (default: data/BTV/<year>/btag_eff_<process>.root)")
    parser.add_argument("-Y", "--year", dest="year", type=str, default="", help="Era name: 2022, 2022EE, 2023, 2023BPix, 2024 (required unless --auto-skim is used)")
    parser.add_argument("-P", "--process", dest="process", type=str, default="ttbar", help="Process group aligned with plotIt (ttbar, singletop, wjets, dyjets, qcd, diboson, ttx, signal, inclusive) [default: ttbar]")
    parser.add_argument("--auto-skim", dest="auto_skim", type=str, nargs="?", const="/data2/common/skimmed_NanoAOD/skim_0812_LFV/muon/mc/", default="", help="Root directory of skimmed MC to automatically scan and compute all process maps across eras [default: /data2/common/skimmed_NanoAOD/skim_0812_LFV/muon/mc/]")
    parser.add_argument("-j", "--jobs", dest="jobs", type=int, default=4, help="Number of concurrent worker processes when using --auto-skim [default: 4]")
    parser.add_argument("--dry-run", dest="dry_run", action="store_true", help="Display discovered jobs and command plan without executing (for --auto-skim)")
    parser.add_argument("--skip-existing", dest="skip_existing", action="store_true", help="Skip jobs whose output ROOT file already exists")
    parser.add_argument("-C", "--channel", dest="channel", type=str, default="muon", choices=["muon", "electron", "both"], help="Lepton channel for selection (analysis default: muon)")
    parser.add_argument("--wp", dest="wp", type=float, default=-1.0, help="Override b-tag Medium WP threshold")
    parser.add_argument("--eta-bins", dest="eta_bins", type=str, default="single", help="Eta binning: 'single' ([0.0, max_eta], default), 'barrel-endcap' ([0.0, 1.4442, max_eta]), or comma-separated edges")
    parser.add_argument("--pt-bins", dest="pt_bins", type=str, default="", help="Custom comma-separated pT bin edges (default: standard BTV boundaries starting from min_pt)")
    parser.add_argument("--muon-pt", dest="muon_pt", type=float, default=50.0, help="Minimum signal muon pT (analysis default: 50.0 GeV)")
    parser.add_argument("--elec-pt", dest="elec_pt", type=float, default=50.0, help="Minimum signal electron pT (analysis default: 50.0 GeV)")
    parser.add_argument("--tau-pt", dest="tau_pt", type=float, default=40.0, help="Minimum signal tau pT (analysis default: 40.0 GeV)")
    parser.add_argument("--min-pt", dest="min_pt", type=float, default=40.0, help="Minimum jet pT (analysis default: 40.0 GeV)")
    parser.add_argument("--max-eta", dest="max_eta", type=float, default=2.5, help="Maximum jet |eta| (analysis default: 2.5)")
    parser.add_argument("--min-njets", dest="min_njets", type=int, default=3, help="Minimum clean jets per event (analysis default: 3)")
    parser.add_argument("--dr-lep", dest="dr_lep", type=float, default=0.4, help="DeltaR overlap removal with leptons/taus (analysis default: 0.4)")
    parser.add_argument("--max-events", dest="max_events", type=int, default=-1, help="Max events to process (default: all)")
    parser.add_argument("--full-baseline", dest="full_baseline", action="store_true", default=True, help="Apply full analysis baseline: N_lep==1, N_veto==0, N_clean_tau==1, OS charge, N_jets>=3 (default: True)")
    parser.add_argument("--no-full-baseline", dest="full_baseline", action="store_false", help="Disable lepton/tau event multiplicity cuts and only apply N_jets>=3")
    args = parser.parse_args()

    # If --auto-skim is requested, delegate execution to run_all_btag_efficiency
    if args.auto_skim:
        this_dir = os.path.dirname(os.path.abspath(__file__))
        if this_dir not in sys.path:
            sys.path.insert(0, this_dir)
        import run_all_btag_efficiency
        sys.argv = [
            "run_all_btag_efficiency.py",
            "-B", args.auto_skim,
            "-j", str(args.jobs),
            "--channel", args.channel,
            "--eta-bins", args.eta_bins,
        ]
        if args.year:
            sys.argv.extend(["-Y", args.year])
        if args.process and args.process != "ttbar":
            sys.argv.extend(["-P", args.process])
        if args.pt_bins:
            sys.argv.extend(["--pt-bins", args.pt_bins])
        if args.max_events > 0:
            sys.argv.extend(["--max-events", str(args.max_events)])
        if args.dry_run:
            sys.argv.append("--dry-run")
        if args.skip_existing:
            sys.argv.append("--skip-existing")
        run_all_btag_efficiency.main()
        return

    if not args.year:
        parser.error("-Y/--year is required unless --auto-skim is specified.")

    # Deferred ROOT import so script can be inspected without ROOT installed
    try:
        import ROOT
    except ImportError:
        print("ERROR: ROOT module not found. Please run this in an environment with ROOT (e.g. LCG view / cvmfs).")
        sys.exit(1)

    ROOT.TH1.SetDefaultSumw2(True)

    # Collect input files recursively
    input_items = []
    if args.infile:
        input_items.extend(args.infile)
    if args.indir:
        input_items.extend(args.indir)

    input_files = find_root_files(input_items)

    if not input_files:
        print(f"ERROR: No input ROOT files found in {input_items}! Specify -I <file/glob> or -D <dir(s)>.")
        sys.exit(1)

    # Determine output file path
    if not args.outfile:
        out_dir = os.path.join("data", "BTV", args.year)
        if args.process.lower() in ["inclusive", "default", "common"]:
            args.outfile = os.path.join(out_dir, "btag_eff.root")
        else:
            args.outfile = os.path.join(out_dir, f"btag_eff_{args.process.lower()}.root")
    elif os.path.isdir(args.outfile) or args.outfile.endswith("/"):
        out_dir = args.outfile
        if args.process.lower() in ["inclusive", "default", "common"]:
            args.outfile = os.path.join(out_dir, "btag_eff.root")
        else:
            args.outfile = os.path.join(out_dir, f"btag_eff_{args.process.lower()}.root")
    else:
        out_dir = os.path.dirname(args.outfile)

    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    print(f">> Era: {args.year}, Process: {args.process}, Channel: {args.channel}")
    print(f">> Total input files found (recursive): {len(input_files)}")
    print(f">> Output file: {args.outfile}")

    tagger_col = get_tagger_col(args.year)
    btag_cut = args.wp if args.wp > 0 else get_wp_threshold(args.year)
    print(f">> Tagger column: {tagger_col}, Medium WP threshold: {btag_cut}")
    print(f">> Jet Selection: pT > {args.min_pt} GeV, |eta| < {args.max_eta}, JetID & LepOverlap dR >= {args.dr_lep}")
    if args.full_baseline:
        print(f">> Event Selection: FULL ANALYSIS BASELINE (N_lep==1, N_veto_lep==0, N_clean_tau==1, leptau_charge < 0, N_clean_jets >= {args.min_njets}, NO b-tag cut)")
    else:
        print(f">> Event Selection: INCLUSIVE JET BASELINE (N_clean_jets >= {args.min_njets}, NO b-tag cut)")

    # Dynamic pT binning
    if args.pt_bins.strip():
        pt_edges = [float(x.strip()) for x in args.pt_bins.split(",") if x.strip()]
    else:
        standard_pt_edges = [20.0, 30.0, 40.0, 50.0, 70.0, 100.0, 140.0, 200.0, 300.0, 600.0, 1000.0]
        pt_edges = [b for b in standard_pt_edges if b >= args.min_pt]
        if not pt_edges or pt_edges[0] != args.min_pt:
            pt_edges.insert(0, args.min_pt)
        if pt_edges[-1] < 1000.0:
            pt_edges.append(1000.0)

    # Dynamic |eta| binning (default: single bin [0.0, max_eta])
    if args.eta_bins.lower() == "single":
        eta_edges = [0.0, args.max_eta]
    elif args.eta_bins.lower() in ["barrel-endcap", "be"]:
        eta_edges = [0.0, 1.4442, args.max_eta]
    elif args.eta_bins.lower() in ["four", "homogeneous"]:
        eta_edges = [0.0, 0.6, 1.2, 2.1, args.max_eta]
    else:
        eta_edges = [float(x.strip()) for x in args.eta_bins.split(",") if x.strip()]
        if eta_edges[0] != 0.0:
            eta_edges.insert(0, 0.0)
        if eta_edges[-1] < args.max_eta:
            eta_edges.append(args.max_eta)

    pt_bins = array.array('d', pt_edges)
    eta_bins = array.array('d', eta_edges)
    n_pt = len(pt_edges) - 1
    n_eta = len(eta_edges) - 1

    print(f">> pT bin edges: {pt_edges}")
    print(f">> |eta| bin edges: {eta_edges}")

    # Book 2D Denominator and Numerator histograms: X = pT, Y = |eta|
    flavors = ["b", "c", "light"]
    h2_denom = {}
    h2_num = {}
    for fl in flavors:
        h2_denom[fl] = ROOT.TH2F(f"h2_denom_{fl}", f"Total jets ({fl});Jet p_{{T}} [GeV];Jet |#eta|", n_pt, pt_bins, n_eta, eta_bins)
        h2_num[fl]   = ROOT.TH2F(f"h2_num_{fl}",   f"Tagged jets ({fl});Jet p_{{T}} [GeV];Jet |#eta|", n_pt, pt_bins, n_eta, eta_bins)

    chain = ROOT.TChain("Events")
    valid_file_count = 0
    for f in input_files:
        try:
            tf = ROOT.TFile.Open(f, "READ")
            if tf and not tf.IsZombie():
                tt = tf.Get("Events")
                if tt and tt.GetEntries() > 0:
                    chain.Add(f)
                    valid_file_count += 1
                tf.Close()
            elif tf:
                tf.Close()
        except Exception:
            pass

    if valid_file_count == 0:
        # Fallback: if TFile.Open failed due to remote paths or permission, try chain.Add directly
        for f in input_files:
            chain.Add(f)

    total_events = chain.GetEntries()
    print(f">> Total entries in TChain: {total_events}")
    if total_events == 0:
        print(f"ERROR: TChain has 0 entries across {len(input_files)} input files.")
        sys.exit(1)

    print(">> Processing events...")
    n_processed = 0
    n_passed_events = 0

    for ev in chain:
        n_processed += 1
        if args.max_events > 0 and n_processed > args.max_events:
            break
        if n_processed % 100000 == 0:
            print(f"   Processed {n_processed}/{total_events} events (passed baseline: {n_passed_events})...")

        w = 1.0
        if hasattr(ev, "genWeight"):
            w = 1.0 if to_float(getattr(ev, "genWeight", 1.0)) >= 0 else -1.0

        # 0. Primary vertex cut (using to_int to handle PyROOT UChar_t -> str/bytes mapping)
        pv_npvs = getattr(ev, "PV_npvsGood", None)
        if args.full_baseline and pv_npvs is not None and to_int(pv_npvs) <= 0:
            continue

        # Fast path 1: Check if skimmed branches exist (nmuonpass, nelepass, nvetomuons, nvetoelepass)
        has_skim_leptons = hasattr(ev, "nmuonpass") and hasattr(ev, "nvetomuons") and hasattr(ev, "nvetoelepass")

        if has_skim_leptons and args.full_baseline:
            # Check single lepton requirement
            if args.channel == "muon":
                if to_int(getattr(ev, "nmuonpass", 0)) != 1:
                    continue
            elif args.channel == "electron":
                if to_int(getattr(ev, "nelepass", 0)) != 1:
                    continue
            else:
                if to_int(getattr(ev, "nmuonpass", 0)) != 1 and to_int(getattr(ev, "nelepass", 0)) != 1:
                    continue

            # Check veto lepton requirement (no additional loose leptons)
            if to_int(getattr(ev, "nvetomuons", 0)) != 0 or to_int(getattr(ev, "nvetoelepass", 0)) != 0:
                continue

        # Fast path 2: Early exit if raw tau or jet count is insufficient to pass baseline
        n_tau_raw = to_int(getattr(ev, "nTau", None))
        n_jet_raw = to_int(getattr(ev, "nJet", None))
        if args.full_baseline:
            if n_tau_raw < 1 or n_jet_raw < args.min_njets:
                continue

        signal_leptons = []  # list of (eta, phi, charge)
        clean_taus = []      # list of (eta, phi, charge)

        if has_skim_leptons:
            # Retrieve signal lepton eta/phi/charge
            if args.channel in ["muon", "both"] and to_int(getattr(ev, "nmuonpass", 0)) == 1:
                mu_etas = getattr(ev, "Muon_eta", None)
                mu_phis = getattr(ev, "Muon_phi", None)
                mu_qs = getattr(ev, "Muon_charge", None)
                n_mu = safe_len(mu_etas, getattr(ev, "nMuon", None))
                if n_mu > 0 and mu_etas is not None and mu_phis is not None:
                    q = to_int(mu_qs[0], 0, signed=True) if mu_qs is not None else 0
                    signal_leptons.append((to_float(mu_etas[0]), to_float(mu_phis[0]), q))
            elif to_int(getattr(ev, "nelepass", 0)) == 1:
                el_etas = getattr(ev, "Electron_eta", None)
                el_phis = getattr(ev, "Electron_phi", None)
                el_qs = getattr(ev, "Electron_charge", None)
                n_el = safe_len(el_etas, getattr(ev, "nElectron", None))
                if n_el > 0 and el_etas is not None and el_phis is not None:
                    q = to_int(el_qs[0], 0, signed=True) if el_qs is not None else 0
                    signal_leptons.append((to_float(el_etas[0]), to_float(el_phis[0]), q))
        else:
            # Fallback path: Full manual object inspection
            n_veto_muons = 0
            n_veto_elecs = 0

            mu_pts = getattr(ev, "Muon_pt", None)
            mu_etas = getattr(ev, "Muon_eta", None)
            mu_phis = getattr(ev, "Muon_phi", None)
            mu_charges = getattr(ev, "Muon_charge", None)
            mu_tight = getattr(ev, "Muon_tightId", None)
            mu_loose = getattr(ev, "Muon_looseId", None)
            mu_iso = getattr(ev, "Muon_pfRelIso04_all", None)

            n_mu_cnt = safe_len(mu_pts, getattr(ev, "nMuon", None))
            sig_mu_indices = set()
            if mu_pts is not None and mu_etas is not None and mu_phis is not None:
                for m_i in range(n_mu_cnt):
                    m_pt = to_float(mu_pts[m_i])
                    m_eta = to_float(mu_etas[m_i])
                    m_phi = to_float(mu_phis[m_i])
                    m_q = to_int(mu_charges[m_i], 0, signed=True) if mu_charges is not None else 0
                    is_tight = bool(to_int(mu_tight[m_i])) if mu_tight is not None else True
                    iso_val = to_float(mu_iso[m_i]) if mu_iso is not None else 0.0

                    if m_pt > args.muon_pt and abs(m_eta) < 2.4 and is_tight and iso_val < 0.15:
                        if args.channel in ["muon", "both"]:
                            signal_leptons.append((m_eta, m_phi, m_q))
                            sig_mu_indices.add(m_i)

                for m_i in range(n_mu_cnt):
                    if m_i in sig_mu_indices:
                        continue
                    m_pt = to_float(mu_pts[m_i])
                    m_eta = to_float(mu_etas[m_i])
                    is_loose = bool(to_int(mu_loose[m_i])) if mu_loose is not None else False
                    iso_val = to_float(mu_iso[m_i]) if mu_iso is not None else 1.0
                    if m_pt > 15.0 and abs(m_eta) < 2.4 and is_loose and iso_val < 0.25:
                        n_veto_muons += 1

            el_pts = getattr(ev, "Electron_pt", None)
            el_etas = getattr(ev, "Electron_eta", None)
            el_phis = getattr(ev, "Electron_phi", None)
            el_charges = getattr(ev, "Electron_charge", None)
            el_mva = getattr(ev, "Electron_mvaIso_WP90", None)
            el_cutbased = getattr(ev, "Electron_cutBased", None)

            n_el_cnt = safe_len(el_pts, getattr(ev, "nElectron", None))
            sig_el_indices = set()
            if el_pts is not None and el_etas is not None and el_phis is not None:
                for e_i in range(n_el_cnt):
                    e_pt = to_float(el_pts[e_i])
                    e_eta = to_float(el_etas[e_i])
                    e_abseta = abs(e_eta)
                    e_phi = to_float(el_phis[e_i])
                    e_q = to_int(el_charges[e_i], 0, signed=True) if el_charges is not None else 0
                    is_mva90 = bool(to_int(el_mva[e_i])) if el_mva is not None else True

                    if e_pt > args.elec_pt and e_abseta < 2.5 and not (1.4442 < e_abseta < 1.566) and is_mva90:
                        if args.channel in ["electron", "both"]:
                            signal_leptons.append((e_eta, e_phi, e_q))
                            sig_el_indices.add(e_i)

                for e_i in range(n_el_cnt):
                    if e_i in sig_el_indices:
                        continue
                    e_pt = to_float(el_pts[e_i])
                    e_abseta = abs(to_float(el_etas[e_i]))
                    cb_val = to_int(el_cutbased[e_i]) if el_cutbased is not None else 0
                    if e_pt > 15.0 and e_abseta < 2.5 and not (1.4442 < e_abseta < 1.566) and cb_val >= 1:
                        n_veto_elecs += 1

            if args.full_baseline:
                if len(signal_leptons) != 1:
                    continue
                if n_veto_muons > 0 or n_veto_elecs > 0:
                    continue

        if args.full_baseline and len(signal_leptons) != 1:
            continue

        # Clean Tau selection
        tau_pts = getattr(ev, "Tau_pt", None)
        tau_etas = getattr(ev, "Tau_eta", None)
        tau_phis = getattr(ev, "Tau_phi", None)
        tau_charges = getattr(ev, "Tau_charge", None)
        tau_dm_new = getattr(ev, "Tau_idDecayModeNewDMs", None)
        tau_dm = getattr(ev, "Tau_decayMode", None)
        tau_vsjet = getattr(ev, "Tau_idDeepTau2018v2p5VSjet", None)
        tau_vsmu = getattr(ev, "Tau_idDeepTau2018v2p5VSmu", None)
        tau_vse = getattr(ev, "Tau_idDeepTau2018v2p5VSe", None)

        vsjet_th = 7
        vsmu_th = 4 if args.channel == "muon" else 1
        vse_th = 2 if args.channel == "muon" else 6

        n_tau_cnt = safe_len(tau_pts, n_tau_raw)
        if tau_pts is not None and tau_etas is not None and tau_phis is not None:
            for t_i in range(n_tau_cnt):
                t_pt = to_float(tau_pts[t_i])
                t_eta = to_float(tau_etas[t_i])
                t_phi = to_float(tau_phis[t_i])
                t_q = to_int(tau_charges[t_i], 0, signed=True) if tau_charges is not None else 0

                if t_pt <= args.tau_pt or abs(t_eta) >= 2.5:
                    continue
                if tau_dm_new is not None and to_int(tau_dm_new[t_i]) == 0:
                    continue
                if tau_dm is not None:
                    dm_val = to_int(tau_dm[t_i])
                    if dm_val == 5 or dm_val == 6:
                        continue
                if tau_vsjet is not None and to_int(tau_vsjet[t_i]) < vsjet_th:
                    continue
                if tau_vsmu is not None and to_int(tau_vsmu[t_i]) < vsmu_th:
                    continue
                if tau_vse is not None and to_int(tau_vse[t_i]) < vse_th:
                    continue

                # Overlap removal with signal leptons
                overlap_lep = False
                for l_eta, l_phi, _ in signal_leptons:
                    if delta_r(t_eta, t_phi, l_eta, l_phi) < args.dr_lep:
                        overlap_lep = True
                        break
                if overlap_lep:
                    continue

                clean_taus.append((t_eta, t_phi, t_q))

        if args.full_baseline:
            if len(clean_taus) != 1:
                continue
            # Opposite-sign charge cut: leptau_charge < 0
            lep_charge = signal_leptons[0][2]
            tau_charge = clean_taus[0][2]
            if lep_charge * tau_charge >= 0:
                continue

        # Clean Jet selection and overlap removal (lepjetoverlap && taujetoverlap)
        pts = getattr(ev, "Jet_pt", None)
        etas = getattr(ev, "Jet_eta", None)
        phis = getattr(ev, "Jet_phi", None)
        hadflav = getattr(ev, "Jet_hadronFlavour", None)
        scores = getattr(ev, tagger_col, None)

        if pts is None or etas is None or hadflav is None or scores is None:
            continue

        n_jet_cnt = safe_len(pts, n_jet_raw)
        if n_jet_cnt < args.min_njets:
            continue

        pass_jet_id = getattr(ev, "Jet_passJetIdTightLepVeto", None)
        jet_id = getattr(ev, "Jet_jetId", None)
        mu_ef = getattr(ev, "Jet_muEF", None)
        ch_em_ef = getattr(ev, "Jet_chEmEF", None)

        clean_jets = []

        for j in range(n_jet_cnt):
            pt = to_float(pts[j])
            eta = to_float(etas[j])
            abseta = abs(eta)
            phi = to_float(phis[j]) if phis is not None else 0.0

            # Kinematic cut
            if pt <= args.min_pt or abseta >= args.max_eta:
                continue

            # Jet ID & energy fraction cuts matching analysis_config.json
            if pass_jet_id is not None and to_float(pass_jet_id[j]) < 1.0:
                continue
            if pass_jet_id is None and jet_id is not None and to_int(jet_id[j]) < 2:
                continue
            if mu_ef is not None and to_float(mu_ef[j]) >= 0.8:
                continue
            if ch_em_ef is not None and to_float(ch_em_ef[j]) >= 0.8:
                continue

            # Overlap removal with signal leptons (lepjetoverlap: dR >= dr_lep)
            is_overlap = False
            if phis is not None:
                for l_eta, l_phi, _ in signal_leptons:
                    if delta_r(eta, phi, l_eta, l_phi) < args.dr_lep:
                        is_overlap = True
                        break
                if is_overlap:
                    continue

                # Overlap removal with clean taus (taujetoverlap: dR >= dr_lep)
                for t_eta, t_phi, _ in clean_taus:
                    if delta_r(eta, phi, t_eta, t_phi) < args.dr_lep:
                        is_overlap = True
                        break
                if is_overlap:
                    continue

            flav = to_int(hadflav[j])
            score = to_float(scores[j])
            clean_jets.append((pt, abseta, flav, score))

        # Event-level multiplicity cut: N_clean_jets >= min_njets (default: 3)
        # Note: BTV principle: NEVER apply b-tag cut here!
        if len(clean_jets) < args.min_njets:
            continue

        n_passed_events += 1

        # Fill denominator and numerator for clean jets
        for pt, abseta, flav, score in clean_jets:
            if flav == 5:
                fl_key = "b"
            elif flav == 4:
                fl_key = "c"
            else:
                fl_key = "light"

            h2_denom[fl_key].Fill(pt, abseta, w)
            if score >= btag_cut:
                h2_num[fl_key].Fill(pt, abseta, w)

    # Compute efficiencies and fill fallbacks for empty bins
    out_dir = os.path.dirname(os.path.abspath(args.outfile))
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    fout = ROOT.TFile(args.outfile, "RECREATE")
    default_effs = {"b": 0.65, "c": 0.15, "light": 0.01}

    for fl in flavors:
        h2_denom[fl].Write()
        h2_num[fl].Write()

        h2_eff = h2_num[fl].Clone(f"h2_eff_{fl}")
        h2_eff.SetTitle(f"B-tagging Efficiency ({fl});Jet p_{{T}} [GeV];Jet |#eta|")
        h2_eff.Divide(h2_num[fl], h2_denom[fl], 1.0, 1.0, "B")

        # Sanity check & empty bin smoothing with default fallback
        for bx in range(1, n_pt + 1):
            for by in range(1, n_eta + 1):
                eff = h2_eff.GetBinContent(bx, by)
                denom_val = h2_denom[fl].GetBinContent(bx, by)
                if denom_val <= 0 or eff <= 0.0 or eff > 1.0:
                    h2_eff.SetBinContent(bx, by, default_effs[fl])
                    h2_eff.SetBinError(bx, by, 0.0)

        h2_eff.Write()
        print(f">> Generated h2_eff_{fl} (mean eff: {h2_eff.Integral() / (n_pt * n_eta):.4f})")

    fout.Close()
    print(f">> Successfully saved efficiency maps to: {args.outfile}")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        import traceback
        print(f"\n[FATAL ERROR in compute_btag_efficiency]: {e}", file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)

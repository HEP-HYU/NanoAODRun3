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

def delta_r(eta1, phi1, eta2, phi2):
    deta = eta1 - eta2
    dphi = math.remainder(phi1 - phi2, 2.0 * math.pi)
    return math.hypot(deta, dphi)

def main():
    parser = argparse.ArgumentParser(description="Compute 2D MC b-tag efficiency maps for Method 1a aligned with analysis selection")
    parser.add_argument("-I", "--infile", dest="infile", type=str, default="", help="Single input ROOT file or comma-separated files")
    parser.add_argument("-D", "--indir", dest="indir", type=str, default="", help="Directory containing skimmed MC ROOT files")
    parser.add_argument("-O", "--outfile", dest="outfile", type=str, required=True, help="Output ROOT file path (e.g. data/BTV/2022_Summer22/btag_eff.root)")
    parser.add_argument("-Y", "--year", dest="year", type=str, required=True, help="Era name: 2022, 2022EE, 2023, 2023BPix, 2024")
    parser.add_argument("--channel", dest="channel", type=str, default="muon", choices=["muon", "electron", "both"], help="Lepton channel for selection (analysis default: muon)")
    parser.add_argument("--wp", dest="wp", type=float, default=-1.0, help="Override b-tag Medium WP threshold")
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

    # Deferred ROOT import so script can be inspected without ROOT installed
    try:
        import ROOT
    except ImportError:
        print("ERROR: ROOT module not found. Please run this in an environment with ROOT (e.g. LCG view / cvmfs).")
        sys.exit(1)

    ROOT.TH1.SetDefaultSumw2(True)

    input_files = []
    if args.infile:
        for f in args.infile.split(","):
            input_files.extend(glob.glob(f.strip()))
    if args.indir:
        input_files.extend(glob.glob(os.path.join(args.indir, "*.root")))

    if not input_files:
        print("ERROR: No input files found! Specify -I <file> or -D <dir>.")
        sys.exit(1)

    print(f">> Total input files: {len(input_files)}")
    print(f">> Era: {args.year}, Channel: {args.channel}")

    tagger_col = get_tagger_col(args.year)
    btag_cut = args.wp if args.wp > 0 else get_wp_threshold(args.year)
    print(f">> Tagger column: {tagger_col}, Medium WP threshold: {btag_cut}")
    print(f">> Jet Selection: pT > {args.min_pt} GeV, |eta| < {args.max_eta}, JetID & LepOverlap dR >= {args.dr_lep}")
    if args.full_baseline:
        print(f">> Event Selection: FULL ANALYSIS BASELINE (N_lep==1, N_veto_lep==0, N_clean_tau==1, leptau_charge < 0, N_clean_jets >= {args.min_njets}, NO b-tag cut)")
    else:
        print(f">> Event Selection: INCLUSIVE JET BASELINE (N_clean_jets >= {args.min_njets}, NO b-tag cut)")

    # Dynamic pT binning matching min_pt and standard BTV boundaries
    standard_pt_edges = [20.0, 30.0, 40.0, 50.0, 70.0, 100.0, 140.0, 200.0, 300.0, 600.0, 1000.0]
    pt_edges = [b for b in standard_pt_edges if b >= args.min_pt]
    if not pt_edges or pt_edges[0] != args.min_pt:
        pt_edges.insert(0, args.min_pt)
    if pt_edges[-1] < 1000.0:
        pt_edges.append(1000.0)

    pt_bins = array.array('d', pt_edges)
    eta_bins = array.array('d', ETA_BINS)
    n_pt = len(pt_edges) - 1
    n_eta = len(ETA_BINS) - 1

    print(f">> pT bin edges: {pt_edges}")
    print(f">> |eta| bin edges: {list(ETA_BINS)}")

    # Book 2D Denominator and Numerator histograms: X = pT, Y = |eta|
    flavors = ["b", "c", "light"]
    h2_denom = {}
    h2_num = {}
    for fl in flavors:
        h2_denom[fl] = ROOT.TH2F(f"h2_denom_{fl}", f"Total jets ({fl});Jet p_{{T}} [GeV];Jet |#eta|", n_pt, pt_bins, n_eta, eta_bins)
        h2_num[fl]   = ROOT.TH2F(f"h2_num_{fl}",   f"Tagged jets ({fl});Jet p_{{T}} [GeV];Jet |#eta|", n_pt, pt_bins, n_eta, eta_bins)

    chain = ROOT.TChain("Events")
    for f in input_files:
        chain.Add(f)

    total_events = chain.GetEntries()
    print(f">> Total entries in TChain: {total_events}")
    if total_events == 0:
        print("ERROR: TChain has 0 entries.")
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
            w = 1.0 if ev.genWeight >= 0 else -1.0

        # 0. Primary vertex cut
        pv_npvs = getattr(ev, "PV_npvsGood", None)
        if args.full_baseline and pv_npvs is not None and pv_npvs <= 0:
            continue

        # 1. Signal & Veto Lepton Selection
        # Analysis baseline:
        #   Muon channel: nmuonpass == 1 && nvetoelepass == 0 && nvetomuons == 0
        #   Electron channel: nelepass == 1 && nvetoelepass == 0 && nvetomuons == 0
        signal_leptons = []  # list of (eta, phi, charge)
        n_veto_muons = 0
        n_veto_elecs = 0

        # Muon inspection
        mu_pts = getattr(ev, "Muon_pt", None)
        mu_etas = getattr(ev, "Muon_eta", None)
        mu_phis = getattr(ev, "Muon_phi", None)
        mu_charges = getattr(ev, "Muon_charge", None)
        mu_tight = getattr(ev, "Muon_tightId", None)
        mu_loose = getattr(ev, "Muon_looseId", None)
        mu_iso = getattr(ev, "Muon_pfRelIso04_all", None)

        sig_mu_indices = set()
        if mu_pts is not None and mu_etas is not None and mu_phis is not None:
            for m_i in range(len(mu_pts)):
                m_pt = mu_pts[m_i]
                m_eta = mu_etas[m_i]
                m_phi = mu_phis[m_i]
                m_q = mu_charges[m_i] if mu_charges is not None else 0
                is_tight = bool(mu_tight[m_i]) if mu_tight is not None else True
                iso_val = mu_iso[m_i] if mu_iso is not None else 0.0

                # Signal muon cut: pt > 50, |eta| < 2.4, tightId, iso < 0.15
                if m_pt > args.muon_pt and abs(m_eta) < 2.4 and is_tight and iso_val < 0.15:
                    if args.channel in ["muon", "both"]:
                        signal_leptons.append((m_eta, m_phi, m_q))
                        sig_mu_indices.add(m_i)

            # Veto muon check: !signal && pt > 15, |eta| < 2.4, looseId, iso < 0.25
            for m_i in range(len(mu_pts)):
                if m_i in sig_mu_indices:
                    continue
                m_pt = mu_pts[m_i]
                m_eta = mu_etas[m_i]
                is_loose = bool(mu_loose[m_i]) if mu_loose is not None else False
                iso_val = mu_iso[m_i] if mu_iso is not None else 1.0
                if m_pt > 15.0 and abs(m_eta) < 2.4 and is_loose and iso_val < 0.25:
                    n_veto_muons += 1

        # Electron inspection
        el_pts = getattr(ev, "Electron_pt", None)
        el_etas = getattr(ev, "Electron_eta", None)
        el_phis = getattr(ev, "Electron_phi", None)
        el_charges = getattr(ev, "Electron_charge", None)
        el_mva = getattr(ev, "Electron_mvaIso_WP90", None)
        el_cutbased = getattr(ev, "Electron_cutBased", None)

        sig_el_indices = set()
        if el_pts is not None and el_etas is not None and el_phis is not None:
            for e_i in range(len(el_pts)):
                e_pt = el_pts[e_i]
                e_eta = el_etas[e_i]
                e_abseta = abs(e_eta)
                e_phi = el_phis[e_i]
                e_q = el_charges[e_i] if el_charges is not None else 0
                is_mva90 = bool(el_mva[e_i]) if el_mva is not None else True

                # Signal electron cut: pt > 50, |eta| < 2.5, !(1.4442 < |eta| < 1.566), mvaIso_WP90
                if e_pt > args.elec_pt and e_abseta < 2.5 and not (1.4442 < e_abseta < 1.566) and is_mva90:
                    if args.channel in ["electron", "both"]:
                        signal_leptons.append((e_eta, e_phi, e_q))
                        sig_el_indices.add(e_i)

            # Veto electron check: !signal && pt > 15, |eta| < 2.5, !(1.4442 < |eta| < 1.566), cutBased >= 1
            for e_i in range(len(el_pts)):
                if e_i in sig_el_indices:
                    continue
                e_pt = el_pts[e_i]
                e_abseta = abs(el_etas[e_i])
                cb_val = el_cutbased[e_i] if el_cutbased is not None else 0
                if e_pt > 15.0 and e_abseta < 2.5 and not (1.4442 < e_abseta < 1.566) and cb_val >= 1:
                    n_veto_elecs += 1

        # Apply lepton baseline cuts if full_baseline enabled
        if args.full_baseline:
            if len(signal_leptons) != 1:
                continue
            if n_veto_muons > 0 or n_veto_elecs > 0:
                continue

        # 2. Signal Clean Tau selection (matching cleantau4vecs in analysis_config.json / impl_Selection.cpp)
        # Tau: pt > 40.0, |eta| < 2.5, idDecayModeNewDMs, decayMode != 5 && decayMode != 6
        # DeepTau WPs:
        #   Muon ch: VSjet >= 7 (VVTight), VSmu >= 4 (Tight), VSe >= 2 (VVLoose)
        #   Elec ch: VSjet >= 7 (VVTight), VSmu >= 1 (VLoose), VSe >= 6 (VTight)
        # leptauoverlap: dR(tau, signal_lepton) >= 0.4
        clean_taus = []  # list of (eta, phi, charge)
        tau_pts = getattr(ev, "Tau_pt", None)
        tau_etas = getattr(ev, "Tau_eta", None)
        tau_phis = getattr(ev, "Tau_phi", None)
        tau_charges = getattr(ev, "Tau_charge", None)
        tau_dm_new = getattr(ev, "Tau_idDecayModeNewDMs", None)
        tau_dm = getattr(ev, "Tau_decayMode", None)
        tau_vsjet = getattr(ev, "Tau_idDeepTau2018v2p5VSjet", None)
        tau_vsmu = getattr(ev, "Tau_idDeepTau2018v2p5VSmu", None)
        tau_vse = getattr(ev, "Tau_idDeepTau2018v2p5VSe", None)

        vsjet_th = 7  # VVTight
        vsmu_th = 4 if args.channel == "muon" else 1  # Tight (muon ch) or VLoose (elec ch)
        vse_th = 2 if args.channel == "muon" else 6   # VVLoose (muon ch) or VTight (elec ch)

        if tau_pts is not None and tau_etas is not None and tau_phis is not None:
            for t_i in range(len(tau_pts)):
                t_pt = tau_pts[t_i]
                t_eta = tau_etas[t_i]
                t_phi = tau_phis[t_i]
                t_q = tau_charges[t_i] if tau_charges is not None else 0
                if t_pt <= args.tau_pt or abs(t_eta) >= 2.5:
                    continue
                if tau_dm_new is not None and not tau_dm_new[t_i]:
                    continue
                if tau_dm is not None and (tau_dm[t_i] == 5 or tau_dm[t_i] == 6):
                    continue
                if tau_vsjet is not None and tau_vsjet[t_i] < vsjet_th:
                    continue
                if tau_vsmu is not None and tau_vsmu[t_i] < vsmu_th:
                    continue
                if tau_vse is not None and tau_vse[t_i] < vse_th:
                    continue

                # Overlap with signal leptons (leptauoverlap: dR >= 0.4)
                overlap_lep = False
                for l_eta, l_phi, _ in signal_leptons:
                    if delta_r(t_eta, t_phi, l_eta, l_phi) < args.dr_lep:
                        overlap_lep = True
                        break
                if overlap_lep:
                    continue

                clean_taus.append((t_eta, t_phi, t_q))

        # Apply tau baseline cuts if full_baseline enabled
        if args.full_baseline:
            if len(clean_taus) != 1:
                continue
            # Opposite-sign charge cut: leptau_charge < 0
            lep_charge = signal_leptons[0][2]
            tau_charge = clean_taus[0][2]
            if lep_charge * tau_charge >= 0:
                continue

        # 3. Clean Jet selection and overlap removal (lepjetoverlap && taujetoverlap)
        # Jet: pt > 40.0, |eta| < 2.5, passJetIdTightLepVeto == 1, muEF < 0.8, chEmEF < 0.8
        # Overlap: dR(jet, signal_lepton) >= 0.4 AND dR(jet, clean_tau) >= 0.4
        pts = getattr(ev, "Jet_pt", None)
        etas = getattr(ev, "Jet_eta", None)
        phis = getattr(ev, "Jet_phi", None)
        hadflav = getattr(ev, "Jet_hadronFlavour", None)
        scores = getattr(ev, tagger_col, None)

        if pts is None or etas is None or hadflav is None or scores is None:
            continue

        pass_jet_id = getattr(ev, "Jet_passJetIdTightLepVeto", None)
        jet_id = getattr(ev, "Jet_jetId", None)
        mu_ef = getattr(ev, "Jet_muEF", None)
        ch_em_ef = getattr(ev, "Jet_chEmEF", None)

        n_jets = min(len(pts), len(etas), len(hadflav), len(scores))
        clean_jets = []

        for j in range(n_jets):
            pt = pts[j]
            eta = etas[j]
            abseta = abs(eta)
            phi = phis[j] if phis is not None else 0.0

            # Kinematic cut
            if pt <= args.min_pt or abseta >= args.max_eta:
                continue

            # Jet ID & energy fraction cuts matching analysis_config.json
            if pass_jet_id is not None and pass_jet_id[j] < 1.0:
                continue
            if pass_jet_id is None and jet_id is not None and jet_id[j] < 2:
                continue
            if mu_ef is not None and mu_ef[j] >= 0.8:
                continue
            if ch_em_ef is not None and ch_em_ef[j] >= 0.8:
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

            clean_jets.append((pt, abseta, hadflav[j], scores[j]))

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
    main()

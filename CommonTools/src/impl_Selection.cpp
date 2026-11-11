#include "NanoAODAnalyzerrdframe.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
#include <random>
#include <chrono>
#include <ctime>
#include <typeinfo>
#include <regex>
#include "utility.h"
#include "correction.h"
#include "GEScaleSyst.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TH2F.h"
#include "TSystem.h"
#include "ROOT/RDFHelpers.hxx"
#include "Math/GenVector/VectorUtil.h"
using namespace std;

void NanoAODAnalyzerrdframe::selectElectrons(string cut, string vetocut) {
    if (cut.find("onlyveto") != std::string::npos){
        cout<<"selectElectrons: only vetoelecuts"<<endl;
        _rlm = _rlm.Define("vetoelecuts", vetocut)
                   .Define("nvetoelepass","Sum(vetoelecuts)");
    } else{
        cout<<"selectElectrons: electroncuts, vetoelecuts"<<endl;
        _rlm = _rlm.Define("elecuts", cut)               
                   .Define("vetoelecuts", vetocut);

        _rlm = _rlm.Define("nvetoelepass","Sum(vetoelecuts)")
                   .Redefine("Electron_pt", "Electron_pt[elecuts]")
                   .Redefine("Electron_eta", "Electron_eta[elecuts]")
                   .Redefine("Electron_phi", "Electron_phi[elecuts]")
                   .Redefine("Electron_mass", "Electron_mass[elecuts]")
                   .Redefine("Electron_charge", "Electron_charge[elecuts]")
                   //.Define("Electron_Id", "Electron_looseId[elecuts]")
                   .Redefine("Electron_pfRelIso03_all", "Electron_pfRelIso03_all[elecuts]")
                   .Define("Sel_eleidx", ::good_idx, {"elecuts"})
                   .Define("nelepass", "int(Electron_pt.size())")
                   .Define("lep4vecs", ::gen4vec, {"Electron_pt", "Electron_eta", "Electron_phi", "Electron_mass"});
        if (isDefined("Electron_deltaEtaSC")) {
            _rlm = _rlm.Redefine("Electron_deltaEtaSC", "Electron_deltaEtaSC[elecuts]");
        }
        if (isDefined("Electron_r9")) {
            _rlm = _rlm.Redefine("Electron_r9", "Electron_r9[elecuts]");
        }
        if (isDefined("Electron_seedGain")) {
            _rlm = _rlm.Redefine("Electron_seedGain", "Electron_seedGain[elecuts]");
        }
        if (isDefined("Electron_pt_uncor")) {
            auto elescaremet = [](const floats &pts_uncor, const floats &pts_corr, const floats &phis, float met, float metphi)->float {
                double metx = met * std::cos(metphi);
                double mety = met * std::sin(metphi);
                for (size_t i = 0; i < pts_uncor.size(); i++) {
                    double dpx = (pts_corr[i] - pts_uncor[i]) * std::cos(phis[i]);
                    double dpy = (pts_corr[i] - pts_uncor[i]) * std::sin(phis[i]);
                    metx -= dpx;
                    mety -= dpy;
                }
                return float(std::sqrt(metx * metx + mety * mety));
            };
            auto elescaremetphi = [](const floats &pts_uncor, const floats &pts_corr, const floats &phis, float met, float metphi)->float {
                double metx = met * std::cos(metphi);
                double mety = met * std::sin(metphi);
                for (size_t i = 0; i < pts_uncor.size(); i++) {
                    double dpx = (pts_corr[i] - pts_uncor[i]) * std::cos(phis[i]);
                    double dpy = (pts_corr[i] - pts_uncor[i]) * std::sin(phis[i]);
                    metx -= dpx;
                    mety -= dpy;
                }
                return float(std::atan2(mety, metx));
            };
            _rlm = _rlm.Redefine("Electron_pt_uncor", "Electron_pt_uncor[elecuts]")
                       .Define("skim_PuppiMET_pt_ele_unsmear", "PuppiMET_pt")
                       .Define("skim_PuppiMET_phi_ele_unsmear", "PuppiMET_phi")
                       .Redefine("PuppiMET_phi", elescaremetphi, {"Electron_pt_uncor", "Electron_pt", "Electron_phi", "skim_PuppiMET_pt_ele_unsmear", "skim_PuppiMET_phi_ele_unsmear"})
                       .Redefine("PuppiMET_pt",  elescaremet,    {"Electron_pt_uncor", "Electron_pt", "Electron_phi", "skim_PuppiMET_pt_ele_unsmear", "skim_PuppiMET_phi_ele_unsmear"});
        }
    }
}

void NanoAODAnalyzerrdframe::selectMuons(string cut, string vetocut) {
    if (cut.find("onlyveto") != std::string::npos){
        cout<<"selectMuons: only vetomuoncuts"<<endl;
        _rlm = _rlm.Define("vetomuoncuts", vetocut)
                   .Define("nvetomuons", "Sum(vetomuoncuts)");
    } else{
        cout<<"selectMuons: muoncuts, vetomuoncuts"<<endl;
        _rlm = _rlm.Define("muoncuts", cut);
        _rlm = _rlm.Define("vetomuoncuts", vetocut);

        _rlm = _rlm.Define("nvetomuons","Sum(vetomuoncuts)")
                   .Redefine("Muon_pt", "Muon_pt[muoncuts]")
                   .Redefine("Muon_eta", "Muon_eta[muoncuts]")
                   .Redefine("Muon_phi", "Muon_phi[muoncuts]")
                   .Redefine("Muon_mass", "Muon_mass[muoncuts]")
                   .Redefine("Muon_charge", "Muon_charge[muoncuts]")
                   .Redefine("Muon_looseId", "Muon_looseId[muoncuts]")
                   .Redefine("Muon_pfRelIso04_all", "Muon_pfRelIso04_all[muoncuts]")
                   .Define("Sel_muonidx", ::good_idx, {"muoncuts"})
                   .Define("nmuonpass", "int(Muon_pt.size())")
                   .Define("lep4vecs", ::gen4vec, {"Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass"});
        if (isDefined("Muon_tightId")) {
            _rlm = _rlm.Redefine("Muon_tightId", "Muon_tightId[muoncuts]");
        }
        if (isDefined("Muon_nTrackerLayers")) {
            _rlm = _rlm.Redefine("Muon_nTrackerLayers", "Muon_nTrackerLayers[muoncuts]");
        }
        if (isDefined("Muon_pt_uncor")) {
            auto muonscaremet = [](const floats &pts_uncor, const floats &pts_corr, const floats &phis, float met, float metphi)->float {
                double metx = met * std::cos(metphi);
                double mety = met * std::sin(metphi);
                for (size_t i = 0; i < pts_uncor.size(); i++) {
                    double dpx = (pts_corr[i] - pts_uncor[i]) * std::cos(phis[i]);
                    double dpy = (pts_corr[i] - pts_uncor[i]) * std::sin(phis[i]);
                    metx -= dpx;
                    mety -= dpy;
                }
                return float(std::sqrt(metx * metx + mety * mety));
            };
            auto muonscaremetphi = [](const floats &pts_uncor, const floats &pts_corr, const floats &phis, float met, float metphi)->float {
                double metx = met * std::cos(metphi);
                double mety = met * std::sin(metphi);
                for (size_t i = 0; i < pts_uncor.size(); i++) {
                    double dpx = (pts_corr[i] - pts_uncor[i]) * std::cos(phis[i]);
                    double dpy = (pts_corr[i] - pts_uncor[i]) * std::sin(phis[i]);
                    metx -= dpx;
                    mety -= dpy;
                }
                return float(std::atan2(mety, metx));
            };
            _rlm = _rlm.Redefine("Muon_pt_uncor", "Muon_pt_uncor[muoncuts]")
                       .Define("skim_PuppiMET_pt_mu_unsmear", "PuppiMET_pt")
                       .Define("skim_PuppiMET_phi_mu_unsmear", "PuppiMET_phi")
                       .Redefine("PuppiMET_phi", muonscaremetphi, {"Muon_pt_uncor", "Muon_pt", "Muon_phi", "skim_PuppiMET_pt_mu_unsmear", "skim_PuppiMET_phi_mu_unsmear"})
                       .Redefine("PuppiMET_pt",  muonscaremet,    {"Muon_pt_uncor", "Muon_pt", "Muon_phi", "skim_PuppiMET_pt_mu_unsmear", "skim_PuppiMET_phi_mu_unsmear"});
        }
    }
}

void NanoAODAnalyzerrdframe::applyBSFs(string btagYear, string btagMap, string btagMapLight, float btagcut) {
    cout << "Loading Btag SF (fixed WP, Method 1a): " << btagYear << "  comb=" << btagMap << "  light=" << btagMapLight << endl;
    // Fixed-WP Method 1a uses the multi-era uncertainty scheme (up/down_correlated and up/down_uncorrelated),
    // where the SF's own JES/JER calibration uncertainties are already included in quadrature inside
    // up/down_uncorrelated (btagSFbc_<era>); JME jet pT shifts propagate automatically via cleanjet4vecs.

    auto bSFreader = loadCorrectionSet("data/BTV/" + btagYear + "/btagging.json.gz");
    // Fixed-WP corrections:
    //   btagMap      = e.g. "particleNet_comb" (2022/23) or "UParTAK4_comb" (2024) — for b/c jets
    //   btagMapLight = e.g. "particleNet_light"           or "UParTAK4_light"       — for light jets
    // Inputs: (systematic, working_point, flavor, abseta, pt)  — no discriminant needed.
    auto _btagSF      = bSFreader->at(btagMap);
    auto _btagSFlight = bSFreader->at(btagMapLight);

    // Cross-validate btagcut against the official POG *_wp_values ("M") in btagging.json.gz
    if (btagcut <= 0.0f || !std::isfinite(btagcut)) {
        std::string msg = "[ERROR] applyBSFs: Invalid btagcut (" + std::to_string(btagcut) + ") for " + btagYear;
        std::cerr << msg << std::endl;
        throw std::runtime_error(msg);
    }
    const std::string wpValName = _isRun24 ? "UParTAK4_wp_values" : "particleNet_wp_values";
    float wpJsonVal = static_cast<float>(bSFreader->at(wpValName)->evaluate({"M"}));
    if (std::abs(btagcut - wpJsonVal) > 1e-5f) {
        std::string msg = "[ERROR] applyBSFs: btagcut from analysis_config.json (" + std::to_string(btagcut) +
                          ") mismatches " + wpValName + "(\"M\") = " + std::to_string(wpJsonVal) +
                          " in data/BTV/" + btagYear + "/btagging.json.gz";
        std::cerr << msg << std::endl;
        throw std::runtime_error(msg);
    }

    // -----------------------------------------------------------------------
    // BTV Method 1a (Event Reweighting) MC Efficiency Map Loading
    // Weight = Prod_{tagged} SF_i * Prod_{untagged} (1 - SF_j * eff_j) / (1 - eff_j)
    // Hierarchical search:
    //   1. Process-specific map: data/BTV/<ch>/<year>/btag_eff_<proc>.root (e.g. ttbar, qcd, wjets)
    //   2. Sample-specific map:  data/BTV/<ch>/<year>/btag_eff_<sample>.root
    //   3. Inclusive map:        data/BTV/<ch>/<year>/btag_eff.root (or btag_eff_inclusive.root)
    //   4. Missing/invalid map:  Explicit [ERROR] + std::runtime_error on evaluation
    // -----------------------------------------------------------------------
    std::shared_ptr<TH2> sp_eff_b = nullptr;
    std::shared_ptr<TH2> sp_eff_c = nullptr;
    std::shared_ptr<TH2> sp_eff_light = nullptr;
    bool has_eff = false;

    std::string baseName = _outfilename;
    size_t lastSlash = baseName.find_last_of("/\\");
    if (lastSlash != std::string::npos) baseName = baseName.substr(lastSlash + 1);
    size_t dotPos = baseName.find_last_of(".");
    if (dotPos != std::string::npos) baseName = baseName.substr(0, dotPos);

    // Strip leading "hist_" if present (as in plotIt files.yml hist_*.root)
    std::string cleanName = baseName;
    if (cleanName.rfind("hist_", 0) == 0) {
        cleanName = cleanName.substr(5);
    }

    std::string lowerName = cleanName;
    std::transform(lowerName.begin(), lowerName.end(), lowerName.begin(), ::tolower);

    // -----------------------------------------------------------------------
    // Process group classification aligned with plotIt files.yml:
    // -----------------------------------------------------------------------
    std::string procGroup = "";

    // 1. Signal (type: signal, groups: GLFVST*, GLFVTT*)
    //    Split into signal_tc (tc coupling: TCMuTau, TCETau, TTtoCMuTau, TTtoCETau)
    //    and signal_tu (tu coupling: TUMuTau, TUETau, TTtoUMuTau, TTtoUETau).
    //    MUST check before ttbar because TTto*LFV contains "ttto".
    if (lowerName.find("tcmutau") != std::string::npos || lowerName.find("tcetau") != std::string::npos ||
        lowerName.find("tttocmutau") != std::string::npos || lowerName.find("tttocetau") != std::string::npos ||
        lowerName.find("signal_tc") != std::string::npos) {
        procGroup = "signal_tc";
    }
    else if (lowerName.find("tumutau") != std::string::npos || lowerName.find("tuetau") != std::string::npos ||
             lowerName.find("tttoumutau") != std::string::npos || lowerName.find("tttouetau") != std::string::npos ||
             lowerName.find("signal_tu") != std::string::npos) {
        procGroup = "signal_tu";
    }
    else if (lowerName.find("lfv") != std::string::npos) {
        procGroup = "signal";
    }
    // 2. tt+X / ttV (group: GttV, yields-group: '7tt+X')
    //    e.g. TTWtoQQ, TTZtoQQ, TTHto2B, TTHtoNon2B
    //    MUST check before ttbar because TTW/TTZ/TTH start with "TT".
    else if (lowerName.find("ttw") != std::string::npos ||
             lowerName.find("ttz") != std::string::npos ||
             lowerName.find("tth") != std::string::npos) {
        procGroup = "ttx";
    }
    // 3. ttbar (groups: Gttll, Gttlj, Gttjj, yields-groups: '1\ttbar LL', '2\ttbar LJ', '3\ttbar Had')
    //    e.g. TTto2L2Nu, TTtoLNu2Q, TTto4Q
    else if (lowerName.find("ttto") != std::string::npos ||
             lowerName.find("ttbar") != std::string::npos ||
             lowerName.find("tt_") != std::string::npos) {
        procGroup = "ttbar";
    }
    // 4. Single Top (group: GSingleT, yields-group: '5Single Top')
    //    e.g. TBbar_s-channel, TbarB_s-channel, TQbarto*-t-channel, TbarQto*-t-channel,
    //         TWminusto*, TbarWplusto*, TBbarQ_t-channel, TbarBQ_t-channel
    else if (lowerName.find("tbbar") != std::string::npos ||
             lowerName.find("tbarb") != std::string::npos ||
             lowerName.find("tqbarto") != std::string::npos ||
             lowerName.find("tbarqto") != std::string::npos ||
             lowerName.find("twminus") != std::string::npos ||
             lowerName.find("tbarwplus") != std::string::npos ||
             lowerName.find("singletop") != std::string::npos ||
             lowerName.find("st_") != std::string::npos) {
        procGroup = "singletop";
    }
    // 5. W+Jets (group: GWJets, yields-group: '6W+Jets')
    //    2022/2023: WtoLNu-2Jets, WtoLNu-2Jets_0J, WtoLNu-2Jets_1J, WtoLNu-2Jets_2J
    //    2024:      WtoENu-4Jets, WtoMuNu-4Jets, WtoTauNu-4Jets, WtoLNu-4Jets_*
    else if (lowerName.find("wtolnu") != std::string::npos ||
             lowerName.find("wtoenu") != std::string::npos ||
             lowerName.find("wtomunu") != std::string::npos ||
             lowerName.find("wtotaunu") != std::string::npos ||
             lowerName.find("wjets") != std::string::npos) {
        procGroup = "wjets";
    }
    // 6. Z+Jets / Drell-Yan (group: GZJets, yields-group: '4Z+Jets')
    //    2022/2023: DYto2L-2Jets_MLL-10to50, DYto2L-2Jets_MLL-50
    //    2024:      DYto2E-2Jets_MLL-10-50, DYto2E-2Jets_MLL-50,
    //               DYto2Mu-2Jets_MLL-10-50, DYto2Mu-2Jets_MLL-50,
    //               DYto2Tau-2Jets_MLL-10-50, DYto2Tau-2Jets_MLL-50
    else if (lowerName.find("dyto2l") != std::string::npos ||
             lowerName.find("dyto2e") != std::string::npos ||
             lowerName.find("dyto2mu") != std::string::npos ||
             lowerName.find("dyto2tau") != std::string::npos ||
             lowerName.find("dyto") != std::string::npos ||
             lowerName.find("dyjets") != std::string::npos ||
             lowerName.find("zjets") != std::string::npos) {
        procGroup = "dyjets";
    }
    // 7. Diboson (group: GVV, yields-group: '8VV')
    //    e.g. WW, WZ, ZZ
    else if (lowerName == "ww" || lowerName == "wz" || lowerName == "zz" ||
             lowerName.find("ww_") != std::string::npos || lowerName.find("wz_") != std::string::npos || lowerName.find("zz_") != std::string::npos ||
             lowerName.find("_ww") != std::string::npos || lowerName.find("_wz") != std::string::npos || lowerName.find("_zz") != std::string::npos ||
             lowerName.find("diboson") != std::string::npos) {
        procGroup = "diboson";
    }
    // 8. QCD (group: GQCD, yields-group: '9QCD')
    //    e.g. QCD_Pt*MuEnriched, QCD_Pt*EMEnriched, QCD_Pt*
    else if (lowerName.find("qcd") != std::string::npos) {
        procGroup = "qcd";
    }

    std::string chName = _isMuonCh ? "muon" : "electron";
    std::vector<std::string> yearList = {btagYear};
    if (!_year.empty() && _year != btagYear) {
        yearList.push_back(_year);
    }

    std::vector<std::string> effCandidates;

    // 1. Channel-specific candidates (highest priority)
    for (const auto &y : yearList) {
        if (!procGroup.empty() && procGroup != "inclusive") {
            effCandidates.push_back("data/BTV/" + chName + "/" + y + "/btag_eff_" + procGroup + ".root");
        }
        if (!cleanName.empty() && cleanName != procGroup) {
            effCandidates.push_back("data/BTV/" + chName + "/" + y + "/btag_eff_" + cleanName + ".root");
        }
        if (procGroup.empty() || procGroup == "inclusive") {
            effCandidates.push_back("data/BTV/" + chName + "/" + y + "/btag_eff.root");
            effCandidates.push_back("data/BTV/" + chName + "/" + y + "/btag_eff_inclusive.root");
        }
    }

    // 2. Legacy / channel-agnostic fallback candidates
    for (const auto &y : yearList) {
        if (!procGroup.empty() && procGroup != "inclusive") {
            effCandidates.push_back("data/BTV/" + y + "/btag_eff_" + procGroup + ".root");
        }
        if (!cleanName.empty() && cleanName != procGroup) {
            effCandidates.push_back("data/BTV/" + y + "/btag_eff_" + cleanName + ".root");
        }
        if (procGroup.empty() || procGroup == "inclusive") {
            effCandidates.push_back("data/BTV/" + y + "/btag_eff.root");
            effCandidates.push_back("data/BTV/" + y + "/btag_eff_inclusive.root");
        }
    }

    std::string chosenEffPath = "";
    for (const auto &candPath : effCandidates) {
        if (gSystem->AccessPathName(candPath.c_str()) != 0) continue; // File does not exist
        TFile *fEff = TFile::Open(candPath.c_str(), "READ");
        if (fEff && !fEff->IsZombie()) {
            TH2 *h_b = dynamic_cast<TH2*>(fEff->Get("h2_eff_b"));
            TH2 *h_c = dynamic_cast<TH2*>(fEff->Get("h2_eff_c"));
            TH2 *h_l = dynamic_cast<TH2*>(fEff->Get("h2_eff_light"));
            if (h_b && h_c && h_l) {
                h_b->SetDirectory(0);
                h_c->SetDirectory(0);
                h_l->SetDirectory(0);
                sp_eff_b = std::shared_ptr<TH2>(h_b);
                sp_eff_c = std::shared_ptr<TH2>(h_c);
                sp_eff_light = std::shared_ptr<TH2>(h_l);
                has_eff = true;
                chosenEffPath = candPath;
                fEff->Close();
                delete fEff;
                break;
            }
            fEff->Close();
            delete fEff;
        }
    }

    if (has_eff) {
        cout << "[INFO] Loaded BTV Method 1a efficiency maps from " << chosenEffPath
             << " (process group: " << (procGroup.empty() ? "unclassified" : procGroup) << ")" << endl;
    } else {
        cerr << "[WARN] No b-tag efficiency map found for " << btagYear
             << " in candidates (btag_eff_" << procGroup << ".root / btag_eff.root). "
             << "Evaluating btagWeight will throw an error unless running BTagEfficiencyAnalyzer." << endl;
    }

    // -----------------------------------------------------------------------
    // BTV POG Multi-era Standard Scheme (correlated + uncorrelated):
    //   0: central
    //   1/2: up/down_correlated   (legacy combined HF+LF, preserved)
    //   3/4: up/down_uncorrelated (legacy combined HF+LF, preserved)
    //   5/6: btagSFbc_correlated up/down     (HF _comb up/down_correlated, LF fixed to central)
    //   7/8: btagSFbc_<era> up/down          (HF _comb up/down_uncorrelated, LF fixed to central)
    //   9/10: btagSFlight_correlated up/down (LF _light up/down_correlated, HF fixed to central)
    //   11/12: btagSFlight_<era> up/down     (LF _light up/down_uncorrelated, HF fixed to central)
    // -----------------------------------------------------------------------
    const std::string wp = "M";  // Medium WP
    const std::vector<std::string> systs_btag = {
        "central",
        "up_correlated",   "down_correlated",    // 1/2
        "up_uncorrelated", "down_uncorrelated"   // 3/4
    };
    const int n_base_systs = static_cast<int>(systs_btag.size()); // 5
    const int n_systs = 13;

    // Helper lambda to query efficiency from 2D map with boundary clamping
    auto get_efficiency = [has_eff, sp_eff_b, sp_eff_c, sp_eff_light, chosenEffPath](int flav, float pt, float abseta) -> float {
        if (!has_eff) {
            std::string msg = "[ERROR] get_efficiency: BTV Method 1a efficiency map is missing.";
            std::cerr << msg << std::endl;
            throw std::runtime_error(msg);
        }
        TH2 *h = nullptr;
        if (flav == 5) {
            h = sp_eff_b.get();
        } else if (flav == 4) {
            h = sp_eff_c.get();
        } else {
            h = sp_eff_light.get();
        }
        if (!h) {
            std::string msg = "[ERROR] get_efficiency: Null TH2 efficiency histogram for flav=" + std::to_string(flav) +
                              " in " + chosenEffPath;
            std::cerr << msg << std::endl;
            throw std::runtime_error(msg);
        }

        float minPt  = h->GetXaxis()->GetXmin();
        float maxPt  = h->GetXaxis()->GetXmax() - 1e-3f;
        float minEta = h->GetYaxis()->GetXmin();
        float maxEta = h->GetYaxis()->GetXmax() - 1e-3f;

        float clampedPt  = std::clamp(pt, minPt, maxPt);
        float clampedEta = std::clamp(abseta, minEta, maxEta);

        int bx = h->GetXaxis()->FindFixBin(clampedPt);
        int by = h->GetYaxis()->FindFixBin(clampedEta);
        float val = h->GetBinContent(bx, by);
        if (!std::isfinite(val) || val <= 0.0f || val >= 1.0f) {
            std::string msg = "[ERROR] get_efficiency: Invalid BTV MC efficiency value (" + std::to_string(val) +
                              ") at flav=" + std::to_string(flav) + ", pt=" + std::to_string(pt) +
                              ", abseta=" + std::to_string(abseta) + " (bin " + std::to_string(bx) + "," +
                              std::to_string(by) + ") in " + chosenEffPath +
                              ". Method 1a requires 0 < eff < 1.";
            std::cerr << msg << std::endl;
            throw std::runtime_error(msg);
        }
        return val;
    };

    auto btagSF_fixedWP = [this, _btagSF, _btagSFlight, wp, systs_btag, n_base_systs, n_systs, btagcut, has_eff, btagYear, procGroup, get_efficiency]
                          (floats &pts, floats &etas, uchars &hadflav, floats &bscores) -> floats {
        // Returns event-level weight vector [n_systs = 13].
        // Under BTV Method 1a:
        //   Tagged jets: multiply SF
        //   Untagged jets: multiply (1 - SF * eff) / (1 - eff)
        if (!has_eff) {
            std::string msg = "[ERROR] btagSF_fixedWP: No b-tag efficiency map found for " + btagYear +
                              " (process group: " + procGroup + "). Method 1a requires efficiency maps.";
            std::cerr << msg << std::endl;
            throw std::runtime_error(msg);
        }
        floats out(n_systs, 1.0f);
        size_t n_jets = std::min({pts.size(), etas.size(), hadflav.size(), bscores.size()});
        for (size_t j = 0; j < n_jets; j++) {
            float pt     = pts[j];
            float abseta = std::abs(etas[j]);
            int   flav   = static_cast<int>(hadflav[j]);
            float score  = bscores[j];

            // Skip jets outside kinematic range
            if (pt < 20.f || abseta >= 2.5f) continue;
            bool is_tagged = (score > btagcut);

            float eff = get_efficiency(flav, pt, abseta);
            float denom = 1.0f - eff;

            float w_jet[5];
            for (int s = 0; s < n_base_systs; s++) {
                float sf = (flav == 0)
                    ? _btagSFlight->evaluate({systs_btag[s], wp, flav, abseta, pt})
                    : _btagSF->evaluate({systs_btag[s], wp, flav, abseta, pt});

                if (is_tagged) {
                    w_jet[s] = sf;
                } else {
                    // Method 1a: Untagged jet factor
                    w_jet[s] = (1.0f - sf * eff) / denom;
                }
            }

            // 0: central, 1/2: combined corr up/dn, 3/4: combined uncorr up/dn
            for (int s = 0; s < n_base_systs; s++) {
                out[s] *= w_jet[s];
            }

            if (flav != 0) {
                // Heavy-flavour (b/c) jet: varies for HF systematics [5..8], central for LF [9..12]
                out[5]  *= w_jet[1];
                out[6]  *= w_jet[2];
                out[7]  *= w_jet[3];
                out[8]  *= w_jet[4];
                out[9]  *= w_jet[0];
                out[10] *= w_jet[0];
                out[11] *= w_jet[0];
                out[12] *= w_jet[0];
            } else {
                // Light-flavour (udsg) jet: central for HF [5..8], varies for LF systematics [9..12]
                out[5]  *= w_jet[0];
                out[6]  *= w_jet[0];
                out[7]  *= w_jet[0];
                out[8]  *= w_jet[0];
                out[9]  *= w_jet[1];
                out[10] *= w_jet[2];
                out[11] *= w_jet[3];
                out[12] *= w_jet[4];
            }
        }
        return out;
    };

    // btagWeight[0] = central, [1/2] = corr up/dn, [3/4] = uncorr up/dn,
    // [5/6] = btagSFbc_correlated up/dn, [7/8] = btagSFbc_<era> up/dn,
    // [9/10] = btagSFlight_correlated up/dn, [11/12] = btagSFlight_<era> up/dn
    const std::string bTagCol = _isRun24 ? "Jet_btagUParTAK4B" : "Jet_btagPNetB";
    if (isDefined("btagWeight")) {
        _rlm = _rlm.Redefine("btagWeight", btagSF_fixedWP,
                             {"Jet_pt", "Jet_eta", "Jet_hadronFlavour", bTagCol});
    } else {
        _rlm = _rlm.Define("btagWeight", btagSF_fixedWP,
                           {"Jet_pt", "Jet_eta", "Jet_hadronFlavour", bTagCol});
    }
    if (isDefined("btagWeight_loose")) {
        _rlm = _rlm.Redefine("btagWeight_loose", btagSF_fixedWP,
                             {"Jet_pt_loose", "Jet_eta_loose", "Jet_hadronFlavour_loose", bTagCol + "_loose"});
    } else {
        _rlm = _rlm.Define("btagWeight_loose", btagSF_fixedWP,
                           {"Jet_pt_loose", "Jet_eta_loose", "Jet_hadronFlavour_loose", bTagCol + "_loose"});
    }
    // Note: btagNormalization is NOT needed for fixed-WP SFs.
    // The direct event weight btagWeight[0] is used in defineWeightVars().

    // =======================================================================
    // SHAPE CORRECTION (Case 3) — disabled; restore when needed
    // Use this block instead of the fixed-WP block above when the analysis
    // uses the b-tag discriminant distribution directly in the fit (e.g. DNN
    // input, template fit on discriminant). Requires:
    //   btagMap = "particleNet_shape" (or "deepJet_shape" / "robustParticleTransformer_shape")
    //   After enabling, also call defineBTagNormalization() in defineWeightVars()
    //   and switch all btagWeight[*] references back to btagWeightNorm[*].
    //
    // Systematic indices for shape correction (17 entries):
    //   0: central
    //   1/2:  up/down_hf          (b-jets, HF)
    //   3/4:  up/down_lf          (light jets, LF)
    //   5/6:  up/down_hfstats1    (b-jets stat1)
    //   7/8:  up/down_hfstats2    (b-jets stat2)
    //   9/10: up/down_lfstats1    (light stat1)
    //  11/12: up/down_lfstats2    (light stat2)
    //  13/14: up/down_cferr1      (c-jets err1)
    //  15/16: up/down_cferr2      (c-jets err2)
    // =======================================================================
    /*
    auto btagSF_shape = [this, _btagSF](floats &pts, floats &etas, uchars &hadflav, floats &btags) -> floatsVec {
        std::vector<std::string> systs = {"central",
            "up_hf",       "down_hf",
            "up_lf",       "down_lf",
            "up_hfstats1", "down_hfstats1",
            "up_hfstats2", "down_hfstats2",
            "up_lfstats1", "down_lfstats1",
            "up_lfstats2", "down_lfstats2",
            "up_cferr1",   "down_cferr1",
            "up_cferr2",   "down_cferr2"};

        floats wVec;
        wVec.reserve(systs.size());
        floatsVec out;
        out.reserve(pts.size());

        for (auto i = 0; i < int(pts.size()); i++) {
            for (auto &syst : systs) {
                float sf = 1.0f;
                if (pts[i] < 30 || std::abs(etas[i]) > 2.5) {
                    sf = 1.0f;
                } else if (btags[i] < 0.0f) {
                    // disc < 0 (e.g. btagPNetB=-1): jet is untaggable/masked.
                    // Correctionlib binning starts at 0; evaluating would throw.
                    sf = 1.0f;
                } else if (syst.find("cferr") != std::string::npos && hadflav[i] != 4) {
                    sf = 1.0f;  // cferr systs apply only to c-jets (flav==4)
                } else if ((syst.find("hf") != std::string::npos || syst.find("lf") != std::string::npos)
                           && hadflav[i] == 4) {
                    sf = 1.0f;  // hf/lf systs do NOT apply to c-jets
                } else {
                    try {
                        sf = _btagSF->evaluate({syst, int(hadflav[i]), std::abs(etas[i]), pts[i], btags[i]});
                    } catch (const std::exception &e) {
                        std::cerr << "[WARN] btagSF shape evaluate failed"
                                  << " syst=" << syst << " flav=" << int(hadflav[i])
                                  << " eta=" << etas[i] << " pt=" << pts[i]
                                  << " disc=" << btags[i] << ": " << e.what() << "\n";
                        sf = 1.0f;
                    }
                }
                wVec.emplace_back(sf);
            }
            out.emplace_back(wVec);
            wVec.clear();
        }
        return out;
    };

    auto btag_evWeight = [this](floatsVec &btagWeights) -> floats {
        const int vars = 17;
        floats out(vars, 1.0f);
        for (auto &jet : btagWeights) {
            for (int i = 0; i < vars; i++) {
                out[i] *= jet[i];
            }
        }
        return out;
    };

    // Shape SF needs the per-jet floatsVec column (Jet_btagSF) as intermediate.
    // The final btagWeight[17] is the per-event product.
    // After enabling, also restore defineBTagNormalization() call and
    // use btagWeightNorm[*] (not btagWeight[*]) in defineWeightVars().
    _rlm = _rlm.Define("btag_var",     []()        { return strings(btag_var); })
               .Define("btag_jes_var", [jes_var]() { return strings(jes_var); })
               .Define("Jet_btagSF",   btagSF_shape, {"Jet_pt", "Jet_eta", "Jet_hadronFlavour", "Jet_btagPNetB"})
               .Define("btagWeight",   btag_evWeight, {"Jet_btagSF"});
    */
}

void NanoAODAnalyzerrdframe::selectTaus(string cut, string tauYear, string vsjet, string vsmu, string vse) {
    auto syst_unc = _syst;

    //TES var.
    _rlm = _rlm.Define("Tau_pt_nom_tes", "Tau_pt");
    if (!_isData) {
        // In the Process stage, Tau_pt has already been corrected by Tau_ES_nom in the skim.
        // Tau_ES_up and Tau_ES_down are stored as flat float vectors in the skim tree.
        // To apply the up/down variation: rescale by Tau_ES_up/Tau_ES_nom or Tau_ES_down/Tau_ES_nom.
        // Guard against Tau_ES_nom == 0 (unshifted jet) and against old skim files that
        // pre-date Tau_ES_* columns.
        const bool hasTauES = isDefined("Tau_ES_nom") && isDefined("Tau_ES_up") && isDefined("Tau_ES_down");
        if (!hasTauES && _syst.find("tes") != std::string::npos) {
            std::cerr << "[ERROR] selectTaus: syst=" << _syst
                      << " requested, but Tau_ES_nom/up/down branches are not found in input skim! "
                      << "Skim must contain Tau_ES_* branches to apply TES variation." << std::endl;
            throw std::runtime_error("[ERROR] [selectTaus] Missing Tau_ES_nom/up/down in skimmed input for " + _syst);
        }
        if (hasTauES && _syst.find("tesup") != std::string::npos) {
            _rlm = _rlm.Define("Tau_pt_unc_toapply",
                               [](floats &es_nom, floats &es_up) -> floats {
                                   floats out;
                                   out.reserve(es_nom.size());
                                   for (size_t i = 0; i < es_nom.size(); i++) {
                                       float ratio = (es_nom[i] > 0.f) ? (es_up[i] / es_nom[i]) : 1.f;
                                       out.emplace_back(ratio);
                                   }
                                   return out;
                               },
                               {"Tau_ES_nom", "Tau_ES_up"})
                       .Redefine("Tau_pt",   "Tau_pt   * Tau_pt_unc_toapply")
                       .Redefine("Tau_mass", "Tau_mass * Tau_pt_unc_toapply");
            std::cout << "[TES] Applied Tau_ES_up/Tau_ES_nom rescaling (tesup)\n";
        } else if (hasTauES && _syst.find("tesdown") != std::string::npos) {
            _rlm = _rlm.Define("Tau_pt_unc_toapply",
                               [](floats &es_nom, floats &es_dn) -> floats {
                                   floats out;
                                   out.reserve(es_nom.size());
                                   for (size_t i = 0; i < es_nom.size(); i++) {
                                       float ratio = (es_nom[i] > 0.f) ? (es_dn[i] / es_nom[i]) : 1.f;
                                       out.emplace_back(ratio);
                                   }
                                   return out;
                               },
                               {"Tau_ES_nom", "Tau_ES_down"})
                       .Redefine("Tau_pt",   "Tau_pt   * Tau_pt_unc_toapply")
                       .Redefine("Tau_mass", "Tau_mass * Tau_pt_unc_toapply");
            std::cout << "[TES] Applied Tau_ES_down/Tau_ES_nom rescaling (tesdown)\n";
        }
        // nominal: Tau_pt already includes Tau_ES_nom from the skim — no further action needed.
    }


    auto overlap_removal_leptau = [](FourVectorVec &lep4vecs, FourVectorVec &tau4vecs) {
        ints out;
        for (auto tau: tau4vecs) {
            int check = 0;
            for (auto lep: lep4vecs) {
                auto dR = ROOT::Math::VectorUtil::DeltaR(lep, tau);
                if( dR >= 0.4 ) check = 1;
            }
            out.emplace_back(check);
        }
        return out;
    };

    _rlm = _rlm.Define("tau4vecs", ::gen4vec, {"Tau_pt", "Tau_eta", "Tau_phi", "Tau_mass"});
    _rlm = _rlm.Define("leptauoverlap", overlap_removal_leptau, {"lep4vecs","tau4vecs"});

    // input vector: vec[pt][vars]
    auto skimCol = [this](floatsVec toSkim, ints cut)->floatsVec {

        floatsVec out;
        size_t n = std::min(toSkim.size(), cut.size());
        for (size_t i=0; i<n; i++) {
            if (cut[i] > 0) out.emplace_back(toSkim[i]);
        }
        if (out.size() > 0) return out;
        else {
            size_t vec_len = (toSkim.size() > 0 && toSkim[0].size() > 0) ? toSkim[0].size() : 39;
            return {floats(vec_len, 1.0f)};
        }
    };

    // Fake factor study - loose but not tight
    // Hadronic Tau Object Selections
    string deeptauidcuts = "Tau_idDeepTau2018v2p5VSmu >= " + vsmu + " && Tau_idDeepTau2018v2p5VSe >= " + vse + " && Tau_idDeepTau2018v2p5VSjet >= " + vsjet;
    string deeptauidcuts_loose = "Tau_idDeepTau2018v2p5VSmu >= " + vsmu + " && Tau_idDeepTau2018v2p5VSe >= " + vse + " && Tau_idDeepTau2018v2p5VSjet >= 4";
    _rlm = _rlm.Define("taucuts", cut)
               .Define("deeptauidcuts", deeptauidcuts)
               .Define("deeptauidcuts_loose", deeptauidcuts_loose);

    // Hadronic Tau Selection
    _rlm = _rlm.Define("seltaucuts_loose","taucuts && deeptauidcuts_loose && leptauoverlap")
               .Define("Tau_pt_loose", "Tau_pt[seltaucuts_loose]")
               .Define("Tau_pt_loose_gen", "Tau_pt[seltaucuts_loose]")
               .Define("Tau_eta_loose", "Tau_eta[seltaucuts_loose]")
               .Define("Tau_phi_loose", "Tau_phi[seltaucuts_loose]")
               .Define("Tau_mass_loose", "Tau_mass[seltaucuts_loose]")
               .Define("Tau_charge_loose", "Tau_charge[seltaucuts_loose]")
               .Define("Tau_decayMode_loose", "Tau_decayMode[seltaucuts_loose]")
               .Define("nloosetaupass", "int(Tau_pt_loose.size())")
               .Define("cleanloosetau4vecs", ::gen4vec, {"Tau_pt_loose", "Tau_eta_loose", "Tau_phi_loose", "Tau_mass_loose"});


    _rlm = _rlm.Define("seltaucuts","taucuts && deeptauidcuts && leptauoverlap")
               .Define("Tau_pt_gen", "Tau_pt[seltaucuts]")
               .Define("sel_Tau_pt_nom", "Tau_pt_nom_tes[seltaucuts]")
               .Redefine("Tau_pt", "Tau_pt[seltaucuts]")
               .Redefine("Tau_eta", "Tau_eta[seltaucuts]")
               .Redefine("Tau_phi", "Tau_phi[seltaucuts]")
               .Redefine("Tau_mass", "Tau_mass[seltaucuts]")
               .Redefine("Tau_charge", "Tau_charge[seltaucuts]")
               .Redefine("Tau_jetIdx", "Tau_jetIdx[seltaucuts]")
               .Redefine("Tau_decayMode", "Tau_decayMode[seltaucuts]")
               .Define("ncleantaupass", "int(Tau_pt.size())")
               .Define("cleantau4vecs", ::gen4vec, {"Tau_pt", "Tau_eta", "Tau_phi", "Tau_mass"});

    if (!_isData && (_syst.find("tesup") != std::string::npos || _syst.find("tesdown") != std::string::npos)) {
        auto tesMetPt = [](const floats &pts_nom, const floats &pts_tes, const floats &phis, float met, float metphi) -> float {
            double metx = met * std::cos(metphi);
            double mety = met * std::sin(metphi);
            size_t n = std::min({pts_nom.size(), pts_tes.size(), phis.size()});
            for (size_t i = 0; i < n; i++) {
                double dpx = (pts_tes[i] - pts_nom[i]) * std::cos(phis[i]);
                double dpy = (pts_tes[i] - pts_nom[i]) * std::sin(phis[i]);
                metx -= dpx;
                mety -= dpy;
            }
            return static_cast<float>(std::sqrt(metx * metx + mety * mety));
        };

        auto tesMetPhi = [](const floats &pts_nom, const floats &pts_tes, const floats &phis, float met, float metphi) -> float {
            double metx = met * std::cos(metphi);
            double mety = met * std::sin(metphi);
            size_t n = std::min({pts_nom.size(), pts_tes.size(), phis.size()});
            for (size_t i = 0; i < n; i++) {
                double dpx = (pts_tes[i] - pts_nom[i]) * std::cos(phis[i]);
                double dpy = (pts_tes[i] - pts_nom[i]) * std::sin(phis[i]);
                metx -= dpx;
                mety -= dpy;
            }
            return static_cast<float>(std::atan2(mety, metx));
        };

        _rlm = _rlm.Define("PuppiMET_pt_tes_uncorr",  "PuppiMET_pt")
                   .Define("PuppiMET_phi_tes_uncorr", "PuppiMET_phi")
                   .Redefine("PuppiMET_phi", tesMetPhi, {"sel_Tau_pt_nom", "Tau_pt", "Tau_phi", "PuppiMET_pt_tes_uncorr", "PuppiMET_phi_tes_uncorr"})
                   .Redefine("PuppiMET_pt",  tesMetPt,  {"sel_Tau_pt_nom", "Tau_pt", "Tau_phi", "PuppiMET_pt_tes_uncorr", "PuppiMET_phi_tes_uncorr"});

        if (isDefined("MET_pt") && isDefined("MET_phi")) {
            _rlm = _rlm.Define("MET_pt_tes_uncorr",  "MET_pt")
                       .Define("MET_phi_tes_uncorr", "MET_phi")
                       .Redefine("MET_phi", tesMetPhi, {"sel_Tau_pt_nom", "Tau_pt", "Tau_phi", "MET_pt_tes_uncorr", "MET_phi_tes_uncorr"})
                       .Redefine("MET_pt",  tesMetPt,  {"sel_Tau_pt_nom", "Tau_pt", "Tau_phi", "MET_pt_tes_uncorr", "MET_phi_tes_uncorr"});
        }
        std::cout << "[TES] Propagated tau energy scale variation to PuppiMET (and MET) for selected taus (" << _syst << ")\n";
    }

    if (!_isData) {
        // Tau SF
        cout << "Loading Tau SF" << endl;

        // Convert WP integer strings ("1".."8" / "1".."4") to CorrectionLib string WPs:
        // VSjet and VSe (1-8): 1=VVVLoose, 2=VVLoose, 3=VLoose, 4=Loose, 5=Medium, 6=Tight, 7=VTight, 8=VVTight
        // VSmu (1-4):          1=VLoose, 2=Loose, 3=Medium, 4=Tight
        auto wpToNameVSjetOrEle = [](const std::string& wp) -> std::string {
            if (wp == "1" || wp == "VVVLoose") return "VVVLoose";
            if (wp == "2" || wp == "VVLoose")  return "VVLoose";
            if (wp == "3" || wp == "VLoose")   return "VLoose";
            if (wp == "4" || wp == "Loose")    return "Loose";
            if (wp == "5" || wp == "Medium")   return "Medium";
            if (wp == "6" || wp == "Tight")    return "Tight";
            if (wp == "7" || wp == "VTight")   return "VTight";
            if (wp == "8" || wp == "VVTight")  return "VVTight";
            return wp;
        };

        auto wpToNameVSmu = [](const std::string& wp) -> std::string {
            if (wp == "1" || wp == "VLoose")  return "VLoose";
            if (wp == "2" || wp == "Loose")   return "Loose";
            if (wp == "3" || wp == "Medium")  return "Medium";
            if (wp == "4" || wp == "Tight")   return "Tight";
            return wp;
        };

        std::string tauid_vsjet = wpToNameVSjetOrEle(vsjet);
        std::string tauid_vsmu  = wpToNameVSmu(vsmu);
        std::string tauid_vse   = wpToNameVSjetOrEle(vse);

        cout << "Tau ID SF WP vsJet : " << tauid_vsjet << endl;
        cout << "Tau ID SF WP vsMuon : " << tauid_vsmu << endl;
        cout << "Tau ID SF WP vsElectron : " << tauid_vse << endl;

        std::string tauFolder = tauYear;
        if (tauYear == "2022_preEE")       tauFolder = "2022_Summer22";
        else if (tauYear == "2022_postEE")  tauFolder = "2022_Summer22EE";
        else if (tauYear == "2023_preBPix") tauFolder = "2023_Summer23";
        else if (tauYear == "2023_postBPix")tauFolder = "2023_Summer23BPix";
        else if (tauYear == "2024")        tauFolder = "2024_Summer24";

        auto tauSFreader = loadCorrectionSet("data/TauIDSFs/" + tauFolder + "/tau.json.gz");
        auto _tauidSFjet = tauSFreader->at("DeepTau2018v2p5VSjet");
        auto _tauidSFele = tauSFreader->at("DeepTau2018v2p5VSe");
        auto _tauidSFmu  = tauSFreader->at("DeepTau2018v2p5VSmu");

        // Dynamic WP fallback for DeepTau2018v2p5VSjet:
        // Current 2022/2023 tau POG JSONs provide only 'Medium' WP, while 2024 has 'Loose', 'Medium', 'Tight', 'VTight'.
        // If the requested WP (e.g. 'VTight') is not present in the JSON, probe throws.
        // For 2022/2023, this is an expected official limitation of the central Tau POG release, so we fall back to 'Medium'.
        // For 2024 (or if 'Medium' also fails), an explicit [ERROR] is thrown per Principle 4.
        std::string active_tauid_vsjet = tauid_vsjet;
        try {
            _tauidSFjet->evaluate({40.0, 0, 5, active_tauid_vsjet, tauid_vse, "nom", "dm"});
            std::cout << "[INFO] DeepTau2018v2p5VSjet: successfully loaded requested WP '" << active_tauid_vsjet << "'" << std::endl;
        } catch (const std::exception& e) {
            if (tauYear != "2024") {
                std::cout << "[INFO] DeepTau2018v2p5VSjet: requested WP '" << active_tauid_vsjet
                          << "' not provided in central " << tauFolder << " release (only 'Medium' is available). "
                          << "Using 'Medium' SF per Tau POG recommendation." << std::endl;
                active_tauid_vsjet = "Medium";
                try {
                    _tauidSFjet->evaluate({40.0, 0, 5, active_tauid_vsjet, tauid_vse, "nom", "dm"});
                } catch (const std::exception& e2) {
                    std::cerr << "[ERROR] DeepTau2018v2p5VSjet: fallback 'Medium' WP also failed in " << tauFolder << ": " << e2.what() << std::endl;
                    throw std::runtime_error("DeepTau2018v2p5VSjet failed for fallback Medium WP in " + tauFolder + ": " + e2.what());
                }
            } else {
                std::cerr << "[ERROR] DeepTau2018v2p5VSjet: requested WP '" << active_tauid_vsjet
                          << "' failed in " << tauFolder << ": " << e.what() << std::endl;
                throw std::runtime_error("DeepTau2018v2p5VSjet failed for requested WP " + active_tauid_vsjet + " in " + tauFolder + ": " + e.what());
            }
        }

        std::vector<std::string> systNamesVsJet = {
            "nom",                                     // 0
            "up", "down",                              // 1, 2
            "stat1_dm0_up", "stat1_dm0_down",          // 3, 4
            "stat1_dm1_up", "stat1_dm1_down",          // 5, 6
            "stat1_dm10_up", "stat1_dm10_down",        // 7, 8
            "stat1_dm11_up", "stat1_dm11_down",        // 9, 10
            "stat2_dm0_up", "stat2_dm0_down",          // 11, 12
            "stat2_dm1_up", "stat2_dm1_down",          // 13, 14
            "stat2_dm10_up", "stat2_dm10_down",        // 15, 16
            "stat2_dm11_up", "stat2_dm11_down",        // 17, 18
            "syst_alleras_up", "syst_alleras_down",    // 19, 20
            "syst_" + tauYear + "_up", "syst_" + tauYear + "_down", // 21, 22
            "syst_TES_" + tauYear + "_dm0_up", "syst_TES_" + tauYear + "_dm0_down",   // 23, 24
            "syst_TES_" + tauYear + "_dm1_up", "syst_TES_" + tauYear + "_dm1_down",   // 25, 26
            "syst_TES_" + tauYear + "_dm10_up", "syst_TES_" + tauYear + "_dm10_down", // 27, 28
            "syst_TES_" + tauYear + "_dm11_up", "syst_TES_" + tauYear + "_dm11_down", // 29, 30
            "stat_highpT_bin1_up", "stat_highpT_bin1_down", // 31, 32
            "stat_highpT_bin2_up", "stat_highpT_bin2_down", // 33, 34
            "syst_highpT_up", "syst_highpT_down",           // 35, 36
            "syst_highpT_extrap_up", "syst_highpT_extrap_down" // 37, 38
        };
        const size_t nSystVsJet = systNamesVsJet.size(); // 39

        // [TAU-02] 2024_Summer24 DeepTau2018v2p5VSjet (flag="dm") schema:
        //   2022/2023: per-DM TFormula with decomposed keys (stat1/stat2_dm*, syst_alleras, syst_<era>, syst_TES_<era>_dm*) + 'default'.
        //   2024     : per-DM pT-binned [20,40,60,200] constants with ONLY 'nom', 'up', 'down' (no decomposed keys, no 'default').
        // For 2024, only 'nom'/'up'/'down' are evaluated. [TAU-01] option A-1: the DM-binned 2024 up/down is routed into the
        // matching stat1_dm<DM> slot (DM 0 -> [3,4], 1 -> [5,6], 10 -> [7,8], 11 -> [9,10]) so that the 4 DMs are treated
        // as uncorrelated nuisances (tauidjetStat1dm{0,1,10,11}2024). All other decomposed slots stay at nom_sf in 2024.
        const bool isVsJetSchema2024 = (tauYear == "2024");
        std::cout << "[INFO] DeepTau2018v2p5VSjet: " << (isVsJetSchema2024
                      ? "2024 schema (nom/up/down per DM -> stat1_dm<DM> slots)"
                      : "2022/2023 schema (decomposed DM-binned sub-sources)") << std::endl;

        // [TAU-01/02] Map HPS decay mode -> first index of its stat1_dm<DM> up/down slot in systNamesVsJet.
        auto stat1SlotForDM = [](int dmv) -> int {
            if (dmv == 0)  return 3;
            if (dmv == 1)  return 5;
            if (dmv == 10) return 7;
            if (dmv == 11) return 9;
            return -1;
        };

        auto tauSFIdVsJet = [this, _tauidSFjet, active_tauid_vsjet, tauid_vse, systNamesVsJet, nSystVsJet, isVsJetSchema2024, stat1SlotForDM](floats &pt, uchars &dm, uchars &genmatch)->floatsVec {
            // CMS tau POG: DeepTauVSjet SF applies ONLY to real hadronic taus (genmatch==5).
            floats uncSources;
            uncSources.reserve(nSystVsJet);
            floatsVec wVec;
            size_t n = std::min({pt.size(), dm.size(), genmatch.size()});
            wVec.reserve(n > 0 ? n : 1);

            if (n > 0) {
                for (size_t i=0; i<n; i++) {
                    int gm = int(genmatch[i]);
                    if (gm == 5) {  // real hadronic tau only
                        int dm_val = int(dm[i]);
                        // [TAU-03] CMS Tau POG DeepTau2018v2p5 only supports HPS decay modes 0, 1, 10, 11.
                        if (dm_val != 0 && dm_val != 1 && dm_val != 10 && dm_val != 11) {
                            std::cerr << "[ERROR] tauSFIdVsJet: unsupported decay mode dm=" << dm_val
                                      << " for genmatch=5 (expected 0, 1, 10, 11; pt=" << pt[i] << ")" << std::endl;
                            throw std::runtime_error("tauSFIdVsJet: unsupported decay mode dm=" + std::to_string(dm_val));
                        }
                        float nom_sf = 1.0f;
                        try {
                            nom_sf = _tauidSFjet->evaluate({pt[i], dm_val, gm, active_tauid_vsjet, tauid_vse, "nom", "dm"});
                        } catch (const std::exception &e) {
                            std::cerr << "[ERROR] tauSFIdVsJet evaluate failed (pt=" << pt[i] << ", dm=" << dm_val
                                      << ", gm=" << gm << ", vsjet=" << active_tauid_vsjet << ", vse=" << tauid_vse << "): " << e.what() << std::endl;
                            throw std::runtime_error("tauSFIdVsJet evaluate failed (pt=" + std::to_string(pt[i]) +
                                                     ", dm=" + std::to_string(dm_val) + ", gm=" + std::to_string(gm) +
                                                     ", vsjet=" + active_tauid_vsjet + ", vse=" + tauid_vse + "): " + e.what());
                        }
                        uncSources.assign(nSystVsJet, nom_sf);
                        uncSources[0] = nom_sf;

                        if (isVsJetSchema2024) {
                            // [TAU-02] 2024: JSON provides only 'nom'/'up'/'down' per (DM, pT bin).
                            float up_sf = nom_sf, dn_sf = nom_sf;
                            try {
                                up_sf = _tauidSFjet->evaluate({pt[i], dm_val, gm, active_tauid_vsjet, tauid_vse, "up",   "dm"});
                                dn_sf = _tauidSFjet->evaluate({pt[i], dm_val, gm, active_tauid_vsjet, tauid_vse, "down", "dm"});
                            } catch (const std::exception &e) {
                                std::cerr << "[ERROR] tauSFIdVsJet (2024) evaluate failed for up/down (pt=" << pt[i]
                                          << ", dm=" << dm_val << ", gm=" << gm << "): " << e.what() << std::endl;
                                throw std::runtime_error("tauSFIdVsJet (2024) evaluate failed for up/down (pt=" + std::to_string(pt[i]) +
                                                         ", dm=" + std::to_string(dm_val) + ", gm=" + std::to_string(gm) + "): " + e.what());
                            }
                            // Inclusive up/down (indices 1, 2) = native 2024 POG total variation.
                            uncSources[1] = up_sf;
                            uncSources[2] = dn_sf;
                            // [TAU-01] option A-1: route into the DM-specific stat1_dm<DM> slot.
                            const int slot = stat1SlotForDM(dm_val);
                            if (slot < 0) {
                                std::cerr << "[ERROR] tauSFIdVsJet (2024): unsupported decay mode dm=" << dm_val
                                          << " (expected 0, 1, 10, 11)" << std::endl;
                                throw std::runtime_error("tauSFIdVsJet (2024): unsupported decay mode dm=" + std::to_string(dm_val));
                            }
                            uncSources[slot]     = up_sf;
                            uncSources[slot + 1] = dn_sf;
                        } else {
                        for (size_t s = 1; s < nSystVsJet; ++s) {
                            try {
                                uncSources[s] = _tauidSFjet->evaluate({pt[i], dm_val, gm, active_tauid_vsjet, tauid_vse, systNamesVsJet[s], "dm"});
                            } catch (const std::exception &e) {
                                std::cerr << "[ERROR] tauSFIdVsJet evaluate failed for systematic '" << systNamesVsJet[s]
                                          << "' (pt=" << pt[i] << ", dm=" << dm_val << ", gm=" << gm << "): " << e.what() << std::endl;
                                throw std::runtime_error("tauSFIdVsJet evaluate failed for systematic '" + systNamesVsJet[s] +
                                                         "' (pt=" + std::to_string(pt[i]) + ", dm=" + std::to_string(dm_val) +
                                                         ", gm=" + std::to_string(gm) + "): " + e.what());
                            }
                        }

                        // Inclusive up/down harmonization:
                        // In 2024, systNamesVsJet[1] ("up") and [2] ("down") are natively provided by POG.
                        // In 2022/2023, POG only provides 32 DM-binned sub-sources (indices 3..38), so evaluate("up")
                        // defaults to nom_sf. If [1] and [2] equal nom_sf, synthesize inclusive up/down via quadrature sum
                        // of the sub-source deviations so that Y<year>_tauidjet is physically consistent across all eras.
                        if (std::abs(uncSources[1] - nom_sf) < 1e-6f && std::abs(uncSources[2] - nom_sf) < 1e-6f) {
                            float sum_sq_up = 0.0f;
                            float sum_sq_dn = 0.0f;
                            for (size_t s = 3; s < nSystVsJet; ++s) {
                                float diff = uncSources[s] - nom_sf;
                                if (diff > 0.0f) sum_sq_up += diff * diff;
                                else             sum_sq_dn += diff * diff;
                            }
                            uncSources[1] = nom_sf + std::sqrt(sum_sq_up);
                            uncSources[2] = nom_sf - std::sqrt(sum_sq_dn);
                        }
                        } // end 2022/2023 schema
                    } else {
                        // SF = 1.0 for non-real-tau objects (jets, electrons, muons faking tau)
                        uncSources.assign(nSystVsJet, 1.0f);
                    }
                    wVec.emplace_back(uncSources);
                    uncSources.clear();
                }
            } else {
                wVec.emplace_back(floats(nSystVsJet, 1.0f));
            }
            return wVec;
        };



        auto tauSFIdVsEl = [this, _tauidSFele, tauid_vse](floats &eta, uchars &dm, uchars &genmatch)->floatsVec {
            // CMS tau POG: DeepTauVSe SF applies to electrons faking tau:
            //   genmatch=1 (prompt e), genmatch=3 (tau→e).
            // Other genmatch values are not in the JSON → exception.
            floats uncSources;
            uncSources.reserve(3);
            floatsVec wVec;
            size_t n = std::min({eta.size(), dm.size(), genmatch.size()});
            wVec.reserve(n > 0 ? n : 1);

            if (n > 0) {
                for (size_t i=0; i<n; i++) {
                    int gm = int(genmatch[i]);
                    if (gm == 1 || gm == 3) {  // prompt e or tau→e
                        int dm_val = int(dm[i]);
                        // [TAU-03] CMS Tau POG DeepTau2018v2p5 only supports HPS decay modes 0, 1, 10, 11.
                        if (dm_val != 0 && dm_val != 1 && dm_val != 10 && dm_val != 11) {
                            std::cerr << "[ERROR] tauSFIdVsEl: unsupported decay mode dm=" << dm_val
                                      << " for genmatch=" << gm << " (expected 0, 1, 10, 11; eta=" << eta[i] << ")" << std::endl;
                            throw std::runtime_error("tauSFIdVsEl: unsupported decay mode dm=" + std::to_string(dm_val));
                        }
                        try {
                            uncSources.emplace_back(_tauidSFele->evaluate({eta[i], dm_val, gm, tauid_vse, "nom"}));
                            uncSources.emplace_back(_tauidSFele->evaluate({eta[i], dm_val, gm, tauid_vse, "up"}));
                            uncSources.emplace_back(_tauidSFele->evaluate({eta[i], dm_val, gm, tauid_vse, "down"}));
                        } catch (const std::exception &e) {
                            std::cerr << "[ERROR] tauSFIdVsEl evaluate failed (eta=" << eta[i] << " dm=" << dm_val
                                      << " gm=" << gm << " vse=" << tauid_vse << "): " << e.what() << std::endl;
                            throw std::runtime_error("tauSFIdVsEl evaluate failed (eta=" + std::to_string(eta[i]) +
                                                     " dm=" + std::to_string(dm_val) + " gm=" + std::to_string(gm) + "): " + e.what());
                        }
                    } else {
                        uncSources.assign(3, 1.0f);
                    }
                    wVec.emplace_back(uncSources);
                    uncSources.clear();
                }
            } else {
                wVec.emplace_back(floats(3, 1.0f));
            }
            return wVec;
        };


        // 2024 DeepTau2018v2p5VSmu schema adds wp_VSe and wp_VSjet inputs:
        // pre-2024: [eta, genmatch, wp, syst]                      → 4 inputs
        // 2024:     [eta, genmatch, wp, wp_VSe, wp_VSjet, syst]    → 6 inputs
        const bool isVsMuSchema6 = (_tauidSFmu->inputs().size() == 6);
        std::string tauid_vsjet_for_vsmu = (tauid_vsjet == "VTight" || tauid_vsjet == "VVTight") ? "Tight" : tauid_vsjet;
        // [TAU-02] Loose-tau (FF, VSjet >= Loose) wp_VSjet for the 2024 6-input VSmu schema.
        // 2024 JSON: wp(VSmu)=VLoose/Tight -> wp_VSjet in {Loose, Medium, Tight}; wp(VSmu)=Loose/Medium -> wp_VSjet = Medium only.
        std::string tauid_vsjet_for_vsmu_loose = (tauid_vsmu == "Loose" || tauid_vsmu == "Medium") ? "Medium" : "Loose";
        if (isVsMuSchema6) {
            std::cout << "[INFO] DeepTau2018v2p5VSmu (6-input schema): wp_VSjet tight='" << tauid_vsjet_for_vsmu
                      << "', loose='" << tauid_vsjet_for_vsmu_loose << "'" << std::endl;
        }
        auto makeTauSFIdVsMu = [this, _tauidSFmu, tauid_vsmu, tauid_vse, isVsMuSchema6](const std::string tauid_vsjet_for_vsmu) {
        return [this, _tauidSFmu, tauid_vsmu, tauid_vse, tauid_vsjet_for_vsmu, isVsMuSchema6](floats &eta, uchars &genmatch)->floatsVec {
            // CMS tau POG: DeepTauVSmu SF applies to muons faking tau:
            //   genmatch=2 (prompt μ), genmatch=4 (tau→μ).
            // Other genmatch values are not in the JSON → exception.
            floats uncSources;
            uncSources.reserve(3);
            floatsVec wVec;
            size_t n = std::min(eta.size(), genmatch.size());
            wVec.reserve(n > 0 ? n : 1);

            if (n > 0) {
                for (size_t i=0; i<n; i++) {
                    int gm = int(genmatch[i]);
                    if (gm == 2 || gm == 4) {  // prompt μ or tau→μ
                        try {
                            if (isVsMuSchema6) {
                                uncSources.emplace_back(_tauidSFmu->evaluate({eta[i], gm, tauid_vsmu, tauid_vse, tauid_vsjet_for_vsmu, "nom"}));
                                uncSources.emplace_back(_tauidSFmu->evaluate({eta[i], gm, tauid_vsmu, tauid_vse, tauid_vsjet_for_vsmu, "up"}));
                                uncSources.emplace_back(_tauidSFmu->evaluate({eta[i], gm, tauid_vsmu, tauid_vse, tauid_vsjet_for_vsmu, "down"}));
                            } else {
                                uncSources.emplace_back(_tauidSFmu->evaluate({eta[i], gm, tauid_vsmu, "nom"}));
                                uncSources.emplace_back(_tauidSFmu->evaluate({eta[i], gm, tauid_vsmu, "up"}));
                                uncSources.emplace_back(_tauidSFmu->evaluate({eta[i], gm, tauid_vsmu, "down"}));
                            }
                        } catch (const std::exception &e) {
                            std::cerr << "[ERROR] tauSFIdVsMu evaluate failed (eta=" << eta[i] << " gm=" << gm << " vsmu=" << tauid_vsmu
                                      << (isVsMuSchema6 ? (" vse=" + tauid_vse + " vsjet=" + tauid_vsjet_for_vsmu) : std::string(""))
                                      << "): " << e.what() << std::endl;
                            throw std::runtime_error("tauSFIdVsMu evaluate failed (eta=" + std::to_string(eta[i]) +
                                                     " gm=" + std::to_string(gm) + "): " + e.what());
                        }
                    } else {
                        uncSources.assign(3, 1.0f);
                    }
                    wVec.emplace_back(uncSources);
                    uncSources.clear();
                }
            } else {
                wVec.emplace_back(floats(3, 1.0f));
            }
            return wVec;
        };
        };
        auto tauSFIdVsMu       = makeTauSFIdVsMu(tauid_vsjet_for_vsmu);
        auto tauSFIdVsMu_loose = makeTauSFIdVsMu(tauid_vsjet_for_vsmu_loose);


        // ---------------------------------------------------------------
        // Tau_genPartFlav is still raw-size here (before any seltaucuts Redefine).
        // First capture Loose-filtered genPartFlav while it is still available at raw size,
        // then Redefine for tight, then compute both sets of SFs from correctly aligned columns.
        // ---------------------------------------------------------------

        // --- Loose genPartFlav: filter raw column with loose mask (same raw size) ---
        _rlm = _rlm.Define("Tau_genPartFlav_loose", "Tau_genPartFlav[seltaucuts_loose]");

        // --- Tight genPartFlav: filter raw column with tight mask ---
        _rlm = _rlm.Redefine("Tau_genPartFlav", "Tau_genPartFlav[seltaucuts]");

        // --- Tight tau SFs: all columns now aligned to Tight selection ---
        _rlm = _rlm.Define("tauWeightIdVsJet", tauSFIdVsJet, {"Tau_pt","Tau_decayMode","Tau_genPartFlav"})
                   .Define("tauWeightIdVsEl",  tauSFIdVsEl,  {"Tau_eta","Tau_decayMode","Tau_genPartFlav"})
                   .Define("tauWeightIdVsMu",  tauSFIdVsMu,  {"Tau_eta","Tau_genPartFlav"});

        _rlm = _rlm.Define("taugencut","Tau_genPartFlav == 5")
                   .Redefine("Tau_pt_gen", "Tau_pt[taugencut]");

        // --- Loose tau SFs: all columns aligned to Loose selection ---
        // A NEW lambda is needed because tauSFIdVsJet captures active_tauid_vsjet by value at
        // definition time. Reassigning the local variable has no effect on it.
        const std::string tauid_vsjet_loose = "Loose";
        std::string active_tauid_vsjet_loose = tauid_vsjet_loose;
        try {
            _tauidSFjet->evaluate({40.0, 0, 5, active_tauid_vsjet_loose, tauid_vse, "nom", "dm"});
            std::cout << "[INFO] DeepTau2018v2p5VSjet (loose): successfully loaded requested WP '" << active_tauid_vsjet_loose << "'" << std::endl;
        } catch (const std::exception& e) {
            if (tauYear != "2024") {
                std::cout << "[INFO] DeepTau2018v2p5VSjet (loose): requested WP '" << active_tauid_vsjet_loose
                          << "' not provided in central " << tauFolder << " release (only 'Medium' is available). "
                          << "Using 'Medium' SF per Tau POG recommendation." << std::endl;
                active_tauid_vsjet_loose = "Medium";
                try {
                    _tauidSFjet->evaluate({40.0, 0, 5, active_tauid_vsjet_loose, tauid_vse, "nom", "dm"});
                } catch (const std::exception& e2) {
                    std::cerr << "[ERROR] DeepTau2018v2p5VSjet (loose): fallback 'Medium' WP also failed in " << tauFolder << ": " << e2.what() << std::endl;
                    throw std::runtime_error("DeepTau2018v2p5VSjet (loose) failed for fallback Medium WP in " + tauFolder + ": " + e2.what());
                }
            } else {
                std::cerr << "[ERROR] DeepTau2018v2p5VSjet (loose): requested WP '" << active_tauid_vsjet_loose
                          << "' failed in " << tauFolder << ": " << e.what() << std::endl;
                throw std::runtime_error("DeepTau2018v2p5VSjet (loose) failed for requested WP " + active_tauid_vsjet_loose + " in " + tauFolder + ": " + e.what());
            }
        }

        auto tauSFIdVsJet_loose = [this, _tauidSFjet, active_tauid_vsjet_loose, tauid_vse, systNamesVsJet, nSystVsJet, isVsJetSchema2024, stat1SlotForDM](floats &pt, uchars &dm, uchars &genmatch)->floatsVec {
            floats uncSources;
            uncSources.reserve(nSystVsJet);
            floatsVec wVec;
            size_t n = std::min({pt.size(), dm.size(), genmatch.size()});
            wVec.reserve(n > 0 ? n : 1);
            if (n > 0) {
                for (size_t i=0; i<n; i++) {
                    int gm = int(genmatch[i]);
                    if (gm == 5) {
                        int dm_val = int(dm[i]);
                        // [TAU-03] CMS Tau POG DeepTau2018v2p5 only supports HPS decay modes 0, 1, 10, 11.
                        if (dm_val != 0 && dm_val != 1 && dm_val != 10 && dm_val != 11) {
                            std::cerr << "[ERROR] tauSFIdVsJetLoose: unsupported decay mode dm=" << dm_val
                                      << " for genmatch=5 (expected 0, 1, 10, 11; pt=" << pt[i] << ")" << std::endl;
                            throw std::runtime_error("tauSFIdVsJetLoose: unsupported decay mode dm=" + std::to_string(dm_val));
                        }
                        float nom_sf = 1.0f;
                        try {
                            nom_sf = _tauidSFjet->evaluate({pt[i], dm_val, gm, active_tauid_vsjet_loose, tauid_vse, "nom", "dm"});
                        } catch (const std::exception &e) {
                            std::cerr << "[ERROR] tauSFIdVsJetLoose evaluate failed (pt=" << pt[i] << ", dm=" << dm_val
                                      << ", gm=" << gm << ", vsjet=" << active_tauid_vsjet_loose << ", vse=" << tauid_vse << "): " << e.what() << std::endl;
                            throw std::runtime_error("tauSFIdVsJetLoose evaluate failed (pt=" + std::to_string(pt[i]) +
                                                     ", dm=" + std::to_string(dm_val) + ", gm=" + std::to_string(gm) +
                                                     ", vsjet=" + active_tauid_vsjet_loose + ", vse=" + tauid_vse + "): " + e.what());
                        }

                        uncSources.assign(nSystVsJet, nom_sf);
                        uncSources[0] = nom_sf;

                        if (isVsJetSchema2024) {
                            // [TAU-02] 2024: JSON provides only 'nom'/'up'/'down' per (DM, pT bin).
                            float up_sf = nom_sf, dn_sf = nom_sf;
                            try {
                                up_sf = _tauidSFjet->evaluate({pt[i], dm_val, gm, active_tauid_vsjet_loose, tauid_vse, "up",   "dm"});
                                dn_sf = _tauidSFjet->evaluate({pt[i], dm_val, gm, active_tauid_vsjet_loose, tauid_vse, "down", "dm"});
                            } catch (const std::exception &e) {
                                std::cerr << "[ERROR] tauSFIdVsJetLoose (2024) evaluate failed for up/down (pt=" << pt[i]
                                          << ", dm=" << dm_val << ", gm=" << gm << "): " << e.what() << std::endl;
                                throw std::runtime_error("tauSFIdVsJetLoose (2024) evaluate failed for up/down (pt=" + std::to_string(pt[i]) +
                                                         ", dm=" + std::to_string(dm_val) + ", gm=" + std::to_string(gm) + "): " + e.what());
                            }
                            uncSources[1] = up_sf;
                            uncSources[2] = dn_sf;
                            // [TAU-01] option A-1: route into the DM-specific stat1_dm<DM> slot.
                            const int slot = stat1SlotForDM(dm_val);
                            if (slot < 0) {
                                std::cerr << "[ERROR] tauSFIdVsJetLoose (2024): unsupported decay mode dm=" << dm_val
                                          << " (expected 0, 1, 10, 11)" << std::endl;
                                throw std::runtime_error("tauSFIdVsJetLoose (2024): unsupported decay mode dm=" + std::to_string(dm_val));
                            }
                            uncSources[slot]     = up_sf;
                            uncSources[slot + 1] = dn_sf;
                        } else {
                        for (size_t s = 1; s < nSystVsJet; ++s) {
                            try {
                                uncSources[s] = _tauidSFjet->evaluate({pt[i], dm_val, gm, active_tauid_vsjet_loose, tauid_vse, systNamesVsJet[s], "dm"});
                            } catch (const std::exception &e) {
                                std::cerr << "[ERROR] tauSFIdVsJetLoose evaluate failed for systematic '" << systNamesVsJet[s]
                                          << "' (pt=" << pt[i] << ", dm=" << dm_val << ", gm=" << gm << "): " << e.what() << std::endl;
                                throw std::runtime_error("tauSFIdVsJetLoose evaluate failed for systematic '" + systNamesVsJet[s] +
                                                         "' (pt=" + std::to_string(pt[i]) + ", dm=" + std::to_string(dm_val) +
                                                         ", gm=" + std::to_string(gm) + "): " + e.what());
                            }
                        }

                        // Inclusive up/down harmonization:
                        if (std::abs(uncSources[1] - nom_sf) < 1e-6f && std::abs(uncSources[2] - nom_sf) < 1e-6f) {
                            float sum_sq_up = 0.0f;
                            float sum_sq_dn = 0.0f;
                            for (size_t s = 3; s < nSystVsJet; ++s) {
                                float diff = uncSources[s] - nom_sf;
                                if (diff > 0.0f) sum_sq_up += diff * diff;
                                else             sum_sq_dn += diff * diff;
                            }
                            uncSources[1] = nom_sf + std::sqrt(sum_sq_up);
                            uncSources[2] = nom_sf - std::sqrt(sum_sq_dn);
                        }
                        } // end 2022/2023 schema
                    } else {
                        uncSources.assign(nSystVsJet, 1.0f);
                    }
                    wVec.emplace_back(uncSources);
                    uncSources.clear();
                }
            } else {
                wVec.emplace_back(floats(nSystVsJet, 1.0f));
            }
            return wVec;
        };

        _rlm = _rlm.Define("tauWeightIdVsJet_loose", tauSFIdVsJet_loose,
                            {"Tau_pt_loose","Tau_decayMode_loose","Tau_genPartFlav_loose"})
                   .Define("tauWeightIdVsEl_loose",  tauSFIdVsEl,
                            {"Tau_eta_loose","Tau_decayMode_loose","Tau_genPartFlav_loose"})
                   .Define("tauWeightIdVsMu_loose",  tauSFIdVsMu_loose,
                            {"Tau_eta_loose","Tau_genPartFlav_loose"});

        _rlm = _rlm.Define("taugencut_loose","Tau_genPartFlav_loose == 5")
                   .Redefine("Tau_pt_loose_gen","Tau_pt_loose[taugencut_loose]");
    }
}

void NanoAODAnalyzerrdframe::selectJets(std::string cut) {
    //https://twiki.cern.ch/twiki/bin/view/CMS/JetID13p6TeV
    //nanoAOD Flags

    // -------------------------------------------------------------------------
    // JES Systematic 4-Vector Variation & Type-1 MET Propagation (MC only)
    // Applied before jetcuts (Jet_pt > 40 GeV) so jets in the skim buffer zone
    // (30 < Jet_pt < 40 GeV) can properly migrate into/out of the selection.
    // -------------------------------------------------------------------------
    if (!_isData && _syst.find("jes") != std::string::npos) {
        std::string jecFile   = "";
        std::string jecPrefix = "";
        std::string eraTag    = "";
        if (_isRun22) {
            jecFile   = "2022_Summer22";
            jecPrefix = "Summer22_22Sep2023_V4_MC_";
            eraTag    = "2022";
        } else if (_isRun22EE) {
            jecFile   = "2022_Summer22EE";
            jecPrefix = "Summer22EE_22Sep2023_V4_MC_";
            eraTag    = "2022EE";
        } else if (_isRun23) {
            jecFile   = "2023_Summer23";
            jecPrefix = "Summer23Prompt23_V4_MC_";
            eraTag    = "2023";
        } else if (_isRun23BPix) {
            jecFile   = "2023_Summer23BPix";
            jecPrefix = "Summer23BPixPrompt23_V4_MC_";
            eraTag    = "2023BPix";
        } else if (_isRun24) {
            jecFile   = "2024_Summer24";
            jecPrefix = "Summer24Prompt24_V5_MC_";
            eraTag    = "2024";
        }

        std::string systClean = _syst;
        while (systClean.rfind("__", 0) == 0) {
            systClean = systClean.substr(2);
        }

        float jesDir = 0.0f;
        std::string jesBase = systClean;
        if (jesBase.size() >= 2 && jesBase.substr(jesBase.size() - 2) == "up") {
            jesDir  = +1.0f;
            jesBase = jesBase.substr(0, jesBase.size() - 2);
        } else if (jesBase.size() >= 4 && jesBase.substr(jesBase.size() - 4) == "down") {
            jesDir  = -1.0f;
            jesBase = jesBase.substr(0, jesBase.size() - 4);
        }

        std::string srcName = jesBase;
        if (srcName.rfind("jes", 0) == 0) {
            srcName = srcName.substr(3);
        }
        for (const char* vtag : {"_v12_", "_v15_"}) {
            size_t vpos = srcName.find(vtag);
            if (vpos != std::string::npos) {
                srcName.erase(vpos, std::string(vtag).size());
            }
        }
        if (srcName.empty()) {
            srcName = "Total";
        }

        // flavMode: 0 = all jets (Regrouped_*), 1 = Gluon (21 or 0),
        //           2 = Light Quark (1, 2, 3), 3 = Charm (4), 4 = Bottom (5)
        int flavMode = 0;
        std::string corrKey = "";
        if (srcName.rfind("FlavorPure", 0) == 0) {
            corrKey = jecPrefix + srcName + "_AK4PFPuppi";
            if      (srcName == "FlavorPureGluon")  flavMode = 1;
            else if (srcName == "FlavorPureQuark")  flavMode = 2;
            else if (srcName == "FlavorPureCharm")  flavMode = 3;
            else if (srcName == "FlavorPureBottom") flavMode = 4;
        } else {
            size_t yearPos = srcName.find("_20");
            if (yearPos != std::string::npos) {
                std::string comp = srcName.substr(0, yearPos);
                corrKey = jecPrefix + "Regrouped_" + comp + "_" + eraTag + "_AK4PFPuppi";
            } else {
                corrKey = jecPrefix + "Regrouped_" + srcName + "_AK4PFPuppi";
            }
        }

        if (!jecFile.empty() && jesDir != 0.0f) {
            std::string jercPath = "data/JME/" + jecFile + "/jet_jerc.json.gz";
            std::cout << "[JES] Applying " << _syst << " -> key=" << corrKey
                      << " from " << jercPath << " (dir=" << (jesDir > 0.0f ? "+1" : "-1") << ")" << std::endl;
            auto jercReader = loadCorrectionSet(jercPath);
            correction::Correction::Ref _jesUnc;
            try {
                _jesUnc = jercReader->at(corrKey);
            } catch (const std::exception &e) {
                std::cerr << "[ERROR] JES correction key '" << corrKey << "' not found in " << jercPath << ": " << e.what() << std::endl;
                throw std::runtime_error("[ERROR] JES correction key '" + corrKey + "' not found in " + jercPath + ": " + e.what());
            }

            if (flavMode == 0) {
                auto calcJesScale = [_jesUnc, jesDir, corrKey](const floats &pts, const floats &etas) -> floats {
                    floats scales;
                    scales.reserve(pts.size());
                    for (size_t i = 0; i < pts.size(); i++) {
                        float unc = 0.0f;
                        try {
                            unc = _jesUnc->evaluate({etas[i], pts[i]});
                        } catch (const std::exception &e) {
                            std::cerr << "[ERROR] JES unc evaluate failed for key '" << corrKey
                                      << "' (eta=" << etas[i] << ", pt=" << pts[i] << "): " << e.what() << std::endl;
                            throw;
                        }
                        float scale = std::max(0.0f, 1.0f + jesDir * unc);
                        scales.emplace_back(scale);
                    }
                    return scales;
                };
                _rlm = _rlm.Define("Jet_jes_scale", calcJesScale, {"Jet_pt", "Jet_eta"});
            } else {
                // Cast flavor/index columns to int (RVec<int>) for uniform v9/v12/v15 compatibility
                const std::string partFlavExpr = isDefined("Jet_partonFlavour")
                    ? "ROOT::VecOps::Cast<int>(Jet_partonFlavour)" : "ROOT::VecOps::RVec<int>{}";
                const std::string genIdxExpr = isDefined("Jet_genJetIdx")
                    ? "ROOT::VecOps::Cast<int>(Jet_genJetIdx)" : "ROOT::VecOps::RVec<int>{}";
                const std::string genPartFlavExpr = isDefined("GenJet_partonFlavour")
                    ? "ROOT::VecOps::Cast<int>(GenJet_partonFlavour)" : "ROOT::VecOps::RVec<int>{}";
                const std::string hadFlavExpr = isDefined("Jet_hadronFlavour")
                    ? "ROOT::VecOps::Cast<int>(Jet_hadronFlavour)" : "ROOT::VecOps::RVec<int>{}";

                _rlm = _rlm.Define("Jet_partonFlavour_jes_i",    partFlavExpr)
                           .Define("Jet_genJetIdx_jes_i",        genIdxExpr)
                           .Define("GenJet_partonFlavour_jes_i", genPartFlavExpr)
                           .Define("Jet_hadronFlavour_jes_i",    hadFlavExpr);

                auto calcJesScaleFlav = [_jesUnc, jesDir, flavMode, corrKey](const floats &pts, const floats &etas,
                                                                    const ints &partflavs, const ints &genIdxs,
                                                                    const ints &genPartflavs, const ints &hadflavs) -> floats {
                    floats scales;
                    scales.reserve(pts.size());
                    const bool alignedPartFlav = (partflavs.size() == pts.size());
                    for (size_t i = 0; i < pts.size(); i++) {
                        int pf = 0;
                        if (alignedPartFlav) {
                            pf = std::abs(partflavs[i]);
                        } else if (i < genIdxs.size() && genIdxs[i] >= 0 && static_cast<size_t>(genIdxs[i]) < genPartflavs.size()) {
                            pf = std::abs(genPartflavs[genIdxs[i]]);
                        } else if (i < hadflavs.size() && (std::abs(hadflavs[i]) == 5 || std::abs(hadflavs[i]) == 4)) {
                            pf = std::abs(hadflavs[i]);
                        }

                        bool isTargetFlav = false;
                        // Order matches Run 2 JME POG definition: Gluon / Quark / Charm / Bottom
                        if      (flavMode == 1 && (pf == 21 || pf == 0))           isTargetFlav = true;
                        else if (flavMode == 2 && (pf == 1 || pf == 2 || pf == 3)) isTargetFlav = true;
                        else if (flavMode == 3 && (pf == 4))                       isTargetFlav = true;
                        else if (flavMode == 4 && (pf == 5))                       isTargetFlav = true;

                        if (isTargetFlav) {
                            float unc = 0.0f;
                            try {
                                unc = _jesUnc->evaluate({etas[i], pts[i]});
                            } catch (const std::exception &e) {
                                std::cerr << "[ERROR] JES flav unc evaluate failed for key '" << corrKey
                                          << "' (eta=" << etas[i] << ", pt=" << pts[i] << "): " << e.what() << std::endl;
                                throw;
                            }
                            float scale = std::max(0.0f, 1.0f + jesDir * unc);
                            scales.emplace_back(scale);
                        } else {
                            scales.emplace_back(1.0f);
                        }
                    }
                    return scales;
                };

                _rlm = _rlm.Define("Jet_jes_scale", calcJesScaleFlav,
                                   {"Jet_pt", "Jet_eta", "Jet_partonFlavour_jes_i",
                                    "Jet_genJetIdx_jes_i", "GenJet_partonFlavour_jes_i", "Jet_hadronFlavour_jes_i"});
            }

            _rlm = _rlm.Define("Jet_pt_jes_shifted",   "Jet_pt * Jet_jes_scale")
                       .Define("Jet_mass_jes_shifted", "Jet_mass * Jet_jes_scale");

            auto jesMetPt = [](const floats &pts_orig, const floats &pts_shifted, const floats &phis, float met_pt, float met_phi) -> float {
                double metx = met_pt * std::cos(met_phi);
                double mety = met_pt * std::sin(met_phi);
                for (size_t i = 0; i < pts_orig.size(); i++) {
                    if (pts_orig[i] > 15.0f || pts_shifted[i] > 15.0f) {
                        double dpt = static_cast<double>(pts_shifted[i]) - static_cast<double>(pts_orig[i]);
                        metx -= dpt * std::cos(phis[i]);
                        mety -= dpt * std::sin(phis[i]);
                    }
                }
                return static_cast<float>(std::sqrt(metx * metx + mety * mety));
            };

            auto jesMetPhi = [](const floats &pts_orig, const floats &pts_shifted, const floats &phis, float met_pt, float met_phi) -> float {
                double metx = met_pt * std::cos(met_phi);
                double mety = met_pt * std::sin(met_phi);
                for (size_t i = 0; i < pts_orig.size(); i++) {
                    if (pts_orig[i] > 15.0f || pts_shifted[i] > 15.0f) {
                        double dpt = static_cast<double>(pts_shifted[i]) - static_cast<double>(pts_orig[i]);
                        metx -= dpt * std::cos(phis[i]);
                        mety -= dpt * std::sin(phis[i]);
                    }
                }
                return static_cast<float>(std::atan2(mety, metx));
            };

            _rlm = _rlm.Define("PuppiMET_pt_jes_uncorr",  "PuppiMET_pt")
                       .Define("PuppiMET_phi_jes_uncorr", "PuppiMET_phi")
                       .Redefine("PuppiMET_phi", jesMetPhi, {"Jet_pt", "Jet_pt_jes_shifted", "Jet_phi", "PuppiMET_pt_jes_uncorr", "PuppiMET_phi_jes_uncorr"})
                       .Redefine("PuppiMET_pt",  jesMetPt,  {"Jet_pt", "Jet_pt_jes_shifted", "Jet_phi", "PuppiMET_pt_jes_uncorr", "PuppiMET_phi_jes_uncorr"});

            if (isDefined("MET_pt") && isDefined("MET_phi")) {
                _rlm = _rlm.Define("MET_pt_jes_uncorr",  "MET_pt")
                           .Define("MET_phi_jes_uncorr", "MET_phi")
                           .Redefine("MET_phi", jesMetPhi, {"Jet_pt", "Jet_pt_jes_shifted", "Jet_phi", "MET_pt_jes_uncorr", "MET_phi_jes_uncorr"})
                           .Redefine("MET_pt",  jesMetPt,  {"Jet_pt", "Jet_pt_jes_shifted", "Jet_phi", "MET_pt_jes_uncorr", "MET_phi_jes_uncorr"});
            }

            _rlm = _rlm.Redefine("Jet_pt",   "Jet_pt_jes_shifted")
                       .Redefine("Jet_mass", "Jet_mass_jes_shifted");
        }
    }

    // -------------------------------------------------------------------------
    // JER Systematic 4-Vector Variation & Type-1 MET Propagation (MC only)
    // Uses relative ratio c_JER(s_var) / c_JER(s_nom) on top of skimmed jets
    // before jetcuts (Jet_pt > 40 GeV) so jets in the skim buffer zone
    // (30 < Jet_pt < 40 GeV) can properly migrate into/out of the selection.
    // -------------------------------------------------------------------------
    if (!_isData && _syst.find("jer") != std::string::npos) {
        std::string jecFile   = "";
        std::string jerPrefix = "";
        if (_isRun22) {
            jecFile   = "2022_Summer22";
            jerPrefix = "Summer22_22Sep2023_JRV2";
        } else if (_isRun22EE) {
            jecFile   = "2022_Summer22EE";
            jerPrefix = "Summer22EE_22Sep2023_JRV2";
        } else if (_isRun23) {
            jecFile   = "2023_Summer23";
            jerPrefix = "Summer23Prompt23_RunCv1234_JRV3";
        } else if (_isRun23BPix) {
            jecFile   = "2023_Summer23BPix";
            jerPrefix = "Summer23BPixPrompt23_RunD_JRV3";
        } else if (_isRun24) {
            jecFile   = "2024_Summer24";
            jerPrefix = "Summer24Prompt24_JRV2";
        }

        float jerDir = 0.0f;
        if (_syst.find("jerup") != std::string::npos) {
            jerDir = +1.0f;
        } else if (_syst.find("jerdown") != std::string::npos) {
            jerDir = -1.0f;
        }

        if (!jecFile.empty() && jerDir != 0.0f) {
            std::string jercPath = "data/JME/" + jecFile + "/jet_jerc.json.gz";
            std::cout << "[JER] Applying " << _syst << " -> prefix=" << jerPrefix
                      << " from " << jercPath << " (dir=" << (jerDir > 0.0f ? "+1" : "-1") << ")" << std::endl;
            auto jercReader = loadCorrectionSet(jercPath);
            auto _jerReso   = jercReader->at(jerPrefix + "_MC_PtResolution_AK4PFPuppi");
            auto _jerSF     = jercReader->at(jerPrefix + "_MC_ScaleFactor_AK4PFPuppi");
            auto _jerUnc    = jercReader->at(jerPrefix + "_MC_SFUncertainty_AK4PFPuppi");
            auto jerReader  = loadCorrectionSet("data/JME/jer_smear.json.gz");
            auto _jerSmear  = jerReader->at("JERSmear");

            _rlm = _rlm.Define("Jet_genJetIdx_jer_i", "ROOT::VecOps::Cast<int>(Jet_genJetIdx)");

            auto calcJerScale = [_jerReso, _jerSF, _jerUnc, _jerSmear, jerDir](
                                    const floats &pts, const floats &etas, const floats &phis,
                                    const ints &genidxs, const floats &genpts,
                                    const floats &genetas, const floats &genphis,
                                    float rho, ULong64_t event) -> floats {
                floats scales;
                scales.reserve(pts.size());
                for (size_t i = 0; i < pts.size(); i++) {
                    float reso   = _jerReso->evaluate({etas[i], pts[i], rho});
                    float sf_nom = _jerSF->evaluate({etas[i], pts[i]});
                    float sf_unc = _jerUnc->evaluate({etas[i], pts[i]});
                    float sf_var = std::max(0.0f, sf_nom + jerDir * sf_unc);

                    float genPtForSmear = -1.0f;
                    int genidx = (i < genidxs.size()) ? genidxs[i] : -1;
                    if (genidx >= 0 && static_cast<size_t>(genidx) < genpts.size()) {
                        float dPhi = std::abs(genphis[genidx] - phis[i]);
                        float dEta = std::abs(genetas[genidx] - etas[i]);
                        float dR2  = dPhi * dPhi + dEta * dEta;
                        if (dR2 < 0.04f && std::abs(pts[i] - genpts[genidx]) < (3.0f * reso * pts[i])) {
                            genPtForSmear = genpts[genidx];
                        }
                    }

                    float smear_nom = _jerSmear->evaluate({pts[i], etas[i], genPtForSmear, rho, int(event), reso, sf_nom});
                    float smear_var = _jerSmear->evaluate({pts[i], etas[i], genPtForSmear, rho, int(event), reso, sf_var});
                    float c_nom = (std::isfinite(smear_nom) && smear_nom > 0.0f) ? smear_nom : 1.0f;
                    float c_var = (std::isfinite(smear_var) && smear_var > 0.0f) ? smear_var : 1.0f;
                    float scale = (c_nom > 0.0f) ? (c_var / c_nom) : 1.0f;
                    scales.emplace_back(scale);
                }
                return scales;
            };

            _rlm = _rlm.Define("Jet_jer_scale", calcJerScale,
                               {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_genJetIdx_jer_i",
                                "GenJet_pt", "GenJet_eta", "GenJet_phi",
                                "Rho_fixedGridRhoFastjetAll", "event"})
                       .Define("Jet_pt_jer_shifted",   "Jet_pt * Jet_jer_scale")
                       .Define("Jet_mass_jer_shifted", "Jet_mass * Jet_jer_scale");

            auto jerMetPt = [](const floats &pts_orig, const floats &pts_shifted, const floats &phis, float met_pt, float met_phi) -> float {
                double metx = met_pt * std::cos(met_phi);
                double mety = met_pt * std::sin(met_phi);
                for (size_t i = 0; i < pts_orig.size(); i++) {
                    if (pts_orig[i] > 15.0f || pts_shifted[i] > 15.0f) {
                        double dpt = static_cast<double>(pts_shifted[i]) - static_cast<double>(pts_orig[i]);
                        metx -= dpt * std::cos(phis[i]);
                        mety -= dpt * std::sin(phis[i]);
                    }
                }
                return static_cast<float>(std::sqrt(metx * metx + mety * mety));
            };

            auto jerMetPhi = [](const floats &pts_orig, const floats &pts_shifted, const floats &phis, float met_pt, float met_phi) -> float {
                double metx = met_pt * std::cos(met_phi);
                double mety = met_pt * std::sin(met_phi);
                for (size_t i = 0; i < pts_orig.size(); i++) {
                    if (pts_orig[i] > 15.0f || pts_shifted[i] > 15.0f) {
                        double dpt = static_cast<double>(pts_shifted[i]) - static_cast<double>(pts_orig[i]);
                        metx -= dpt * std::cos(phis[i]);
                        mety -= dpt * std::sin(phis[i]);
                    }
                }
                return static_cast<float>(std::atan2(mety, metx));
            };

            _rlm = _rlm.Define("PuppiMET_pt_jer_uncorr",  "PuppiMET_pt")
                       .Define("PuppiMET_phi_jer_uncorr", "PuppiMET_phi")
                       .Redefine("PuppiMET_phi", jerMetPhi, {"Jet_pt", "Jet_pt_jer_shifted", "Jet_phi", "PuppiMET_pt_jer_uncorr", "PuppiMET_phi_jer_uncorr"})
                       .Redefine("PuppiMET_pt",  jerMetPt,  {"Jet_pt", "Jet_pt_jer_shifted", "Jet_phi", "PuppiMET_pt_jer_uncorr", "PuppiMET_phi_jer_uncorr"});

            if (isDefined("MET_pt") && isDefined("MET_phi")) {
                _rlm = _rlm.Define("MET_pt_jer_uncorr",  "MET_pt")
                           .Define("MET_phi_jer_uncorr", "MET_phi")
                           .Redefine("MET_phi", jerMetPhi, {"Jet_pt", "Jet_pt_jer_shifted", "Jet_phi", "MET_pt_jer_uncorr", "MET_phi_jer_uncorr"})
                           .Redefine("MET_pt",  jerMetPt,  {"Jet_pt", "Jet_pt_jer_shifted", "Jet_phi", "MET_pt_jer_uncorr", "MET_phi_jer_uncorr"});
            }

            _rlm = _rlm.Redefine("Jet_pt",   "Jet_pt_jer_shifted")
                       .Redefine("Jet_mass", "Jet_mass_jer_shifted");
        }
    }

    _rlm = _rlm.Define("jetcuts", cut);
    _rlm = _rlm.Redefine("Jet_pt", "Jet_pt[jetcuts]")
               .Redefine("Jet_eta", "Jet_eta[jetcuts]")
               .Redefine("Jet_phi", "Jet_phi[jetcuts]")
               .Redefine("Jet_mass", "Jet_mass[jetcuts]")
               .Redefine("Jet_btagPNetB", "Jet_btagPNetB[jetcuts]")
               .Redefine("Jet_btagPNetCvB", "Jet_btagPNetCvB[jetcuts]")
               .Redefine("Jet_btagPNetCvL", "Jet_btagPNetCvL[jetcuts]")
               .Define("jet4vecs", ::gen4vec, {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass"});
    if (_isRun24){
        _rlm = _rlm.Redefine("Jet_btagUParTAK4B", "Jet_btagUParTAK4B[jetcuts]");
        _rlm = _rlm.Redefine("Jet_btagUParTAK4CvB", "Jet_btagUParTAK4CvB[jetcuts]");
        _rlm = _rlm.Redefine("Jet_btagUParTAK4CvL", "Jet_btagUParTAK4CvL[jetcuts]");
    }

    if (!_isData) {
        _rlm = _rlm.Redefine("Jet_hadronFlavour", "Jet_hadronFlavour[jetcuts]");
        //_rlm = _rlm.Redefine("btagWeight_PNetB_perJet", skimCol, {"btagWeight_PNetB_perJet", "jetcuts"})
        //           .Redefine("btagWeight_PNetB_jes_perJet", skimCol, {"btagWeight_PNetB_jes_perJet", "jetcuts"});
    }

    // for checking overlapped jets with leptons
    auto checkoverlap = [](FourVectorVec &seljets, FourVectorVec &sellep) {

        ints mindrlepton;
        for (auto ajet: seljets) {
            auto mindr = 6.0;
            for ( auto alepton : sellep ) {
                auto dr = ROOT::Math::VectorUtil::DeltaR(ajet, alepton);
                if (dr < mindr) mindr = dr;
            }
            int out = mindr >= 0.4 ? 1 : 0;
            mindrlepton.emplace_back(out);
        }
        return mindrlepton;
    };

    // Overlap removal with muon / electron (used for btagging SF)
    if (isDefined("cleantau4vecs")){
        _rlm = _rlm.Define("lepjetoverlap", checkoverlap, {"jet4vecs","lep4vecs"});
        _rlm = _rlm.Define("taujetoverlap", checkoverlap, {"jet4vecs","cleantau4vecs"})
                   .Define("loosetaujetoverlap", checkoverlap, {"jet4vecs","cleanloosetau4vecs"})
                   .Define("jetoverlap","lepjetoverlap && taujetoverlap")
                   .Define("jetoverlaploose","lepjetoverlap && loosetaujetoverlap");

    } else{
        _rlm = _rlm.Define("jetoverlap", checkoverlap, {"jet4vecs","lep4vecs"});
        _rlm = _rlm.Define("jetoverlaploose", checkoverlap, {"jet4vecs","lep4vecs"});
    }
    _rlm = _rlm.Define("Jet_pt_loose", "Jet_pt[jetoverlaploose]")
               .Define("Jet_eta_loose", "Jet_eta[jetoverlaploose]")
               .Define("Jet_btagPNetB_loose", "Jet_btagPNetB[jetoverlaploose]")
               .Define("Jet_btagPNetCvB_loose", "Jet_btagPNetCvB[jetoverlaploose]")
               .Define("Jet_btagPNetCvL_loose", "Jet_btagPNetCvL[jetoverlaploose]")
               .Define("ncleanjetsloosepass", "int(Jet_pt_loose.size())");
    if (_isRun24){
        _rlm = _rlm.Define("Jet_btagUParTAK4B_loose", "Jet_btagUParTAK4B[jetoverlaploose]");
        _rlm = _rlm.Define("Jet_btagUParTAK4CvB_loose", "Jet_btagUParTAK4CvB[jetoverlaploose]");
        _rlm = _rlm.Define("Jet_btagUParTAK4CvL_loose", "Jet_btagUParTAK4CvL[jetoverlaploose]");
    }

    _rlm = _rlm.Redefine("Jet_pt", "Jet_pt[jetoverlap]")
               .Redefine("Jet_eta", "Jet_eta[jetoverlap]")
               .Redefine("Jet_phi", "Jet_phi[jetoverlap]")
               .Redefine("Jet_mass", "Jet_mass[jetoverlap]")
               .Redefine("Jet_btagPNetB", "Jet_btagPNetB[jetoverlap]")
               .Redefine("Jet_btagPNetCvB", "Jet_btagPNetCvB[jetoverlap]")
               .Redefine("Jet_btagPNetCvL", "Jet_btagPNetCvL[jetoverlap]")
               .Define("ncleanjetspass", "int(Jet_pt.size())")
               .Define("cleanjet4vecs", ::gen4vec, {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass"})
               .Define("Jet_HT", "Sum(Jet_pt)");
    if (_isRun24){
        _rlm = _rlm.Redefine("Jet_btagUParTAK4B", "Jet_btagUParTAK4B[jetoverlap]");
        _rlm = _rlm.Redefine("Jet_btagUParTAK4CvB", "Jet_btagUParTAK4CvB[jetoverlap]");
        _rlm = _rlm.Redefine("Jet_btagUParTAK4CvL", "Jet_btagUParTAK4CvL[jetoverlap]");
    }

    if (!_isData) {
        _rlm = _rlm.Define("Jet_hadronFlavour_loose", "Jet_hadronFlavour[jetoverlaploose]")
                   .Redefine("Jet_hadronFlavour", "Jet_hadronFlavour[jetoverlap]");
        //int nbsf_var = btag_var.size();
        //int njes_var = jes_var.size();
        //_rlm = _rlm.Define("nbsf_var", [nbsf_var](){return int(nbsf_var);})
        //           .Define("njes_var", [njes_var](){return int(njes_var);})
                   //.Define("btagWeight_PNetB_perJet_loose", skimCol, {"btagWeight_PNetB_perJet", "jetoverlaploose"})
                   //.Define("btagWeight_PNetB_loose", calcBSF, {"btagWeight_PNetB_perJet_loose", "nbsf_var"})
                   //.Redefine("btagWeight_PNetB_perJet", skimCol, {"btagWeight_PNetB_perJet", "jetoverlap"})
                   //.Redefine("btagWeight_PNetB_jes_perJet", skimCol, {"btagWeight_PNetB_jes_perJet", "jetoverlap"})
                   //.Define("btagWeight", calcBSF, {"btagWeight", "nbsf_var"});
                   //.Define("btagWeight_PNetB_jes", calcBSF, {"btagWeight_PNetB_jes_perJet", "njes_var"});
    }


    // b-tagging (Medium WP for both tight-tau-cleaned jets and loose-tau-cleaned jets)
    // Threshold _btagcut is centrally loaded from data/config/analysis_config.json.
    if (_btagcut <= 0.0f || !std::isfinite(_btagcut)) {
        std::string msg = "[ERROR] selectJets: Invalid or uninitialized _btagcut (" + std::to_string(_btagcut) +
                          ") for year=" + _year + ". Ensure _btagcut is loaded from analysis_config.json before calling selectJets.";
        std::cerr << msg << std::endl;
        throw std::runtime_error(msg);
    }
    const std::string btagCutStr = std::to_string(_btagcut);
    if (_isRun22 || _isRun22EE || _isRun23 || _isRun23BPix) {
        _rlm = _rlm.Define("btagcuts", "Jet_btagPNetB > " + btagCutStr)
                   .Define("btagcuts_loose", "Jet_btagPNetB_loose > " + btagCutStr);
    } else if (_isRun24) {
        _rlm = _rlm.Define("btagcuts", "Jet_btagUParTAK4B > " + btagCutStr)
                   .Define("btagcuts_loose", "Jet_btagUParTAK4B_loose > " + btagCutStr);
    } else {
        std::string msg = "[ERROR] selectJets: Unrecognized Run 3 era flag for year=" + _year;
        std::cerr << msg << std::endl;
        throw std::runtime_error(msg);
    }

    _rlm = _rlm.Define("bJet_pt", "Jet_pt[btagcuts]")
               .Define("bJet_eta", "Jet_eta[btagcuts]")
               .Define("bJet_phi", "Jet_phi[btagcuts]")
               .Define("bJet_mass", "Jet_mass[btagcuts]")
               .Define("bJet_btagPNetB", "Jet_btagPNetB[btagcuts]")
               .Define("bJet_btagPNetCvB", "Jet_btagPNetCvB[btagcuts]")
               .Define("bJet_btagPNetCvL", "Jet_btagPNetCvL[btagcuts]")
               .Define("ncleanbjetspass", "int(bJet_pt.size())")
               .Define("bJet_HT", "Sum(bJet_pt)")
               .Define("cleanbjet4vecs", ::gen4vec, {"bJet_pt", "bJet_eta", "bJet_phi", "bJet_mass"})
               .Define("bJet_pt_loose", "Jet_pt_loose[btagcuts_loose]")
               .Define("ncleanbjetsloosepass", "int(bJet_pt_loose.size())");
    if (_isRun24){
        _rlm = _rlm.Define("bJet_btagUParTAK4B", "Jet_btagUParTAK4B[btagcuts]");
        _rlm = _rlm.Define("bJet_btagUParTAK4CvB", "Jet_btagUParTAK4CvB[btagcuts]");
        _rlm = _rlm.Define("bJet_btagUParTAK4CvL", "Jet_btagUParTAK4CvL[btagcuts]");
    }
}

void NanoAODAnalyzerrdframe::selectFatJets() {
    _rlm = _rlm.Define("fatjetcuts", "FatJet_pt>400.0 && abs(FatJet_eta)<2.4 && FatJet_tau1>0.0 && FatJet_tau2>0.0 && FatJet_tau3>0.0 && FatJet_tau3/FatJet_tau2<0.5")
               .Define("Sel_fatjetpt", "FatJet_pt[fatjetcuts]")
               .Define("Sel_fatjeteta", "FatJet_eta[fatjetcuts]")
               .Define("Sel_fatjetphi", "FatJet_phi[fatjetcuts]")
               .Define("Sel_fatjetmass", "FatJet_mass[fatjetcuts]")
               .Define("nfatjetspass", "int(Sel_fatjetpt.size())")
               //.Define("Sel_fatjetweight", "std::vector<double>(nfatjetspass, evWeight)")
               .Define("Sel_fatjet4vecs", ::gen4vec, {"Sel_fatjetpt", "Sel_fatjeteta", "Sel_fatjetphi", "Sel_fatjetmass"});
}


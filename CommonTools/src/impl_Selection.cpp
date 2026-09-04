#include "NanoAODAnalyzerrdframe.h"
#include <iostream>
#include <algorithm>
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
#include "TH2F.h"
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
    }
}

void NanoAODAnalyzerrdframe::applyBSFs(std::vector<string> jes_var, string btagYear, string btagMap, string btagMapLight, float btagcut) {
    cout << "Loading Btag SF (fixed WP, Method 1a): " << btagYear << "  comb=" << btagMap << "  light=" << btagMapLight << endl;
    // btag_var / btag_jes_var columns not needed for fixed-WP SF (shape-only).

    auto bSFreader = loadCorrectionSet("data/BTV/" + btagYear + "/btagging.json.gz");
    // Fixed-WP corrections:
    //   btagMap      = e.g. "particleNet_comb" (2022/23) or "UParTAK4_comb" (2024) — for b/c jets
    //   btagMapLight = e.g. "particleNet_light"           or "UParTAK4_light"       — for light jets
    // Inputs: (systematic, working_point, flavor, abseta, pt)  — no discriminant needed.
    auto _btagSF      = bSFreader->at(btagMap);
    auto _btagSFlight = bSFreader->at(btagMapLight);

    // -----------------------------------------------------------------------
    // BTV Method 1a (Event Reweighting) MC Efficiency Map Loading
    // Weight = Prod_{tagged} SF_i * Prod_{untagged} (1 - SF_j * eff_j) / (1 - eff_j)
    // If btag_eff.root does not exist, gracefully fall back to Case 1 (tagged only).
    // -----------------------------------------------------------------------
    std::string effPath = "data/BTV/" + btagYear + "/btag_eff.root";
    std::shared_ptr<TH2F> sp_eff_b = nullptr;
    std::shared_ptr<TH2F> sp_eff_c = nullptr;
    std::shared_ptr<TH2F> sp_eff_light = nullptr;
    bool has_eff = false;

    TFile *fEff = TFile::Open(effPath.c_str(), "READ");
    if (fEff && !fEff->IsZombie()) {
        TH2F *h_b = dynamic_cast<TH2F*>(fEff->Get("h2_eff_b"));
        TH2F *h_c = dynamic_cast<TH2F*>(fEff->Get("h2_eff_c"));
        TH2F *h_l = dynamic_cast<TH2F*>(fEff->Get("h2_eff_light"));
        if (h_b && h_c && h_l) {
            h_b->SetDirectory(0);
            h_c->SetDirectory(0);
            h_l->SetDirectory(0);
            sp_eff_b = std::shared_ptr<TH2F>(h_b);
            sp_eff_c = std::shared_ptr<TH2F>(h_c);
            sp_eff_light = std::shared_ptr<TH2F>(h_l);
            has_eff = true;
            cout << "[INFO] Loaded BTV Method 1a efficiency maps from " << effPath << endl;
        } else {
            cout << "[WARN] Could not retrieve all 3 efficiency histograms from " << effPath
                 << ". Falling back to simple method." << endl;
        }
        fEff->Close();
        delete fEff;
    } else {
        cout << "[INFO] b-tag efficiency file not found at " << effPath
             << ". Falling back to simple method (tagged jets only)." << endl;
    }

    // -----------------------------------------------------------------------
    // Systematic indices (same order for all eras, only common systs used):
    //   0: central
    //   1/2:  up/down_correlated
    //   3/4:  up/down_uncorrelated
    //   5/6:  up/down_statistic
    //   7/8:  up/down_type3
    //   9/10: up/down_bfragmentation
    //
    // Light-jet systs (applied only to hadronFlavour==0):
    //   0: central  1/2: up/down_correlated  3/4: up/down_uncorrelated
    // -----------------------------------------------------------------------
    const std::string wp = "M";  // Medium WP
    const std::vector<std::string> systs_hfbc = {
        "central",
        "up_correlated",   "down_correlated",    // 1/2
        "up_uncorrelated", "down_uncorrelated",  // 3/4
        "up_statistic",    "down_statistic",      // 5/6
        "up_type3",        "down_type3",          // 7/8
        "up_bfragmentation","down_bfragmentation" // 9/10
    };
    const std::vector<std::string> systs_light = {
        "central",
        "up_correlated",   "down_correlated",    // 1/2
        "up_uncorrelated", "down_uncorrelated"   // 3/4
    };
    const int n_systs = static_cast<int>(systs_hfbc.size()); // 11

    // Helper lambda to query efficiency from 2D map with boundary clamping
    auto get_efficiency = [has_eff, sp_eff_b, sp_eff_c, sp_eff_light](int flav, float pt, float abseta) -> float {
        if (!has_eff) return 1.0f;
        TH2F *h = nullptr;
        float def_eff = 0.65f;
        if (flav == 5) {
            h = sp_eff_b.get();
            def_eff = 0.65f;
        } else if (flav == 4) {
            h = sp_eff_c.get();
            def_eff = 0.15f;
        } else {
            h = sp_eff_light.get();
            def_eff = 0.01f;
        }
        if (!h) return def_eff;

        float minPt  = h->GetXaxis()->GetXmin();
        float maxPt  = h->GetXaxis()->GetXmax() - 1e-3f;
        float minEta = h->GetYaxis()->GetXmin();
        float maxEta = h->GetYaxis()->GetXmax() - 1e-3f;

        float clampedPt  = std::clamp(pt, minPt, maxPt);
        float clampedEta = std::clamp(abseta, minEta, maxEta);

        int bx = h->GetXaxis()->FindFixBin(clampedPt);
        int by = h->GetYaxis()->FindFixBin(clampedEta);
        float val = h->GetBinContent(bx, by);
        if (val <= 0.0f || val > 1.0f) val = def_eff;
        return val;
    };

    auto btagSF_fixedWP = [this, _btagSF, _btagSFlight, wp, systs_hfbc, systs_light, n_systs, btagcut, has_eff, get_efficiency]
                          (floats &pts, floats &etas, uchars &hadflav, floats &bscores) -> floats {
        // Returns event-level weight vector [n_systs].
        // Under Method 1a:
        //   Tagged jets: multiply SF
        //   Untagged jets: multiply (1 - SF * eff) / (1 - eff)
        // Under simple fallback (has_eff == false):
        //   Tagged jets: multiply SF
        //   Untagged jets: weight = 1.0
        floats out(n_systs, 1.0f);
        size_t n_jets = std::min({pts.size(), etas.size(), hadflav.size(), bscores.size()});
        for (size_t j = 0; j < n_jets; j++) {
            float pt    = pts[j];
            float abseta = std::abs(etas[j]);
            int   flav  = static_cast<int>(hadflav[j]);
            float score = bscores[j];

            // Skip jets outside kinematic range
            if (pt < 20.f || abseta > 2.5f) continue;
            bool is_tagged = (score >= btagcut);

            // If efficiency map not available, untagged jets do not contribute (SF=1)
            if (!has_eff && !is_tagged) continue;

            float eff = get_efficiency(flav, pt, abseta);

            for (int s = 0; s < n_systs; s++) {
                float sf = 1.0f;
                try {
                    if (flav == 0) {
                        // Light jets: use _light correction; only 5 systs available
                        const auto &syst_light = (s < static_cast<int>(systs_light.size()))
                                                  ? systs_light[s] : systs_light[0];
                        sf = _btagSFlight->evaluate({syst_light, wp, flav, abseta, pt});
                    } else {
                        // b/c jets: use _comb correction
                        sf = _btagSF->evaluate({systs_hfbc[s], wp, flav, abseta, pt});
                    }
                } catch (const std::exception &e) {
                    std::cerr << "[WARN] btagSF_fixedWP evaluate failed"
                              << " syst=" << systs_hfbc[s]
                              << " wp=" << wp << " flav=" << flav
                              << " eta=" << abseta << " pt=" << pt
                              << ": " << e.what() << "\n";
                    sf = 1.0f;
                }

                if (is_tagged) {
                    out[s] *= sf;
                } else {
                    // Method 1a: Untagged jet factor
                    float denom = std::max(1e-4f, 1.0f - eff);
                    float w_fail = (1.0f - sf * eff) / denom;
                    out[s] *= std::clamp(w_fail, 0.0f, 10.0f);
                }
            }
        }
        return out;
    };

    // btagWeight[0] = central, [1/2] = corr up/dn, [3/4] = uncorr up/dn,
    // [5/6] = stat up/dn, [7/8] = type3 up/dn, [9/10] = bfrag up/dn
    const std::string bTagCol = _isRun24 ? "Jet_btagUParTAK4B" : "Jet_btagPNetB";
    _rlm = _rlm.Define("btagWeight", btagSF_fixedWP,
                        {"Jet_pt", "Jet_eta", "Jet_hadronFlavour", bTagCol});
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
    if (!_isData) {
        // In the Process stage, Tau_pt has already been corrected by Tau_ES_nom in the skim.
        // Tau_ES_up and Tau_ES_down are stored as flat float vectors in the skim tree.
        // To apply the up/down variation: rescale by Tau_ES_up/Tau_ES_nom or Tau_ES_down/Tau_ES_nom.
        // Guard against Tau_ES_nom == 0 (unshifted jet) and against old skim files that
        // pre-date Tau_ES_* columns.
        const bool hasTauES = isDefined("Tau_ES_nom");
        if (!hasTauES && _syst.find("tes") != std::string::npos) {
            std::cerr << "[WARN] selectTaus: syst=" << _syst
                      << " but Tau_ES_nom not found in input — skim may be outdated. "
                      << "TES variation will NOT be applied.\n";
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
               .Redefine("Tau_pt", "Tau_pt[seltaucuts]")
               .Redefine("Tau_eta", "Tau_eta[seltaucuts]")
               .Redefine("Tau_phi", "Tau_phi[seltaucuts]")
               .Redefine("Tau_mass", "Tau_mass[seltaucuts]")
               .Redefine("Tau_charge", "Tau_charge[seltaucuts]")
               .Redefine("Tau_jetIdx", "Tau_jetIdx[seltaucuts]")
               .Redefine("Tau_decayMode", "Tau_decayMode[seltaucuts]")
               .Define("ncleantaupass", "int(Tau_pt.size())")
               .Define("cleantau4vecs", ::gen4vec, {"Tau_pt", "Tau_eta", "Tau_phi", "Tau_mass"});

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

        auto tauSFIdVsJet = [this, _tauidSFjet, tauid_vsjet, tauid_vse, systNamesVsJet, nSystVsJet](floats &pt, uchars &dm, uchars &genmatch)->floatsVec {
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
                        float nom_sf = 1.0f;
                        try {
                            nom_sf = _tauidSFjet->evaluate({pt[i], int(dm[i]), gm, tauid_vsjet, tauid_vse, "nom", "dm"});
                        } catch (...) {
                            nom_sf = 1.0f;
                        }
                        for (const auto &sname : systNamesVsJet) {
                            try {
                                float val = _tauidSFjet->evaluate({pt[i], int(dm[i]), gm, tauid_vsjet, tauid_vse, sname, "dm"});
                                uncSources.emplace_back(val);
                            } catch (...) {
                                uncSources.emplace_back(nom_sf);
                            }
                        }
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
                        try {
                            uncSources.emplace_back(_tauidSFele->evaluate({eta[i], int(dm[i]), gm, tauid_vse, "nom"}));
                            uncSources.emplace_back(_tauidSFele->evaluate({eta[i], int(dm[i]), gm, tauid_vse, "up"}));
                            uncSources.emplace_back(_tauidSFele->evaluate({eta[i], int(dm[i]), gm, tauid_vse, "down"}));
                        } catch (const std::exception &e) {
                            std::cerr << "[WARN] tauSFIdVsEl evaluate failed (eta=" << eta[i]
                                      << " dm=" << int(dm[i]) << " gm=" << gm << "): " << e.what() << "\n";
                            uncSources.assign(3, 1.0f);
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
        auto tauSFIdVsMu = [this, _tauidSFmu, tauid_vsmu, tauid_vse, tauid_vsjet_for_vsmu, isVsMuSchema6](floats &eta, uchars &genmatch)->floatsVec {
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
                            std::cerr << "[WARN] tauSFIdVsMu evaluate failed (eta=" << eta[i]
                                      << " gm=" << gm << "): " << e.what() << "\n";
                            uncSources.assign(3, 1.0f);
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
        // A NEW lambda is needed because tauSFIdVsJet captures tauid_vsjet by value at
        // definition time (Tight WP). Reassigning the local variable has no effect on it.
        const std::string tauid_vsjet_loose = "Loose";
        auto tauSFIdVsJet_loose = [this, _tauidSFjet, tauid_vsjet_loose, tauid_vse, systNamesVsJet, nSystVsJet](floats &pt, uchars &dm, uchars &genmatch)->floatsVec {
            floats uncSources;
            uncSources.reserve(nSystVsJet);
            floatsVec wVec;
            size_t n = std::min({pt.size(), dm.size(), genmatch.size()});
            wVec.reserve(n > 0 ? n : 1);
            if (n > 0) {
                for (size_t i=0; i<n; i++) {
                    int gm = int(genmatch[i]);
                    if (gm == 5) {
                        float nom_sf = 1.0f;
                        try {
                            nom_sf = _tauidSFjet->evaluate({pt[i], int(dm[i]), gm, tauid_vsjet_loose, tauid_vse, "nom", "dm"});
                        } catch (...) { nom_sf = 1.0f; }
                        for (const auto &sname : systNamesVsJet) {
                            try {
                                uncSources.emplace_back(_tauidSFjet->evaluate({pt[i], int(dm[i]), gm, tauid_vsjet_loose, tauid_vse, sname, "dm"}));
                            } catch (...) {
                                uncSources.emplace_back(nom_sf);
                            }
                        }
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
                   .Define("tauWeightIdVsMu_loose",  tauSFIdVsMu,
                            {"Tau_eta_loose","Tau_genPartFlav_loose"});

        _rlm = _rlm.Define("taugencut_loose","Tau_genPartFlav_loose == 5")
                   .Redefine("Tau_pt_loose_gen","Tau_pt_loose[taugencut_loose]");
    }
}

void NanoAODAnalyzerrdframe::selectJets(std::vector<std::string> jes_var, std::vector<std::string> jes_var_flav, std::string cut) {
    //https://twiki.cern.ch/twiki/bin/view/CMS/JetID13p6TeV
    //nanoAOD Flags
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
        _rlm = _rlm.Redefine("Jet_hadronFlavour", "Jet_hadronFlavour[jetoverlap]");
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


    // b-tagging
    // https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer22/#ak4-b-tagging
    if (_isRun22) { 
        _rlm = _rlm.Define("btagcuts", "Jet_btagPNetB>0.245") //l: 0.047, m: 0.245, t: 0.6734
                   .Define("btagcuts_loose", "Jet_btagPNetB_loose>0.047");
    } else if (_isRun22EE) { 
        _rlm = _rlm.Define("btagcuts", "Jet_btagPNetB>0.2605") //l: 0.0499, m: 0.2605, t: 0.6915
                   .Define("btagcuts_loose", "Jet_btagPNetB_loose>0.0499");
    } else if (_isRun23) { 
        _rlm = _rlm.Define("btagcuts", "Jet_btagPNetB>0.1917") //l: 0.0358, m: 0.1917, t: 0.6172
                   .Define("btagcuts_loose", "Jet_btagPNetB_loose>0.0358");
    } else if (_isRun23BPix) { 
        _rlm = _rlm.Define("btagcuts", "Jet_btagPNetB>0.1919") //l: 0.0359, m: 0.1919, t: 0.6133
                   .Define("btagcuts_loose", "Jet_btagPNetB_loose>0.0359");
    } else if (_isRun24){
        // TODO: copy from 23BPix
        _rlm = _rlm.Define("btagcuts", "Jet_btagUParTAK4B>0.1272") //l: 0.0246, m: 0.1272, t: 0.4648
                   .Define("btagcuts_loose", "Jet_btagUParTAK4B_loose>0.0246");
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


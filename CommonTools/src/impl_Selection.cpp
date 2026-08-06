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
    cout << "Loading Btag SF: " << btagYear << endl;

    _rlm = _rlm.Define("btag_var", [](){return strings(btag_var);})
               .Define("btag_jes_var", [jes_var](){return strings(jes_var);});

    auto bSFreader = loadCorrectionSet("data/BTV/" + btagYear + "/btagging.json.gz");
    auto _btagSF = bSFreader->at(btagMap);
    auto _btagSFlight = bSFreader->at(btagMapLight);

    

    //case3 - Shape correction
    //If you are interested in using the whole b-tagging discriminant distribution in your analysis,
    //e.g. using b-tagging variables to separate signal and background, then this method is for you
    auto btagSF_shape = [this, _btagSF](floats &pts, floats &etas, uchars &hadflav, floats &btags)->floatsVec{
        std::vector<std::string> systs = {"central",
            "up_hf", "down_hf", "up_lf", "down_lf",
            "up_hfstats1", "down_hfstats1", "up_hfstats2", "down_hfstats2",
            "up_lfstats1", "down_lfstats1", "up_lfstats2", "down_lfstats2",
            "up_cferr1", "down_cferr1", "up_cferr2", "down_cferr2"};
        
        floats wVec;
        wVec.reserve(systs.size());
        floatsVec out;
        out.reserve(pts.size());

        if (int(pts.size()) != int(etas.size())) cout << "eta size hmmmmmmmmmmmm" << endl;
        if (int(pts.size()) != int(hadflav.size())) cout << "hadflav size hmmmmmmmmmmmm" << endl;
        if (int(pts.size()) != int(btags.size())) cout << "btag size hmmmmmmmmmmmm" << endl;
        for (auto i=0; i<int(pts.size()); i++){
            for (auto &syst: systs){
                float sf = 1.0;
                if (pts[i] < 30 || etas[i] > 2.5) {
                    sf = 1.0;
                } else if (btags[i] < 0.0f) {
                    // disc < 0 (e.g. btagPNetB=-1): jet is untaggable/masked.
                    // Correctionlib binning starts at 0; evaluating would throw.
                    // SF = 1.0 by convention (no correction applied).
                    sf = 1.0;
                } else if (syst.find("cferr")!=std::string::npos && hadflav[i]!=4) {
                    sf = 1.0;
                } else if ((syst.find("hf")!=std::string::npos || syst.find("lf")!=std::string::npos) && hadflav[i]==4) {
                    sf = 1.0;
                } else {
                    try {
                        sf = _btagSF->evaluate({syst, int(hadflav[i]), abs(etas[i]), pts[i], btags[i]});
                    } catch (const std::exception &e) {
                        std::cerr << "[WARN] btagSF evaluate failed (syst=" << syst
                                  << " flav=" << int(hadflav[i]) << " eta=" << etas[i]
                                  << " pt=" << pts[i] << " disc=" << btags[i]
                                  << "): " << e.what() << "\n";
                        sf = 1.0;
                    }
                }
                wVec.emplace_back(sf);
            }
            out.emplace_back(wVec);
            wVec.clear();
        }
        return out;
    };
    auto btag_evWeight = [this](floatsVec &btagWeights)->floats{
        const int vars = 17;
        floats out(vars, 1.0);

        for (auto &jet: btagWeights){
            for (int i=0; i<vars; i++){
                out[i] *= jet[i];
            }
        }
        return out;
    };

    _rlm = _rlm.Define("Jet_btagSF", btagSF_shape, {"Jet_pt", "Jet_eta", "Jet_hadronFlavour", "Jet_btagPNetB"})
               .Define("btagWeight", btag_evWeight, {"Jet_btagSF"});


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
        for (size_t i=0; i<toSkim.size(); i++) {
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
        std::string tauid_vsjet = "VTight";
        std::string tauid_vse   = "Tight";
        std::string tauid_vsmu  = "Tight";

        cout << "Tau ID WP vsJet : " << tauid_vsjet << endl;
        cout << "Tau ID WP vsMuon : " << tauid_vsmu << endl;
        cout << "Tau ID WP vsElectron : " << tauid_vse << endl;

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
            wVec.reserve(pt.size());

            if (pt.size() > 0) {
                for (size_t i=0; i<pt.size(); i++) {
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
            wVec.reserve(eta.size());

            if (eta.size() > 0) {
                for (size_t i=0; i<eta.size(); i++) {
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
        auto tauSFIdVsMu = [this, _tauidSFmu, tauid_vsmu, tauid_vse, tauid_vsjet](floats &eta, uchars &genmatch)->floatsVec {
            // CMS tau POG: DeepTauVSmu SF applies to muons faking tau:
            //   genmatch=2 (prompt μ), genmatch=4 (tau→μ).
            // Other genmatch values are not in the JSON → exception.
            floats uncSources;
            uncSources.reserve(3);
            floatsVec wVec;
            wVec.reserve(eta.size());

            if (eta.size() > 0) {
                for (size_t i=0; i<eta.size(); i++) {
                    int gm = int(genmatch[i]);
                    if (gm == 2 || gm == 4) {  // prompt μ or tau→μ
                        try {
                            if (_isRun24) {
                                uncSources.emplace_back(_tauidSFmu->evaluate({eta[i], gm, tauid_vsmu, tauid_vse, tauid_vsjet, "nom"}));
                                uncSources.emplace_back(_tauidSFmu->evaluate({eta[i], gm, tauid_vsmu, tauid_vse, tauid_vsjet, "up"}));
                                uncSources.emplace_back(_tauidSFmu->evaluate({eta[i], gm, tauid_vsmu, tauid_vse, tauid_vsjet, "down"}));
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


        _rlm = _rlm.Define("tauWeightIdVsJet", tauSFIdVsJet, {"Tau_pt","Tau_decayMode","Tau_genPartFlav"})
                   .Define("tauWeightIdVsEl", tauSFIdVsEl, {"Tau_eta","Tau_decayMode", "Tau_genPartFlav"})
                   .Define("tauWeightIdVsMu", tauSFIdVsMu, {"Tau_eta","Tau_genPartFlav"});

        tauid_vsjet = "Loose";

        _rlm = _rlm.Define("tauWeightIdVsJet_loose", tauSFIdVsJet, {"Tau_pt","Tau_decayMode","Tau_genPartFlav"});
        _rlm = _rlm.Define("Tau_genPartFlav_loose","Tau_genPartFlav[seltaucuts_loose]")
                   .Define("taugencut_loose","Tau_genPartFlav_loose == 5")
                   .Redefine("Tau_pt_loose_gen", "Tau_pt_loose_gen[taugencut_loose]")
                   .Redefine("tauWeightIdVsJet_loose", skimCol, {"tauWeightIdVsJet_loose", "seltaucuts_loose"})
                   .Define("tauWeightIdVsEl_loose", skimCol, {"tauWeightIdVsEl", "seltaucuts_loose"})
                   .Define("tauWeightIdVsMu_loose", skimCol, {"tauWeightIdVsMu", "seltaucuts_loose"});

        _rlm = _rlm.Redefine("Tau_genPartFlav","Tau_genPartFlav[seltaucuts]")
                   .Define("taugencut","Tau_genPartFlav == 5")
                   .Redefine("Tau_pt_gen", "Tau_pt[taugencut]")
                   .Redefine("tauWeightIdVsJet", skimCol, {"tauWeightIdVsJet", "seltaucuts"})
                   .Redefine("tauWeightIdVsEl", skimCol, {"tauWeightIdVsEl", "seltaucuts"})
                   .Redefine("tauWeightIdVsMu", skimCol, {"tauWeightIdVsMu", "seltaucuts"});
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
                   .Define("btagcuts_loose", "Jet_btagPNetB>0.047");
    } else if (_isRun22EE) { 
        _rlm = _rlm.Define("btagcuts", "Jet_btagPNetB>0.2605") //l: 0.0499, m: 0.2605, t: 0.6915
                   .Define("btagcuts_loose", "Jet_btagPNetB>0.0499");
    } else if (_isRun23) { 
        _rlm = _rlm.Define("btagcuts", "Jet_btagPNetB>0.1917") //l: 0.0358, m: 0.1917, t: 0.6172
                   .Define("btagcuts_loose", "Jet_btagPNetB>0.0358");
    } else if (_isRun23BPix) { 
        _rlm = _rlm.Define("btagcuts", "Jet_btagPNetB>0.1919") //l: 0.0359, m: 0.1919, t: 0.6133
                   .Define("btagcuts_loose", "Jet_btagPNetB>0.0359");
    } else if (_isRun24){
        // TODO: copy from 23BPix
        _rlm = _rlm.Define("btagcuts", "Jet_btagUParTAK4B>0.1272") //l: 0.0246, m: 0.1272, t: 0.4648
                   .Define("btagcuts_loose", "Jet_btagUParTAK4B>0.0246");
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


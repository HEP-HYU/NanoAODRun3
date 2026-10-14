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
#include "MuonScaReHelper.h"
#include "ElectronSSHelper.h"
#include "TCanvas.h"
#include "ROOT/RDFHelpers.hxx"
#include "Math/GenVector/VectorUtil.h"
using namespace std;

void NanoAODAnalyzerrdframe::calculateMuonSF(string muonid, string muoniso, string muonhlt){
    //// Muon SF
    cout<<"Loading Muon SF"<<endl;
    std::string muonEraFolder = "";

    if (_isRun22) {
        muonEraFolder = "2022_Summer22";
    } else if (_isRun22EE) {
        muonEraFolder = "2022_Summer22EE";
    } else if (_isRun23) {
        muonEraFolder = "2023_Summer23";
    } else if (_isRun23BPix) {
        muonEraFolder = "2023_Summer23BPix";
    } else if (_isRun24){
        muonEraFolder = "2024_Summer24";
    }

    auto muonSFreader = loadCorrectionSet("data/MuonSF/" + muonEraFolder + "/muon_Z.json.gz");
    auto _muonid = muonSFreader->at(muonid);
    auto _muoniso = muonSFreader->at(muoniso);
    auto _muontrg = muonSFreader->at(muonhlt);

    // We have only one muon!
    auto muonSFId = [this, _muonid, muonEraFolder](floats &pt, floats &eta)->floats {
        floats wVec;
        wVec.reserve(7); // [0]=nom, [1]=tot_up, [2]=tot_dn, [3]=stat_up, [4]=stat_dn, [5]=syst_up, [6]=syst_dn

        if (pt.size() == 1){
            float tmp_eta = (muonEraFolder.find("2022") != std::string::npos) ? std::abs(eta[0]) : eta[0];
            float sf = _muonid->evaluate({tmp_eta, pt[0], "nominal"});
            float sf_up = _muonid->evaluate({tmp_eta, pt[0], "systup"});
            float sf_down = _muonid->evaluate({tmp_eta, pt[0], "systdown"});
            float sf_stat = _muonid->evaluate({tmp_eta, pt[0], "stat"});
            float sf_syst = _muonid->evaluate({tmp_eta, pt[0], "syst"});
            wVec.emplace_back(sf);
            wVec.emplace_back(sf_up);
            wVec.emplace_back(sf_down);
            wVec.emplace_back(sf + sf_stat);
            wVec.emplace_back(sf - sf_stat);
            wVec.emplace_back(sf + sf_syst);
            wVec.emplace_back(sf - sf_syst);
        } else wVec = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
        return wVec;
    };

    // We have only one muon!
    auto muonSFIso = [this, _muoniso, muonEraFolder](floats &pt, floats &eta)->floats {
        floats wVec;
        wVec.reserve(7); // [0]=nom, [1]=tot_up, [2]=tot_dn, [3]=stat_up, [4]=stat_dn, [5]=syst_up, [6]=syst_dn

        if (pt.size() == 1){
            float tmp_eta = (muonEraFolder.find("2022") != std::string::npos) ? std::abs(eta[0]) : eta[0];
            float sf = _muoniso->evaluate({tmp_eta, pt[0], "nominal"});
            float sf_up = _muoniso->evaluate({tmp_eta, pt[0], "systup"});
            float sf_down = _muoniso->evaluate({tmp_eta, pt[0], "systdown"});
            float sf_stat = _muoniso->evaluate({tmp_eta, pt[0], "stat"});
            float sf_syst = _muoniso->evaluate({tmp_eta, pt[0], "syst"});
            wVec.emplace_back(sf);
            wVec.emplace_back(sf_up);
            wVec.emplace_back(sf_down);
            wVec.emplace_back(sf + sf_stat);
            wVec.emplace_back(sf - sf_stat);
            wVec.emplace_back(sf + sf_syst);
            wVec.emplace_back(sf - sf_syst);
        } else wVec = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
        return wVec;
    };

    auto muonSFTrg = [this, _muontrg, muonEraFolder](floats &pt, floats &eta)->floats {
        floats wVec;
        wVec.reserve(7); // [0]=nom, [1]=tot_up, [2]=tot_dn, [3]=stat_up, [4]=stat_dn, [5]=syst_up, [6]=syst_dn

        if (pt.size() == 1) {
            float tmp_eta = (muonEraFolder.find("2022") != std::string::npos) ? std::abs(eta[0]) : eta[0];
            float sf = _muontrg->evaluate({tmp_eta, pt[0], "nominal"});
            float sf_up = _muontrg->evaluate({tmp_eta, pt[0], "systup"});
            float sf_down = _muontrg->evaluate({tmp_eta, pt[0], "systdown"});
            float sf_stat = _muontrg->evaluate({tmp_eta, pt[0], "stat"});
            float sf_syst = _muontrg->evaluate({tmp_eta, pt[0], "syst"});
            wVec.emplace_back(sf);
            wVec.emplace_back(sf_up);
            wVec.emplace_back(sf_down);
            wVec.emplace_back(sf + sf_stat);
            wVec.emplace_back(sf - sf_stat);
            wVec.emplace_back(sf + sf_syst);
            wVec.emplace_back(sf - sf_syst);
        }
        else wVec = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
        return wVec;
    };

    if (isDefined("muonWeightId")) {
        _rlm = _rlm.Redefine("muonWeightId", muonSFId, {"Muon_pt","Muon_eta"})
                   .Redefine("muonWeightIso", muonSFIso, {"Muon_pt","Muon_eta"})
                   .Redefine("muonWeightTrg", muonSFTrg, {"Muon_pt","Muon_eta"});
    } else {
        _rlm = _rlm.Define("muonWeightId", muonSFId, {"Muon_pt","Muon_eta"})
                   .Define("muonWeightIso", muonSFIso, {"Muon_pt","Muon_eta"})
                   .Define("muonWeightTrg", muonSFTrg, {"Muon_pt","Muon_eta"});
    }
    
    // GEScaleSyst only has data for Run 2 UL eras (2016_UL_HIPM/2016_UL/2017_UL/2018_UL).
    // For Run 3, no official high-pT recommendation exists yet — return identity.
    // TODO: Replace the guard below with the Run 3 GE scale map once available.
    const bool isRun2GEera = (_year == "2016_UL_HIPM" || _year == "2016_UL" ||
                              _year == "2017_UL"       || _year == "2018_UL");
    std::string muonYear = _year;

    auto muonhighscaleup = [muonYear, isRun2GEera](floats &pts, floats &etas, floats &phis, ints &charges)->floats {
        floats out;
        out.reserve(pts.size());
        for (size_t i=0; i<pts.size(); i++) {
            float pt_tmp = pts[i];
            float pt_out = pt_tmp;
            if (isRun2GEera && pt_tmp > 200) {
                float eta_tmp = etas[i];
                float phi_tmp = phis[i];
                int charge_tmp = charges[i];
                GEScaleSyst GE(muonYear);
                GE.SetVerbose(0);
                pt_out = GE.GEScaleCorrPt(pt_tmp, eta_tmp, phi_tmp, charge_tmp, 0, 1);
            }
            out.emplace_back(pt_out);
        }
        return out;
    };

    auto muonhighscaledn = [muonYear, isRun2GEera](floats &pts, floats &etas, floats &phis, ints &charges)->floats {
        floats out;
        out.reserve(pts.size());
        for (size_t i=0; i<pts.size(); i++) {
            float pt_tmp = pts[i];
            float pt_out = pt_tmp;
            if (isRun2GEera && pt_tmp > 200) {
                float eta_tmp = etas[i];
                float phi_tmp = phis[i];
                int charge_tmp = charges[i];
                GEScaleSyst GE(muonYear);
                GE.SetVerbose(0);
                pt_out = GE.GEScaleCorrPt(pt_tmp, eta_tmp, phi_tmp, charge_tmp, 0, 2);
            }
            out.emplace_back(pt_out);
        }
        return out;
    };

    auto muonhighscalemet = [](floats &pts, floats &ptcors, float met)->float {
        float out = met;
        for (size_t i=0; i<pts.size(); i++) {
            out -= ptcors[i] - pts[i];
        }
        return out;
    };

    auto muonhighscalemetphi = [](floats &pts, floats &ptcors, floats &phis, float met, float metphi)->float {
        float out = 0.;
        auto metx = met * cos(metphi);
        auto mety = met * sin(metphi);
        for (size_t i=0; i<pts.size(); i++) {
            metx -= (ptcors[i] - pts[i]) * cos(phis[i]);
            mety -= (ptcors[i] - pts[i]) * sin(phis[i]);
        }
        out = float(atan2(mety, metx));
        return out;
    };

    if (_syst.find("muonhighscaleup") != std::string::npos) {
        _rlm = _rlm.Define("Muon_pt_scale", muonhighscaleup, {"Muon_pt", "Muon_eta", "Muon_phi", "Muon_charge"})
                   .Redefine("MET_phi",      muonhighscalemetphi, {"Muon_pt", "Muon_pt_scale", "Muon_phi", "MET_pt",      "MET_phi"})
                   .Redefine("MET_pt",       muonhighscalemet,    {"Muon_pt", "Muon_pt_scale",             "MET_pt"})
                   .Redefine("PuppiMET_phi", muonhighscalemetphi, {"Muon_pt", "Muon_pt_scale", "Muon_phi", "PuppiMET_pt", "PuppiMET_phi"})
                   .Redefine("PuppiMET_pt",  muonhighscalemet,    {"Muon_pt", "Muon_pt_scale",             "PuppiMET_pt"})
                   .Redefine("Muon_pt", "Muon_pt_scale"); //order matters
    } else if (_syst.find("muonhighscaledown") != std::string::npos) {
        _rlm = _rlm.Define("Muon_pt_scale", muonhighscaledn, {"Muon_pt", "Muon_eta", "Muon_phi", "Muon_charge"})
                   .Redefine("MET_phi",      muonhighscalemetphi, {"Muon_pt", "Muon_pt_scale", "Muon_phi", "MET_pt",      "MET_phi"})
                   .Redefine("MET_pt",       muonhighscalemet,    {"Muon_pt", "Muon_pt_scale",             "MET_pt"})
                   .Redefine("PuppiMET_phi", muonhighscalemetphi, {"Muon_pt", "Muon_pt_scale", "Muon_phi", "PuppiMET_pt", "PuppiMET_phi"})
                   .Redefine("PuppiMET_pt",  muonhighscalemet,    {"Muon_pt", "Muon_pt_scale",             "PuppiMET_pt"})
                   .Redefine("Muon_pt", "Muon_pt_scale"); //order matters
    }

}

void NanoAODAnalyzerrdframe::applyMuonScaleSmearing() {
    std::string muonEraFolder = "";
    if (_isRun22) {
        muonEraFolder = "2022_Summer22";
    } else if (_isRun22EE) {
        muonEraFolder = "2022_Summer22EE";
    } else if (_isRun23) {
        muonEraFolder = "2023_Summer23";
    } else if (_isRun23BPix) {
        muonEraFolder = "2023_Summer23BPix";
    } else if (_isRun24){
        muonEraFolder = "2024_Summer24";
    }

    if (muonEraFolder.empty()) {
        std::cout << "[ScaReKIT] Era not recognized for Run 3 muon scale/smearing, skipping." << std::endl;
        return;
    }

    std::string scareFile = "data/MuonSF/" + muonEraFolder + "/muon_scalesmearing.json.gz";
    std::cout << "Loading Muon Scale and Smearing (ScaReKIT): " << scareFile << " [syst=" << _syst << "]" << std::endl;

    auto scareReader = loadCorrectionSet(scareFile);
    auto helper = std::make_shared<MuonScaReHelper>(scareReader);
    const bool isData = _isData;
    const std::string systStr = _syst;
    const bool hasNTrackerLayers = isDefined("Muon_nTrackerLayers");

    auto muonScareLambda = [helper, isData, systStr, hasNTrackerLayers](floats &pts, floats &etas, floats &phis, ints &charges,
                                                                        unsigned long long evt, unsigned int lumi,
                                                                        floats &nLayers)->floats {
        floats out;
        out.reserve(pts.size());
        for (size_t i = 0; i < pts.size(); i++) {
            float pt = pts[i];
            float eta = etas[i];
            float phi = phis[i];
            int charge = charges[i];
            float nL = (hasNTrackerLayers && i < nLayers.size()) ? nLayers[i] : 11.0f;

            float pt_final = pt;
            if (isData) {
                pt_final = helper->pt_scale(true, pt, eta, phi, charge);
            } else {
                double pt_scale_corr = helper->pt_scale(false, pt, eta, phi, charge);
                double pt_corr = helper->pt_resol(pt_scale_corr, eta, phi, nL, evt, lumi);

                if (systStr.find("muonscaleup") != std::string::npos) {
                    pt_final = helper->pt_scale_var(pt_corr, eta, phi, charge, "up");
                } else if (systStr.find("muonscaledown") != std::string::npos) {
                    pt_final = helper->pt_scale_var(pt_corr, eta, phi, charge, "dn");
                } else if (systStr.find("muonsmearup") != std::string::npos) {
                    pt_final = helper->pt_resol_var(pt_scale_corr, pt_corr, eta, "up");
                } else if (systStr.find("muonsmeardown") != std::string::npos) {
                    pt_final = helper->pt_resol_var(pt_scale_corr, pt_corr, eta, "dn");
                } else {
                    pt_final = pt_corr;
                }
            }
            out.emplace_back(pt_final);
        }
        return out;
    };

    auto muonScareLambdaNoNLayers = [helper, isData, systStr](floats &pts, floats &etas, floats &phis, ints &charges,
                                                             unsigned long long evt, unsigned int lumi)->floats {
        floats out;
        out.reserve(pts.size());
        for (size_t i = 0; i < pts.size(); i++) {
            float pt = pts[i];
            float eta = etas[i];
            float phi = phis[i];
            int charge = charges[i];
            float nL = 11.0f; // median fallback for tight muons

            float pt_final = pt;
            if (isData) {
                pt_final = helper->pt_scale(true, pt, eta, phi, charge);
            } else {
                double pt_scale_corr = helper->pt_scale(false, pt, eta, phi, charge);
                double pt_corr = helper->pt_resol(pt_scale_corr, eta, phi, nL, evt, lumi);

                if (systStr.find("muonscaleup") != std::string::npos) {
                    pt_final = helper->pt_scale_var(pt_corr, eta, phi, charge, "up");
                } else if (systStr.find("muonscaledown") != std::string::npos) {
                    pt_final = helper->pt_scale_var(pt_corr, eta, phi, charge, "dn");
                } else if (systStr.find("muonsmearup") != std::string::npos) {
                    pt_final = helper->pt_resol_var(pt_scale_corr, pt_corr, eta, "up");
                } else if (systStr.find("muonsmeardown") != std::string::npos) {
                    pt_final = helper->pt_resol_var(pt_scale_corr, pt_corr, eta, "dn");
                } else {
                    pt_final = pt_corr;
                }
            }
            out.emplace_back(pt_final);
        }
        return out;
    };

    auto muonscaremet = [](floats &pts, floats &ptcors, floats &phis, float met, float metphi)->float {
        double metx = met * std::cos(metphi);
        double mety = met * std::sin(metphi);
        for (size_t i = 0; i < pts.size(); i++) {
            double dpx = (ptcors[i] - pts[i]) * std::cos(phis[i]);
            double dpy = (ptcors[i] - pts[i]) * std::sin(phis[i]);
            metx -= dpx;
            mety -= dpy;
        }
        return float(std::sqrt(metx * metx + mety * mety));
    };

    auto muonscaremetphi = [](floats &pts, floats &ptcors, floats &phis, float met, float metphi)->float {
        double metx = met * std::cos(metphi);
        double mety = met * std::sin(metphi);
        for (size_t i = 0; i < pts.size(); i++) {
            double dpx = (ptcors[i] - pts[i]) * std::cos(phis[i]);
            double dpy = (ptcors[i] - pts[i]) * std::sin(phis[i]);
            metx -= dpx;
            mety -= dpy;
        }
        return float(std::atan2(mety, metx));
    };

    if (hasNTrackerLayers) {
        _rlm = _rlm.Define("Muon_nTrackerLayers_f", "ROOT::VecOps::Cast<float>(Muon_nTrackerLayers)")
                   .Define("Muon_pt_scare", muonScareLambda, {"Muon_pt", "Muon_eta", "Muon_phi", "Muon_charge", "event", "luminosityBlock", "Muon_nTrackerLayers_f"});
    } else {
        _rlm = _rlm.Define("Muon_pt_scare", muonScareLambdaNoNLayers, {"Muon_pt", "Muon_eta", "Muon_phi", "Muon_charge", "event", "luminosityBlock"});
    }

    _rlm = _rlm.Redefine("PuppiMET_phi", muonscaremetphi, {"Muon_pt", "Muon_pt_scare", "Muon_phi", "PuppiMET_pt", "PuppiMET_phi"})
               .Redefine("PuppiMET_pt",  muonscaremet,    {"Muon_pt", "Muon_pt_scare", "Muon_phi", "PuppiMET_pt", "PuppiMET_phi"})
               .Redefine("Muon_pt", "Muon_pt_scare");

    if (isDefined("lep4vecs")) {
        _rlm = _rlm.Redefine("lep4vecs", ::gen4vec, {"Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass"});
    }
}
 
void NanoAODAnalyzerrdframe::applyElectronScaleSmearing() {
    std::string elecEraFolder = "";
    if (_isRun22) {
        elecEraFolder = "2022_Summer22";
    } else if (_isRun22EE) {
        elecEraFolder = "2022_Summer22EE";
    } else if (_isRun23) {
        elecEraFolder = "2023_Summer23";
    } else if (_isRun23BPix) {
        elecEraFolder = "2023_Summer23BPix";
    } else if (_isRun24){
        elecEraFolder = "2024_Summer24";
    }

    if (elecEraFolder.empty()) {
        std::cout << "[EGM S&S] Era not recognized for Run 3 electron scale/smearing, skipping." << std::endl;
        return;
    }

    std::string ssFile = "data/ElectronSF/" + elecEraFolder + "/electronSS_EtDependent.json.gz";
    std::cout << "Loading Electron Scale and Smearing: " << ssFile << " [syst=" << _syst << "]" << std::endl;

    auto ssReader = loadCorrectionSet(ssFile);
    auto helper = std::make_shared<ElectronSSHelper>(ssReader);
    const bool isData = _isData;
    const std::string systStr = _syst;
    const bool hasR9 = isDefined("Electron_r9");
    const bool hasDeltaEtaSC = isDefined("Electron_deltaEtaSC");

    auto eleScaleSmearLambda = [helper, isData, systStr](floats &pts, floats &etas, floats &phis, floats &deltaEtas, floats &r9s,
                                                         unsigned long long evt, unsigned int lumi)->floats {
        floats out;
        out.reserve(pts.size());
        for (size_t i = 0; i < pts.size(); i++) {
            float pt = pts[i];
            float eta = etas[i];
            float phi = phis[i];
            float deltaSC = (i < deltaEtas.size()) ? deltaEtas[i] : 0.0f;
            float r9 = (i < r9s.size()) ? r9s[i] : 0.98f;

            double pt_corr = helper->pt_corr(isData, systStr, pt, eta, phi, deltaSC, r9, evt, lumi, static_cast<int>(i));
            out.emplace_back(static_cast<float>(pt_corr));
        }
        return out;
    };

    auto eleScaleSmearLambdaNoR9 = [helper, isData, systStr](floats &pts, floats &etas, floats &phis, floats &deltaEtas,
                                                             unsigned long long evt, unsigned int lumi)->floats {
        floats out;
        out.reserve(pts.size());
        for (size_t i = 0; i < pts.size(); i++) {
            float pt = pts[i];
            float eta = etas[i];
            float phi = phis[i];
            float deltaSC = (i < deltaEtas.size()) ? deltaEtas[i] : 0.0f;
            float r9 = 0.98f; // graceful fallback for existing skims (falls in High-R9 bin >= 0.975)

            double pt_corr = helper->pt_corr(isData, systStr, pt, eta, phi, deltaSC, r9, evt, lumi, static_cast<int>(i));
            out.emplace_back(static_cast<float>(pt_corr));
        }
        return out;
    };

    auto eleScaleSmearLambdaNoDeltaNoR9 = [helper, isData, systStr](floats &pts, floats &etas, floats &phis,
                                                                   unsigned long long evt, unsigned int lumi)->floats {
        floats out;
        out.reserve(pts.size());
        for (size_t i = 0; i < pts.size(); i++) {
            float pt = pts[i];
            float eta = etas[i];
            float phi = phis[i];
            float deltaSC = 0.0f;
            float r9 = 0.98f;

            double pt_corr = helper->pt_corr(isData, systStr, pt, eta, phi, deltaSC, r9, evt, lumi, static_cast<int>(i));
            out.emplace_back(static_cast<float>(pt_corr));
        }
        return out;
    };

    auto elescaremet = [](floats &pts, floats &ptcors, floats &phis, float met, float metphi)->float {
        double metx = met * std::cos(metphi);
        double mety = met * std::sin(metphi);
        for (size_t i = 0; i < pts.size(); i++) {
            double dpx = (ptcors[i] - pts[i]) * std::cos(phis[i]);
            double dpy = (ptcors[i] - pts[i]) * std::sin(phis[i]);
            metx -= dpx;
            mety -= dpy;
        }
        return float(std::sqrt(metx * metx + mety * mety));
    };

    auto elescaremetphi = [](floats &pts, floats &ptcors, floats &phis, float met, float metphi)->float {
        double metx = met * std::cos(metphi);
        double mety = met * std::sin(metphi);
        for (size_t i = 0; i < pts.size(); i++) {
            double dpx = (ptcors[i] - pts[i]) * std::cos(phis[i]);
            double dpy = (ptcors[i] - pts[i]) * std::sin(phis[i]);
            metx -= dpx;
            mety -= dpy;
        }
        return float(std::atan2(mety, metx));
    };

    if (hasR9 && hasDeltaEtaSC) {
        _rlm = _rlm.Define("Electron_pt_scare", eleScaleSmearLambda, {"Electron_pt", "Electron_eta", "Electron_phi", "Electron_deltaEtaSC", "Electron_r9", "event", "luminosityBlock"});
    } else if (hasDeltaEtaSC) {
        _rlm = _rlm.Define("Electron_pt_scare", eleScaleSmearLambdaNoR9, {"Electron_pt", "Electron_eta", "Electron_phi", "Electron_deltaEtaSC", "event", "luminosityBlock"});
    } else {
        _rlm = _rlm.Define("Electron_pt_scare", eleScaleSmearLambdaNoDeltaNoR9, {"Electron_pt", "Electron_eta", "Electron_phi", "event", "luminosityBlock"});
    }

    _rlm = _rlm.Redefine("PuppiMET_phi", elescaremetphi, {"Electron_pt", "Electron_pt_scare", "Electron_phi", "PuppiMET_pt", "PuppiMET_phi"})
               .Redefine("PuppiMET_pt",  elescaremet,    {"Electron_pt", "Electron_pt_scare", "Electron_phi", "PuppiMET_pt", "PuppiMET_phi"})
               .Redefine("Electron_pt", "Electron_pt_scare");

    if (isDefined("lep4vecs")) {
        _rlm = _rlm.Redefine("lep4vecs", ::gen4vec, {"Electron_pt", "Electron_eta", "Electron_phi", "Electron_mass"});
    }
}

void NanoAODAnalyzerrdframe::calculateElectronSF(string elecFile, string elecYear){
    // Electron SF
    cout << "Loading Electron SF" << endl;

    auto elecSFreader = loadCorrectionSet("data/ElectronSF/"+elecFile+"/electron.json.gz");
    auto _elecid = elecSFreader->at("Electron-ID-SF");

    auto elecSFId = [this, _elecid, elecYear](std::string var){
        return [this, _elecid, elecYear, var](floats &pt, floats &eta, floats &phi)->floats {
            floats wVec;
            wVec.reserve(7); // [0]: nom, [1]: tot_up, [2]: tot_down, [3]: stat_up, [4]: stat_down, [5]: syst_up, [6]: syst_down
            std::string wp = var;
            if (pt.size() == 1){
                float sf=1.0f; float sf_up=1.0f; float sf_down=1.0f;
                float stat_up=1.0f; float stat_down=1.0f;
                float syst_up=1.0f; float syst_down=1.0f;
                if (wp.find("Reco") != std::string::npos) {
                    if (pt[0] < 20) wp = "RecoBelow20";
                    else if (pt[0] < 75) wp = "Reco20to75";
                    else wp = "RecoAbove75";
                }
                if (elecYear.find("2022")!=std::string::npos || elecYear.find("2024")!=std::string::npos){
                    sf = _elecid->evaluate({elecYear, "sf", wp, eta[0], pt[0]});
                    sf_up = _elecid->evaluate({elecYear, "sfup", wp, eta[0], pt[0]});
                    sf_down = _elecid->evaluate({elecYear, "sfdown", wp, eta[0], pt[0]});
                    float err_stat = _elecid->evaluate({elecYear, "err_stat", wp, eta[0], pt[0]});
                    float err_syst = _elecid->evaluate({elecYear, "err_syst", wp, eta[0], pt[0]});
                    stat_up = sf + err_stat;
                    stat_down = sf - err_stat;
                    syst_up = sf + err_syst;
                    syst_down = sf - err_syst;
                } else if (elecYear.find("2023")!=std::string::npos){
                    sf = _elecid->evaluate({elecYear, "sf", wp, eta[0], pt[0], phi[0]});
                    sf_up = _elecid->evaluate({elecYear, "sfup", wp, eta[0], pt[0], phi[0]});
                    sf_down = _elecid->evaluate({elecYear, "sfdown", wp, eta[0], pt[0], phi[0]});
                    float err_stat = _elecid->evaluate({elecYear, "err_stat", wp, eta[0], pt[0], phi[0]});
                    float err_syst = _elecid->evaluate({elecYear, "err_syst", wp, eta[0], pt[0], phi[0]});
                    stat_up = sf + err_stat;
                    stat_down = sf - err_stat;
                    syst_up = sf + err_syst;
                    syst_down = sf - err_syst;
                }
                wVec.emplace_back(sf);
                wVec.emplace_back(sf_up);
                wVec.emplace_back(sf_down);
                wVec.emplace_back(stat_up);
                wVec.emplace_back(stat_down);
                wVec.emplace_back(syst_up);
                wVec.emplace_back(syst_down);
            } else {
                wVec = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
            }
            return wVec;
        };
    };

    if (isDefined("elecWeightReco")) {
        _rlm = _rlm.Redefine("elecWeightReco", elecSFId("RecoAbove75"), {"Electron_pt","Electron_eta","Electron_phi"})
                   .Redefine("elecWeightId", elecSFId("wp90iso"), {"Electron_pt","Electron_eta","Electron_phi"});
    } else {
        _rlm = _rlm.Define("elecWeightReco", elecSFId("RecoAbove75"), {"Electron_pt","Electron_eta","Electron_phi"})
                   .Define("elecWeightId", elecSFId("wp90iso"), {"Electron_pt","Electron_eta","Electron_phi"});
    }

    auto elecHltSFreader = loadCorrectionSet("data/ElectronSF/"+elecFile+"/electronHlt.json.gz");
    auto _elechlt = elecHltSFreader->at("Electron-HLT-SF");

    auto elecSFTrg = [this, _elechlt, elecYear](floats &pt, floats &eta)->floats {
        floats wVec;
        wVec.reserve(3); // cent, up, down
        string hltpath = "HLT_SF_Ele30_MVAiso90ID";
        if (pt.size() == 1){
            float sf = _elechlt->evaluate({elecYear, "sf", hltpath, eta[0], pt[0]});
            float sf_up = _elechlt->evaluate({elecYear, "sfup", hltpath, eta[0], pt[0]});
            float sf_down = _elechlt->evaluate({elecYear, "sfdown", hltpath, eta[0], pt[0]});
            wVec.emplace_back(sf);
            wVec.emplace_back(sf_up);
            wVec.emplace_back(sf_down);
        } else {
            wVec = {1.0f, 1.0f, 1.0f};
        }
        return wVec;
    };

    if (isDefined("elecWeightTrg")) {
        _rlm = _rlm.Redefine("elecWeightTrg", elecSFTrg, {"Electron_pt","Electron_eta"});
    } else {
        _rlm = _rlm.Define("elecWeightTrg", elecSFTrg, {"Electron_pt","Electron_eta"});
    }
}


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

    if (muonEraFolder.empty() || muonid.empty() || muoniso.empty() || muonhlt.empty()) {
        std::cerr << "[ERROR] calculateMuonSF: muonEraFolder ('" << muonEraFolder
                  << "'), muonid ('" << muonid << "'), muoniso ('" << muoniso
                  << "'), or muonhlt ('" << muonhlt << "') is empty!" << std::endl;
        return;
    }
    if (!isDefined("Muon_pt") || !isDefined("Muon_eta")) {
        std::cerr << "[ERROR] calculateMuonSF: Missing required Muon_pt or Muon_eta branch!" << std::endl;
        return;
    }

    auto muonSFreader = loadCorrectionSet("data/MuonSF/" + muonEraFolder + "/muon_Z.json.gz");
    auto _muonid = muonSFreader->at(muonid);
    auto _muoniso = muonSFreader->at(muoniso);
    auto _muontrg = muonSFreader->at(muonhlt);

    // We have only one muon!
    auto muonSFId = [this, _muonid](const floats &pt, const floats &eta)->floats {
        floats wVec;
        wVec.reserve(7); // [0]=nom, [1]=tot_up, [2]=tot_dn, [3]=stat_up, [4]=stat_dn, [5]=syst_up, [6]=syst_dn

        if (pt.size() == 1){
            if (eta.size() != pt.size()) {
                std::cerr << "[ERROR] calculateMuonSF (ID): Size mismatch between Muon_pt (" << pt.size()
                          << ") and Muon_eta (" << eta.size() << ")!" << std::endl;
                return {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
            }
            float sf      = _muonid->evaluate({eta[0], pt[0], "nominal"});
            float sf_up   = _muonid->evaluate({eta[0], pt[0], "systup"});
            float sf_down = _muonid->evaluate({eta[0], pt[0], "systdown"});
            float sf_stat = _muonid->evaluate({eta[0], pt[0], "stat"});
            float sf_syst = _muonid->evaluate({eta[0], pt[0], "syst"});
            wVec.emplace_back(sf);
            wVec.emplace_back(sf_up);
            wVec.emplace_back(sf_down);
            wVec.emplace_back(sf + sf_stat);
            wVec.emplace_back(sf - sf_stat);
            wVec.emplace_back(sf + sf_syst);
            wVec.emplace_back(sf - sf_syst);
        } else {
            wVec = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
        }
        return wVec;
    };

    // We have only one muon!
    auto muonSFIso = [this, _muoniso](const floats &pt, const floats &eta)->floats {
        floats wVec;
        wVec.reserve(7); // [0]=nom, [1]=tot_up, [2]=tot_dn, [3]=stat_up, [4]=stat_dn, [5]=syst_up, [6]=syst_dn

        if (pt.size() == 1){
            if (eta.size() != pt.size()) {
                std::cerr << "[ERROR] calculateMuonSF (Iso): Size mismatch between Muon_pt (" << pt.size()
                          << ") and Muon_eta (" << eta.size() << ")!" << std::endl;
                return {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
            }
            float sf      = _muoniso->evaluate({eta[0], pt[0], "nominal"});
            float sf_up   = _muoniso->evaluate({eta[0], pt[0], "systup"});
            float sf_down = _muoniso->evaluate({eta[0], pt[0], "systdown"});
            float sf_stat = _muoniso->evaluate({eta[0], pt[0], "stat"});
            float sf_syst = _muoniso->evaluate({eta[0], pt[0], "syst"});
            wVec.emplace_back(sf);
            wVec.emplace_back(sf_up);
            wVec.emplace_back(sf_down);
            wVec.emplace_back(sf + sf_stat);
            wVec.emplace_back(sf - sf_stat);
            wVec.emplace_back(sf + sf_syst);
            wVec.emplace_back(sf - sf_syst);
        } else {
            wVec = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
        }
        return wVec;
    };

    auto muonSFTrg = [this, _muontrg](const floats &pt, const floats &eta)->floats {
        floats wVec;
        wVec.reserve(7); // [0]=nom, [1]=tot_up, [2]=tot_dn, [3]=stat_up, [4]=stat_dn, [5]=syst_up, [6]=syst_dn

        if (pt.size() == 1) {
            if (eta.size() != pt.size()) {
                std::cerr << "[ERROR] calculateMuonSF (Trg): Size mismatch between Muon_pt (" << pt.size()
                          << ") and Muon_eta (" << eta.size() << ")!" << std::endl;
                return {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
            }
            float sf      = _muontrg->evaluate({eta[0], pt[0], "nominal"});
            float sf_up   = _muontrg->evaluate({eta[0], pt[0], "systup"});
            float sf_down = _muontrg->evaluate({eta[0], pt[0], "systdown"});
            float sf_stat = _muontrg->evaluate({eta[0], pt[0], "stat"});
            float sf_syst = _muontrg->evaluate({eta[0], pt[0], "syst"});
            wVec.emplace_back(sf);
            wVec.emplace_back(sf_up);
            wVec.emplace_back(sf_down);
            wVec.emplace_back(sf + sf_stat);
            wVec.emplace_back(sf - sf_stat);
            wVec.emplace_back(sf + sf_syst);
            wVec.emplace_back(sf - sf_syst);
        } else {
            wVec = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
        }
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
        std::cerr << "[ERROR] [ScaReKIT] Era not recognized for Run 3 muon scale/smearing!" << std::endl;
        return;
    }

    const bool isPostSkim = isDefined("nmuonpass");
    const bool hasPtUncor = isDefined("Muon_pt_uncor");

    if (isPostSkim) {
        if (!hasPtUncor) {
            std::cerr << "[ERROR] [ScaReKIT] Muon_pt_uncor branch is missing in skimmed input! Applying ScaReKIT on Muon_pt in Process stage for backward compatibility." << std::endl;
        } else if (_isData || (_syst.find("muonscale") == std::string::npos && _syst.find("muonsmear") == std::string::npos)) {
            std::cout << "[ScaReKIT] Nominal Muon_pt and PuppiMET already corrected in Skim [syst=" << _syst << "]." << std::endl;
            return;
        }
    }

    if (!isDefined("Muon_pt") || !isDefined("Muon_eta") || !isDefined("Muon_phi") ||
        !isDefined("Muon_charge") || !isDefined("Muon_nTrackerLayers") ||
        !isDefined("event") || !isDefined("luminosityBlock")) {
        std::cerr << "[ERROR] [ScaReKIT] Missing required Muon branches (Muon_pt, Muon_eta, Muon_phi, Muon_charge, Muon_nTrackerLayers, event, or luminosityBlock)!" << std::endl;
        return;
    }

    std::string scareFile = "data/MuonSF/" + muonEraFolder + "/muon_scalesmearing.json.gz";
    std::cout << "Loading Muon Scale and Smearing (ScaReKIT): " << scareFile << " [syst=" << _syst << "]" << std::endl;

    auto scareReader = loadCorrectionSet(scareFile);
    auto helper = std::make_shared<MuonScaReHelper>(scareReader);
    const bool isData = _isData;
    const std::string systStr = _syst;

    auto muonScareLambda = [helper, isData, systStr](const floats &pts, const floats &etas, const floats &phis, const ints &charges,
                                                     ULong64_t evt, UInt_t lumi,
                                                     const floats &nLayers)->floats {
        floats out;
        out.reserve(pts.size());
        if (etas.size() != pts.size() || phis.size() != pts.size() ||
            charges.size() != pts.size() || nLayers.size() != pts.size()) {
            std::cerr << "[ERROR] [ScaReKIT] Size mismatch among Muon branches: pt=" << pts.size()
                      << ", eta=" << etas.size() << ", phi=" << phis.size()
                      << ", charge=" << charges.size() << ", nTrackerLayers=" << nLayers.size() << std::endl;
            return pts;
        }
        for (size_t i = 0; i < pts.size(); i++) {
            float pt = pts[i];
            float eta = etas[i];
            float phi = phis[i];
            int charge = charges[i];
            float nL = nLayers[i];

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

    auto muonscaremet = [](const floats &pts, const floats &ptcors, const floats &phis, float met, float metphi)->float {
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

    auto muonscaremetphi = [](const floats &pts, const floats &ptcors, const floats &phis, float met, float metphi)->float {
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

    if (!isPostSkim) {
        // Skim stage (before selectMuons): preserve Muon_pt_uncor and apply nominal ScaReKIT to Muon_pt.
        // PuppiMET propagation is performed inside selectMuons() after [muoncuts] so only selected muons shift MET.
        _rlm = _rlm.Define("Muon_pt_uncor", "Muon_pt")
                   .Define("skim_Muon_nTrackerLayers_f", "ROOT::VecOps::Cast<float>(Muon_nTrackerLayers)")
                   .Define("skim_Muon_pt_scare", muonScareLambda, {"Muon_pt_uncor", "Muon_eta", "Muon_phi", "Muon_charge", "event", "luminosityBlock", "skim_Muon_nTrackerLayers_f"})
                   .Redefine("Muon_pt", "skim_Muon_pt_scare");
    } else {
        // Process stage: evaluate from Muon_pt_uncor (or Muon_pt for legacy skims) and shift PuppiMET relative to current Muon_pt.
        const std::string inputPtCol = hasPtUncor ? "Muon_pt_uncor" : "Muon_pt";
        _rlm = _rlm.Define("Muon_nTrackerLayers_f", "ROOT::VecOps::Cast<float>(Muon_nTrackerLayers)")
                   .Define("Muon_pt_scare", muonScareLambda, {inputPtCol, "Muon_eta", "Muon_phi", "Muon_charge", "event", "luminosityBlock", "Muon_nTrackerLayers_f"})
                   .Define("PuppiMET_pt_mu_unsmear", "PuppiMET_pt")
                   .Define("PuppiMET_phi_mu_unsmear", "PuppiMET_phi")
                   .Redefine("PuppiMET_phi", muonscaremetphi, {"Muon_pt", "Muon_pt_scare", "Muon_phi", "PuppiMET_pt_mu_unsmear", "PuppiMET_phi_mu_unsmear"})
                   .Redefine("PuppiMET_pt",  muonscaremet,    {"Muon_pt", "Muon_pt_scare", "Muon_phi", "PuppiMET_pt_mu_unsmear", "PuppiMET_phi_mu_unsmear"})
                   .Redefine("Muon_pt", "Muon_pt_scare");

        if (isDefined("lep4vecs")) {
            _rlm = _rlm.Redefine("lep4vecs", ::gen4vec, {"Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass"});
        }
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
        std::cerr << "[ERROR] [EGM S&S] Era not recognized for Run 3 electron scale/smearing!" << std::endl;
        return;
    }

    const bool isPostSkim = isDefined("nelepass");
    const bool hasPtUncor = isDefined("Electron_pt_uncor");

    if (isPostSkim) {
        if (!hasPtUncor) {
            std::cerr << "[ERROR] [EGM S&S] Electron_pt_uncor branch is missing in skimmed input!" << std::endl;
            return;
        }
        // In Process stage, Data and nominal/non-EGM-syst MC already have S&S applied to Electron_pt and PuppiMET in Skim.
        if (_isData || (_syst.find("elescale") == std::string::npos && _syst.find("elesmear") == std::string::npos)) {
            std::cout << "[EGM S&S] Nominal Electron_pt and PuppiMET already corrected in Skim [syst=" << _syst << "]." << std::endl;
            return;
        }
    }

    const bool hasR9 = isDefined("Electron_r9");
    const bool hasDeltaEtaSC = isDefined("Electron_deltaEtaSC");
    const bool hasSeedGain = isDefined("Electron_seedGain");
    if (!hasR9 || !hasDeltaEtaSC || !hasSeedGain) {
        std::cerr << "[ERROR] [EGM S&S] Missing required Electron S&S branch(es): "
                  << "Electron_r9=" << hasR9
                  << ", Electron_deltaEtaSC=" << hasDeltaEtaSC
                  << ", Electron_seedGain=" << hasSeedGain << std::endl;
        return;
    }

    std::string ssFile = "data/ElectronSF/" + elecEraFolder + "/electronSS_EtDependent.json.gz";
    std::cout << "Loading Electron Scale and Smearing: " << ssFile << " [syst=" << _syst << "]" << std::endl;

    auto ssReader = loadCorrectionSet(ssFile);
    auto helper = std::make_shared<ElectronSSHelper>(ssReader);
    const bool isData = _isData;
    const std::string systStr = _syst;

    auto eleScaleSmearLambda = [helper, isData, systStr](const floats &pts, const floats &etas, const floats &phis, const floats &deltaEtas, const floats &r9s, const ints &seedGains,
                                                         UInt_t run, ULong64_t evt, UInt_t lumi)->floats {
        floats out;
        out.reserve(pts.size());
        if (deltaEtas.size() != pts.size() || r9s.size() != pts.size() || seedGains.size() != pts.size()) {
            std::cerr << "[ERROR] [EGM S&S] Branch size mismatch with Electron_pt (" << pts.size()
                      << "): deltaEtaSC=" << deltaEtas.size()
                      << ", r9=" << r9s.size()
                      << ", seedGain=" << seedGains.size() << std::endl;
            return pts;
        }
        for (size_t i = 0; i < pts.size(); i++) {
            float pt = pts[i];
            float eta = etas[i];
            float phi = phis[i];
            float deltaSC = deltaEtas[i];
            float r9 = r9s[i];
            int seedGain = seedGains[i];

            double pt_corr = helper->pt_corr(isData, systStr, pt, eta, phi, deltaSC, r9, run, seedGain, evt, lumi);
            out.emplace_back(static_cast<float>(pt_corr));
        }
        return out;
    };

    auto elescaremet = [](const floats &pts, const floats &ptcors, const floats &phis, float met, float metphi)->float {
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

    auto elescaremetphi = [](const floats &pts, const floats &ptcors, const floats &phis, float met, float metphi)->float {
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

    if (!isPostSkim) {
        // Skim stage (before selectElectrons): preserve Electron_pt_uncor and apply nominal S&S to Electron_pt.
        // PuppiMET propagation is performed inside selectElectrons() after [elecuts] so only selected electrons shift MET.
        _rlm = _rlm.Define("Electron_pt_uncor", "Electron_pt")
                   .Define("Electron_seedGain_i", "ROOT::VecOps::Cast<int>(Electron_seedGain)")
                   .Define("Electron_pt_scare", eleScaleSmearLambda, {"Electron_pt_uncor", "Electron_eta", "Electron_phi", "Electron_deltaEtaSC", "Electron_r9", "Electron_seedGain_i", "run", "event", "luminosityBlock"})
                   .Redefine("Electron_pt", "Electron_pt_scare");
    } else {
        // Process stage (elescaleup/down or elesmearup/down): evaluate from Electron_pt_uncor and shift PuppiMET relative to nominal Electron_pt.
        _rlm = _rlm.Define("Electron_seedGain_i", "ROOT::VecOps::Cast<int>(Electron_seedGain)")
                   .Define("Electron_pt_scare", eleScaleSmearLambda, {"Electron_pt_uncor", "Electron_eta", "Electron_phi", "Electron_deltaEtaSC", "Electron_r9", "Electron_seedGain_i", "run", "event", "luminosityBlock"})
                   .Define("PuppiMET_pt_ele_unsmear", "PuppiMET_pt")
                   .Define("PuppiMET_phi_ele_unsmear", "PuppiMET_phi")
                   .Redefine("PuppiMET_phi", elescaremetphi, {"Electron_pt", "Electron_pt_scare", "Electron_phi", "PuppiMET_pt_ele_unsmear", "PuppiMET_phi_ele_unsmear"})
                   .Redefine("PuppiMET_pt",  elescaremet,    {"Electron_pt", "Electron_pt_scare", "Electron_phi", "PuppiMET_pt_ele_unsmear", "PuppiMET_phi_ele_unsmear"})
                   .Redefine("Electron_pt", "Electron_pt_scare");

        if (isDefined("lep4vecs")) {
            _rlm = _rlm.Redefine("lep4vecs", ::gen4vec, {"Electron_pt", "Electron_eta", "Electron_phi", "Electron_mass"});
        }
    }
}

void NanoAODAnalyzerrdframe::calculateElectronSF(string elecFile, string elecYear){
    // Electron SF
    cout << "Loading Electron SF" << endl;
    if (elecFile.empty() || elecYear.empty()) {
        std::cerr << "[ERROR] calculateElectronSF: elecFile ('" << elecFile << "') or elecYear ('" << elecYear << "') is empty!" << std::endl;
        return;
    }
    if (!isDefined("Electron_pt") || !isDefined("Electron_eta") || !isDefined("Electron_deltaEtaSC") || !isDefined("Electron_phi")) {
        std::cerr << "[ERROR] calculateElectronSF: Missing required Electron_pt, Electron_eta, Electron_deltaEtaSC, or Electron_phi branch!" << std::endl;
        return;
    }

    auto elecSFreader = loadCorrectionSet("data/ElectronSF/"+elecFile+"/electron.json.gz");
    auto _elecid = elecSFreader->at("Electron-ID-SF");

    auto elecSFId = [this, _elecid, elecYear](std::string var){
        return [this, _elecid, elecYear, var](const floats &pt, const floats &eta, const floats &deltaEtaSC, const floats &phi)->floats {
            floats wVec;
            wVec.reserve(7); // [0]: nom, [1]: tot_up, [2]: tot_down, [3]: stat_up, [4]: stat_down, [5]: syst_up, [6]: syst_down
            std::string wp = var;
            if (pt.size() == 1){
                if (eta.size() != pt.size() || deltaEtaSC.size() != pt.size() || phi.size() != pt.size()) {
                    std::cerr << "[ERROR] calculateElectronSF: Size mismatch among Electron_pt (" << pt.size()
                              << "), Electron_eta (" << eta.size()
                              << "), Electron_deltaEtaSC (" << deltaEtaSC.size()
                              << "), Electron_phi (" << phi.size() << ")!" << std::endl;
                    return {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
                }
                float scEta = eta[0] + deltaEtaSC[0];
                float sf = 1.0f, sf_up = 1.0f, sf_down = 1.0f;
                float stat_up = 1.0f, stat_down = 1.0f;
                float syst_up = 1.0f, syst_down = 1.0f;
                if (wp.find("Reco") != std::string::npos) {
                    if (pt[0] < 20) wp = "RecoBelow20";
                    else if (pt[0] < 75) wp = "Reco20to75";
                    else wp = "RecoAbove75";
                }
                if (elecYear.find("2022")!=std::string::npos || elecYear.find("2024")!=std::string::npos){
                    sf = _elecid->evaluate({elecYear, "sf", wp, scEta, pt[0]});
                    sf_up = _elecid->evaluate({elecYear, "sfup", wp, scEta, pt[0]});
                    sf_down = _elecid->evaluate({elecYear, "sfdown", wp, scEta, pt[0]});
                    float err_stat = _elecid->evaluate({elecYear, "err_stat", wp, scEta, pt[0]});
                    float err_syst = _elecid->evaluate({elecYear, "err_syst", wp, scEta, pt[0]});
                    stat_up = sf + err_stat;
                    stat_down = sf - err_stat;
                    syst_up = sf + err_syst;
                    syst_down = sf - err_syst;
                } else if (elecYear.find("2023")!=std::string::npos){
                    sf = _elecid->evaluate({elecYear, "sf", wp, scEta, pt[0], phi[0]});
                    sf_up = _elecid->evaluate({elecYear, "sfup", wp, scEta, pt[0], phi[0]});
                    sf_down = _elecid->evaluate({elecYear, "sfdown", wp, scEta, pt[0], phi[0]});
                    float err_stat = _elecid->evaluate({elecYear, "err_stat", wp, scEta, pt[0], phi[0]});
                    float err_syst = _elecid->evaluate({elecYear, "err_syst", wp, scEta, pt[0], phi[0]});
                    stat_up = sf + err_stat;
                    stat_down = sf - err_stat;
                    syst_up = sf + err_syst;
                    syst_down = sf - err_syst;
                } else {
                    std::cerr << "[ERROR] calculateElectronSF: Unrecognized elecYear '" << elecYear << "' for Electron-ID-SF!" << std::endl;
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
        _rlm = _rlm.Redefine("elecWeightReco", elecSFId("RecoAbove75"), {"Electron_pt","Electron_eta","Electron_deltaEtaSC","Electron_phi"})
                   .Redefine("elecWeightId", elecSFId("wp90iso"), {"Electron_pt","Electron_eta","Electron_deltaEtaSC","Electron_phi"});
    } else {
        _rlm = _rlm.Define("elecWeightReco", elecSFId("RecoAbove75"), {"Electron_pt","Electron_eta","Electron_deltaEtaSC","Electron_phi"})
                   .Define("elecWeightId", elecSFId("wp90iso"), {"Electron_pt","Electron_eta","Electron_deltaEtaSC","Electron_phi"});
    }

    auto elecHltSFreader = loadCorrectionSet("data/ElectronSF/"+elecFile+"/electronHlt.json.gz");
    auto _elechlt = elecHltSFreader->at("Electron-HLT-SF");

    auto elecSFTrg = [this, _elechlt, elecYear](const floats &pt, const floats &eta)->floats {
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


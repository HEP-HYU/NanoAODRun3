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

void NanoAODAnalyzerrdframe::calculateMuonSF(string muonid, string muoniso, string muonhlt){
    //// Muon SF
    cout<<"Loading Muon SF"<<endl;
    std::string muonFile = _year;

    if (_isRun22) {
        muonFile = "2022";
    } else if (_isRun22EE) {
        muonFile = "2022_EE";
    } else if (_isRun23) {
        muonFile = "2023";
    } else if (_isRun23BPix) {
        muonFile = "2023_BPix";
    } else if (_isRun24){
        muonFile = "2024";
    }

    auto muonSFreader = loadCorrectionSet("data/MuonSF/ScaleFactors_Muon_Z_ID_ISO_"+muonFile+"_schemaV2.json");
    auto _muonid = muonSFreader->at(muonid);
    auto _muoniso = muonSFreader->at(muoniso);
    if (_isRun24) muonFile = "2023_BPix";
    auto muonHltSFreader = loadCorrectionSet("data/MuonSF/ScaleFactors_Muon_Z_HLT_"+muonFile+"_eta_pt_schemaV2.json");
    auto _muontrg = muonHltSFreader->at(muonhlt);

    // We have only one muon!
    auto muonSFId = [this, _muonid, muonFile](floats &pt, floats &eta)->floats {
        floats wVec;
        wVec.reserve(3); //cent, up, down

        float tmp_eta = -1000;
        if (pt.size() == 1){
            tmp_eta = (muonFile.find("2022") != std::string::npos) ? std::abs(eta[0]) : eta[0];
            float sf = _muonid->evaluate({tmp_eta, pt[0], "nominal"});
            float sf_up = _muonid->evaluate({tmp_eta, pt[0], "systup"});
            float sf_down = _muonid->evaluate({tmp_eta, pt[0], "systdown"});
            wVec.emplace_back(sf);
            wVec.emplace_back(sf_up);
            wVec.emplace_back(sf_down);
        } else wVec = {1.0, 1.0, 1.0};
        return wVec;
    };

    // We have only one muon!
    auto muonSFIso = [this, _muoniso, muonFile](floats &pt, floats &eta)->floats {
        floats wVec;
        wVec.reserve(3); //cent, up, down

        if (pt.size() == 1){
            float tmp_eta = (muonFile.find("2022") != std::string::npos) ? std::abs(eta[0]) : eta[0];
            float sf = _muoniso->evaluate({tmp_eta, pt[0], "nominal"});
            float sf_up = _muoniso->evaluate({tmp_eta, pt[0], "systup"});
            float sf_down = _muoniso->evaluate({tmp_eta, pt[0], "systdown"});
            wVec.emplace_back(sf);
            wVec.emplace_back(sf_up);
            wVec.emplace_back(sf_down);
        } else wVec = {1.0, 1.0, 1.0};
        return wVec;
    };

    auto muonSFTrg = [this, _muontrg](floats &pt, floats &eta)->floats {

        floats wVec;
        wVec.reserve(3); //cent, up, down

        if (pt.size() == 1) {
            float sf = _muontrg->evaluate({eta[0], pt[0], "nominal"});
            float sf_up = _muontrg->evaluate({eta[0], pt[0], "systup"});
            float sf_down = _muontrg->evaluate({eta[0], pt[0], "systdown"});
            wVec.emplace_back(sf);
            wVec.emplace_back(sf_up);
            wVec.emplace_back(sf_down);
        }
        else wVec = {1.0, 1.0, 1.0};
        return wVec;
    };

    _rlm = _rlm.Define("muonWeightId", muonSFId, {"Muon_pt","Muon_eta"})
               .Define("muonWeightIso", muonSFIso, {"Muon_pt","Muon_eta"});
    if (_isRun24){
        _rlm = _rlm.Define("muonWeightTrg", "std::vector<double>{1.0, 1.0, 1.0}");
    }
    else{
        _rlm = _rlm.Define("muonWeightTrg", muonSFTrg, {"Muon_pt","Muon_eta"});
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
                   .Redefine("MET_phi", muonhighscalemetphi, {"Muon_pt", "Muon_pt_scale", "Muon_phi", "MET_pt", "MET_phi"})
                   .Redefine("MET_pt", muonhighscalemet, {"Muon_pt", "Muon_pt_scale", "MET_pt"})
                   .Redefine("Muon_pt", "Muon_pt_scale"); //order matters
    } else if (_syst.find("muonhighscaledown") != std::string::npos) {
        _rlm = _rlm.Define("Muon_pt_scale", muonhighscaledn, {"Muon_pt", "Muon_eta", "Muon_phi", "Muon_charge"})
                   .Redefine("MET_phi", muonhighscalemetphi, {"Muon_pt", "Muon_pt_scale", "Muon_phi", "MET_pt", "MET_phi"})
                   .Redefine("MET_pt", muonhighscalemet, {"Muon_pt", "Muon_pt_scale", "MET_pt"})
                   .Redefine("Muon_pt", "Muon_pt_scale"); //order matters
    }
}
 
void NanoAODAnalyzerrdframe::calculateElectronSF(string elecFile, string elecYear){
    //Electron SF
    cout << "Loading Electron SF" << endl;

    auto elecSFreader = loadCorrectionSet("data/ElectronSF/"+elecFile+"/electron.json.gz");
    auto _elecid = elecSFreader->at("Electron-ID-SF");

    auto elecSFId = [this, _elecid, elecYear](std::string var){
        return [this, _elecid, elecYear, var](floats &pt, floats &eta, floats &phi)->floats {
            floats wVec;
            wVec.reserve(3); //cent, up, down
            std::string wp = var;
            if (pt.size() == 1){
                float sf=1.0; float sf_up=1.0; float sf_down=1.0;
                if (pt[0] >= 20 && pt[0] < 75 && wp.find("Reco")!=std::string::npos) wp = "Reco20to75";
                if (elecYear.find("2022")!=std::string::npos || elecYear.find("2024")!=std::string::npos){
                    sf = _elecid->evaluate({elecYear, "sf", wp, eta[0], pt[0]});
                    sf_up = _elecid->evaluate({elecYear, "sfup", wp, eta[0], pt[0]});
                    sf_down = _elecid->evaluate({elecYear, "sfdown", wp, eta[0], pt[0]});
                } else if (elecYear.find("2023")!=std::string::npos){
                    sf = _elecid->evaluate({elecYear, "sf", wp, eta[0], pt[0], phi[0]});
                    sf_up = _elecid->evaluate({elecYear, "sfup", wp, eta[0], pt[0], phi[0]});
                    sf_down = _elecid->evaluate({elecYear, "sfdown", wp, eta[0], pt[0], phi[0]});
                }
                wVec.emplace_back(sf);
                wVec.emplace_back(sf_up);
                wVec.emplace_back(sf_down);
            } else wVec = {1.0, 1.0, 1.0};
            return wVec;
        };
    };

    _rlm = _rlm.Define("elecWeightReco", elecSFId("RecoAbove75"), {"Electron_pt","Electron_eta","Electron_phi"})
               .Define("elecWeightId", elecSFId("wp90iso"), {"Electron_pt","Electron_eta","Electron_phi"});

    auto elecHltSFreader = loadCorrectionSet("data/ElectronSF/"+elecFile+"/electronHlt.json.gz");
    auto _elechlt = elecHltSFreader->at("Electron-HLT-SF");

    auto elecSFTrg = [this, _elechlt, elecYear](floats &pt, floats &eta)->floats {
        floats wVec;
        wVec.reserve(3); //cent, up, down
        string hltpath = "HLT_SF_Ele30_MVAiso90ID";
        if (pt.size() == 1){
            float sf = _elechlt->evaluate({elecYear, "sf", hltpath, eta[0], pt[0]});
            float sf_up = _elechlt->evaluate({elecYear, "sfup", hltpath, eta[0], pt[0]});
            float sf_down = _elechlt->evaluate({elecYear, "sfdown", hltpath, eta[0], pt[0]});
            wVec.emplace_back(sf);
            wVec.emplace_back(sf_up);
            wVec.emplace_back(sf_down);
        } else wVec = {1.0, 1.0, 1.0};
        return wVec;
    };

    _rlm = _rlm.Define("elecWeightTrg", elecSFTrg, {"Electron_pt","Electron_eta"});
}


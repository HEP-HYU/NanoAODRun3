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
#include <sys/stat.h>
using namespace std;

void NanoAODAnalyzerrdframe::setupJetMETCorrection(string jecFile, string jecYear, string jerMap, string jecVersion, bool dataMc, string jecYearData, string jerVersion){
    auto jercReader = loadCorrectionSet("data/JME/"+jecFile+"/jet_jerc.json.gz");
    string datamcflag = "";
    if (dataMc) datamcflag = "DATA";
    else datamcflag = "MC";

    // DATA JEC keys include a run-era label (e.g. _RunCD_) that MC keys do not.
    // jecYearData should be set to the data-specific year string (e.g. "Summer22_22Sep2023_RunCD").
    // If not provided, fall back to jecYear (preserves MC behavior).
    const string jecYearLookup = (dataMc && !jecYearData.empty()) ? jecYearData : jecYear;

    string tmp = jecYearLookup + jecVersion + datamcflag + "_L1FastJet_AK4PFPuppi";
    std::cout << tmp << std::endl;
    auto _L1FastJet   = jercReader->at(jecYearLookup + jecVersion + datamcflag + "_L1FastJet_AK4PFPuppi");
    auto _L2Relative  = jercReader->at(jecYearLookup + jecVersion + datamcflag + "_L2Relative_AK4PFPuppi");
    auto _L2L3Residual= jercReader->at(jecYearLookup + jecVersion + datamcflag + "_L2L3Residual_AK4PFPuppi");

    auto applyJes = [this, _L1FastJet, _L2Relative, _L2L3Residual, dataMc](floats jetpts, floats jetetas, floats jetphis, floats jetAreas, floats jetrawf, float rho, unsigned int run, floats toCorr)->floats {
        floats corrfactors;
        corrfactors.reserve(jetpts.size());

        for (size_t i=0; i<jetpts.size(); i++) {
            float rawFact = (1.0-jetrawf[i]);
            float rawjetpt = jetpts[i] * rawFact;
            float L1corr = _L1FastJet->evaluate({jetAreas[i], jetetas[i], rawjetpt, rho});
            float L1pt = rawjetpt * L1corr;
            float L2corr = 1; float L2L3corr = 1;
            // L2Relative: Run3Summer24 AND Run3Summer23BPix add JetPhi as 2nd input → [JetEta, JetPhi, JetPt].
            // Earlier eras (22, 22EE, 23 preBPix) use only [JetEta, JetPt].
            // Passing the wrong number of args corrupts the heap.
            if (_isRun24 || _isRun23BPix)
                L2corr = _L2Relative->evaluate({jetetas[i], jetphis[i], L1pt});
            else
                L2corr = _L2Relative->evaluate({jetetas[i], L1pt});
            float L2pt = L1pt * L2corr;
            // L2L3Residual: Run3Summer24 stores the run number as an additional first input → [run, JetEta, JetPt].
            // All pre-24 eras use only [JetEta, JetPt].
            // Passing 3 args to a 2-input correction corrupts the heap → munmap_chunk for 22/23 DATA.
            if (dataMc) {
                if (_isRun24)
                    L2L3corr = _L2L3Residual->evaluate({float(run), jetetas[i], L2pt});
                else
                    L2L3corr = _L2L3Residual->evaluate({jetetas[i], L2pt});
            }
            float corrfactor = toCorr[i] * rawFact * L1corr * L2corr * L2L3corr; 
            corrfactors.emplace_back(corrfactor);
        }
        return corrfactors;
    };
    _rlm = _rlm.Define("Jet_pt_uncorr", "Jet_pt");
    _rlm = _rlm.Define("Jet_pt_corr", applyJes, {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_area", "Jet_rawFactor", "Rho_fixedGridRhoFastjetAll", "run", "Jet_pt"})
               .Redefine("Jet_mass", applyJes, {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_area", "Jet_rawFactor", "Rho_fixedGridRhoFastjetAll", "run", "Jet_mass"});
    _rlm = _rlm.Redefine("Jet_pt", "Jet_pt_corr");

    // JER: only for MC. jercReader is a local unique_ptr; _jerReso/_jerSF are raw pointers
    // into jercReader's memory. On data this block is skipped so the pointers never dangle.
    if (!_isData) {
        auto _jerReso = jercReader->at(jerMap + "_" + jerVersion + "_MC_PtResolution_AK4PFPuppi");
        auto _jerSF   = jercReader->at(jerMap + "_" + jerVersion + "_MC_ScaleFactor_AK4PFPuppi");
        auto jerReader = loadCorrectionSet("data/JME/jer_smear.json.gz");
        auto _jerSmear = jerReader->at("JERSmear");

        auto applyJer = [this, _jerReso, _jerSF, _jerSmear](floats jetpts, floats jetetas, floats jetphis, shorts genidxs, floats genpts, floats genetas, floats genphis, float rho, ULong64_t event, floats toCorr)->floats {
            floats corrfactors;
            corrfactors.reserve(jetpts.size());

            for (size_t i=0; i<jetpts.size(); i++){
                float reso = _jerReso->evaluate({jetetas[i], jetpts[i], rho});
                // Run3Summer24 ScaleFactor has no systematic input: [JetEta, JetPt] (2 inputs).
                // Earlier eras have [JetEta, JetPt, systematic] (3 inputs).
                // Passing 3 args to a 2-input correction corrupts the heap.
                float sf = _isRun24 ? _jerSF->evaluate({jetetas[i], jetpts[i]})
                                    : _jerSF->evaluate({jetetas[i], jetpts[i], "nom"});
                float genPtForSmear = -1.0;

                int genidx = genidxs[i];
                if (genidx >= 0){
                    float dPhi = std::abs(genphis[genidx] - jetphis[i]);
                    float dEta = std::abs(genetas[genidx] - jetetas[i]);
                    float dR2 = dPhi*dPhi + dEta*dEta;
                    if (dR2 < 0.04 && std::abs(jetpts[i]-genpts[genidx]) < (3*reso*jetpts[i])) {
                        genPtForSmear = genpts[genidx];
                    }
                }
                float smear = _jerSmear->evaluate({jetpts[i], jetetas[i], genPtForSmear, rho, int(event), reso, sf});
                float corr = (std::isfinite(smear) && smear > 0.0) ? smear : 1.0;
                float corrfactor = toCorr[i] * corr;
                corrfactors.emplace_back(corrfactor);
            }
            return corrfactors;
        };
        _rlm = _rlm.Redefine("Jet_pt_corr", applyJer, {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_genJetIdx", "GenJet_pt", "GenJet_eta", "GenJet_phi", "Rho_fixedGridRhoFastjetAll", "event", "Jet_pt"})
                   .Redefine("Jet_mass",    applyJer, {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_genJetIdx", "GenJet_pt", "GenJet_eta", "GenJet_phi", "Rho_fixedGridRhoFastjetAll", "event", "Jet_mass"})
                   .Redefine("Jet_pt", "Jet_pt_corr");
    }
}

void NanoAODAnalyzerrdframe::JetVetoMap(string jetFile, string map) {
    cout << "apply JetVetoMap" << endl;
    

    //TODO: _isRun24 is implementation for NanoAODv15
    if (_isRun24){
        auto jetIDreader = loadCorrectionSet("data/JME/2024_Summer24/jetid.json.gz");
        auto _jetid = jetIDreader->at("AK4PUPPI_TightLeptonVeto");
        // NanoAODv15: Jet_chMultiplicity and Jet_neMultiplicity are UChar_t
        auto getJetId1 = [this, _jetid](floats &etas, floats &chHEF, floats &neHEF, floats &chEmEF, floats &neEmEF, floats &muEF, uchars &Nch, uchars &Nne)->floats{
            floats jetIds;
            jetIds.reserve(etas.size());
            for (size_t i=0; i<etas.size(); i++){
                float jetId = _jetid->evaluate({float(etas[i]), float(chHEF[i]), float(neHEF[i]), float(chEmEF[i]), float(neEmEF[i]), float(muEF[i]), int(Nch[i]), int(Nne[i]), int(Nch[i])+int(Nne[i])});
                jetIds.emplace_back(jetId);
            }
            return jetIds;
        };
        _rlm = _rlm.Define("Jet_passJetIdTightLepVeto", getJetId1, {"Jet_eta", "Jet_chHEF", "Jet_neHEF", "Jet_chEmEF", "Jet_neEmEF", "Jet_muEF", "Jet_chMultiplicity", "Jet_neMultiplicity"});
    } 
    //TODO: below is for NanoAODv12
    //If the others are going to v15, this should be removed
    else{
        // NanoAODv12: Jet_jetId is UChar_t (bit-field stored as unsigned char)
        auto getJetId2 = [this](floats &etas, uchars &jetIds, floats &muEF, floats &chEmEF, floats &neHEF, floats &neEmEF)->floats{
            floats idVec;
            idVec.reserve(etas.size());
            for (size_t i=0; i<etas.size(); i++){
                float id=0;
                if (std::abs(etas[i]) <= 2.7) id  = (jetIds[i] & (1 << 1)) && (muEF[i] < 0.8) && (chEmEF[i] < 0.8) ;
                else if (std::abs(etas[i]) > 2.7 && std::abs(etas[i]) <= 3.0) id = (jetIds[i] & (1 << 1)) && (neHEF[i] < 0.99);
                else if (std::abs(etas[i]) > 3.0) id = (jetIds[i] & (1 << 1)) && (neEmEF[i] < 0.4);
                idVec.emplace_back(id);
            }
            return idVec;
        };
        _rlm = _rlm.Define("Jet_passJetIdTightLepVeto", getJetId2, {"Jet_eta", "Jet_jetId", "Jet_muEF", "Jet_chEmEF", "Jet_neHEF", "Jet_neEmEF"});


    }

    
    _rlm = _rlm.Define("loosejetcuts", "Jet_pt>15 && (Jet_passJetIdTightLepVeto==1.0) && (Jet_neEmEF+Jet_chEmEF)<0.9");
    
    _rlm = _rlm.Define("Jet_pt_loosejet", "Jet_pt[loosejetcuts]")
               .Define("Jet_phi_loosejet", "Jet_phi[loosejetcuts]")
               .Define("Jet_eta_loosejet", "Jet_eta[loosejetcuts]")
               .Define("Jet_mass_loosejet", "Jet_mass[loosejetcuts]")
               .Define("njetloose", "int(Jet_pt_loosejet.size())")
               .Define("loosejet4vecs", ::gen4vec, {"Jet_pt_loosejet", "Jet_eta_loosejet", "Jet_phi_loosejet", "Jet_mass_loosejet"});
    
    auto vetoMapreader = loadCorrectionSet("data/JME/"+jetFile+"/jetvetomaps.json.gz");
    
    auto _vetomap = vetoMapreader->at(map);
 
    auto vetomap = [_vetomap](const std::string &key) {
        return [_vetomap, key](floats &eta, floats &phi)->floats {
            floats xout;
            xout.reserve(eta.size());
            for (size_t i=0; i<eta.size(); i++) {
                xout.emplace_back(_vetomap->evaluate({key, eta[i], phi[i]}));
                //xout.emplace_back(0.0);
            }
            return xout;
        };
    };
    
    _rlm = _rlm.Define("Jet_isVeto_loose", vetomap("jetvetomap"), {"Jet_eta_loosejet","Jet_phi_loosejet"})
               .Define("events_isVeto","Sum(Jet_isVeto_loose)");

    if (_isRun24){
        _rlm = _rlm.Define("Jet_isVeto_fpix", vetomap("jetvetomap_fpix"), {"Jet_eta_loosejet", "Jet_phi_loosejet"})
                   .Redefine("events_isVeto", "Sum(Jet_isVeto_loose)+Sum(Jet_isVeto_fpix)");
    }
}

// Apply MET XY-shift corrections using correctionlib.
// The met_xyCorrections JSON provides linear corrections to PuppiMET (and PFMET)
// to remove the azimuthal modulation caused by detector non-uniformities.
// Inputs: met_pt, met_phi, npvGood, data/MC flag, era epoch.
// The corrected pt and phi are obtained by evaluating the correction twice
// (once for 'pt', once for 'phi'), then redefining the MET branches.
// NOTE: 2024_Summer24 has no met_xyCorrections file yet — correction is skipped.
void NanoAODAnalyzerrdframe::applyMETXYCorrection(const std::string &jecFile, const std::string &epoch) {
    // Only 2022 and 2023 eras have correction files available.
    const std::string metXYFile = "data/JME/" + jecFile + "/met_xyCorrections_" +
                                   epoch.substr(0, 4) + "_" + epoch + ".json.gz";

    // Check file exists before loading (2024 has no correction file).
    struct stat buf;
    if (::stat(metXYFile.c_str(), &buf) != 0) {
        std::cout << "[MET XY] No correction file for " << jecFile
                  << " (‘" << metXYFile << "' not found). Skipping XY correction." << std::endl;
        return;
    }

    std::cout << "[MET XY] Applying XY correction from " << metXYFile << std::endl;
    auto metXYreader = loadCorrectionSet(metXYFile);
    auto _metXY = metXYreader->at("met_xy_corrections");

    const std::string dtmc    = _isData ? "DATA" : "MC";
    const std::string ep      = epoch;   // e.g. "2022", "2022EE", "2023", "2023BPix"

    // Lambda: evaluate corrected PuppiMET pt/phi.
    // The correction returns the corrected value directly (pt or phi).
    // PV_npvsGood is UChar_t (unsigned char) in NanoAOD.
    auto applyXY = [this, _metXY, dtmc, ep](float met_pt, float met_phi, unsigned char npvGood) -> std::pair<float,float> {
        float npv = float(npvGood);
        float new_pt  = _metXY->evaluate({"pt",  "PuppiMET", ep, dtmc, "nom", met_pt, met_phi, npv});
        float new_phi = _metXY->evaluate({"phi", "PuppiMET", ep, dtmc, "nom", met_pt, met_phi, npv});
        return {new_pt, new_phi};
    };

    auto getCorr_pt = [applyXY](float met_pt, float met_phi, unsigned char npvGood) -> float {
        return applyXY(met_pt, met_phi, npvGood).first;
    };
    auto getCorr_phi = [applyXY](float met_pt, float met_phi, unsigned char npvGood) -> float {
        return applyXY(met_pt, met_phi, npvGood).second;
    };

    // Redefine PuppiMET branches with XY-corrected values.
    // MET_phi must be redefined before MET_pt because phi depends on old pt.
    _rlm = _rlm
        .Define("PuppiMET_phi_xyCorr", getCorr_phi, {"PuppiMET_pt", "PuppiMET_phi", "PV_npvsGood"})
        .Define("PuppiMET_pt_xyCorr",  getCorr_pt,  {"PuppiMET_pt", "PuppiMET_phi", "PV_npvsGood"})
        .Redefine("PuppiMET_phi", "PuppiMET_phi_xyCorr")
        .Redefine("PuppiMET_pt",  "PuppiMET_pt_xyCorr");

    // Also apply to PFMET for completeness.
    auto getCorr_pfpt = [this, _metXY, dtmc, ep](float met_pt, float met_phi, unsigned char npvGood) -> float {
        float npv = float(npvGood);
        return _metXY->evaluate({"pt",  "MET", ep, dtmc, "nom", met_pt, met_phi, npv});
    };
    auto getCorr_pfphi = [this, _metXY, dtmc, ep](float met_pt, float met_phi, unsigned char npvGood) -> float {
        float npv = float(npvGood);
        return _metXY->evaluate({"phi", "MET", ep, dtmc, "nom", met_pt, met_phi, npv});
    };
    _rlm = _rlm
        .Define("MET_phi_xyCorr", getCorr_pfphi, {"MET_pt", "MET_phi", "PV_npvsGood"})
        .Define("MET_pt_xyCorr",  getCorr_pfpt,  {"MET_pt", "MET_phi", "PV_npvsGood"})
        .Redefine("MET_phi", "MET_phi_xyCorr")
        .Redefine("MET_pt",  "MET_pt_xyCorr");
}

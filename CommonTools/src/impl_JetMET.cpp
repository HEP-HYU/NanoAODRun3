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
            // L2L3Residual: DATA keys in ALL Run3 eras require run as 1st input → [run, JetEta, JetPt].
            if (dataMc) {
                L2L3corr = _L2L3Residual->evaluate({float(run), jetetas[i], L2pt});
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
                // JER ScaleFactor across ALL Run3 eras (2022, 2022EE, 2023, 2023BPix, 2024)
                // takes 2 inputs: [JetEta, JetPt].
                float sf = _jerSF->evaluate({jetetas[i], jetpts[i]});
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

    // -------------------------------------------------------------------------
    // [JME-01] Type-1 MET Recomputation (CMS_jmeType1METCorrection.pdf & CMS_jmeFAQ.pdf)
    // Propagate nominal JEC (and for MC, nominal JER smearing) to PuppiMET (and MET).
    // Baseline: RawPuppiMET_pt, RawPuppiMET_phi.
    // Jets: standard AK4 jets (Jet_) and low-pT MET jets (CorrT1METJet_).
    // Muon-subtracted raw pT used for L1FastJet; Delta pT = pT_full - pT_L1.
    // -------------------------------------------------------------------------

    // 1. Check mandatory EM fraction branches
    if (!isDefined("Jet_chEmEF") || !isDefined("Jet_neEmEF")) {
        std::cerr << "[ERROR] setupJetMETCorrection: 'Jet_chEmEF' or 'Jet_neEmEF' branch not found! Cannot evaluate EM fraction < 0.9 for Type-1 MET." << std::endl;
        throw std::runtime_error("setupJetMETCorrection: missing Jet EM fraction branches");
    }

    // 2. Safe wrappers for Jet muon subtraction branches
    const bool hasJetMuSub = isDefined("Jet_muonSubtrFactor");
    if (!hasJetMuSub) {
        std::cerr << "[ERROR] setupJetMETCorrection: 'Jet_muonSubtrFactor' branch not found! Falling back to raw jet pT without muon subtraction." << std::endl;
        _rlm = _rlm.Define("Jet_muonSubtrFactor_safe", [](const floats &pts) { return floats(pts.size(), 0.0f); }, {"Jet_pt"});
    } else {
        _rlm = _rlm.Define("Jet_muonSubtrFactor_safe", "Jet_muonSubtrFactor");
    }

    const bool hasJetMuDeltaPhi = isDefined("Jet_muonSubtrDeltaPhi");
    if (!hasJetMuDeltaPhi) {
        // Per CMS_jmeType1METCorrection.pdf Section 1.5 note: fallback to phi_jet
        _rlm = _rlm.Define("Jet_muonSubtrDeltaPhi_safe", [](const floats &pts) { return floats(pts.size(), 0.0f); }, {"Jet_pt"});
    } else {
        _rlm = _rlm.Define("Jet_muonSubtrDeltaPhi_safe", "Jet_muonSubtrDeltaPhi");
    }

    // 3. Safe wrappers for CorrT1METJet branches (low-pT MET-specific jets)
    const bool hasCorrT1 = isDefined("CorrT1METJet_rawPt") && isDefined("CorrT1METJet_eta") &&
                          isDefined("CorrT1METJet_phi") && isDefined("CorrT1METJet_area");
    if (!hasCorrT1) {
        std::cerr << "[ERROR] setupJetMETCorrection: 'CorrT1METJet' branches not found! Skipping low-pT CorrT1METJet contribution to Type-1 MET." << std::endl;
        _rlm = _rlm.Define("CorrT1METJet_rawPt_safe",             []() { return floats{}; })
                   .Define("CorrT1METJet_eta_safe",               []() { return floats{}; })
                   .Define("CorrT1METJet_phi_safe",               []() { return floats{}; })
                   .Define("CorrT1METJet_area_safe",              []() { return floats{}; })
                   .Define("CorrT1METJet_muonSubtrFactor_safe",   []() { return floats{}; })
                   .Define("CorrT1METJet_muonSubtrDeltaPhi_safe", []() { return floats{}; })
                   .Define("CorrT1METJet_EmEF_safe",              []() { return floats{}; });
    } else {
        _rlm = _rlm.Define("CorrT1METJet_rawPt_safe", "CorrT1METJet_rawPt")
                   .Define("CorrT1METJet_eta_safe",   "CorrT1METJet_eta")
                   .Define("CorrT1METJet_phi_safe",   "CorrT1METJet_phi")
                   .Define("CorrT1METJet_area_safe",  "CorrT1METJet_area");

        if (!isDefined("CorrT1METJet_muonSubtrFactor")) {
            std::cerr << "[ERROR] setupJetMETCorrection: 'CorrT1METJet_muonSubtrFactor' not found! Falling back to rawPt for CorrT1METJet." << std::endl;
            _rlm = _rlm.Define("CorrT1METJet_muonSubtrFactor_safe", [](const floats &pts) { return floats(pts.size(), 0.0f); }, {"CorrT1METJet_rawPt"});
        } else {
            _rlm = _rlm.Define("CorrT1METJet_muonSubtrFactor_safe", "CorrT1METJet_muonSubtrFactor");
        }

        if (!isDefined("CorrT1METJet_muonSubtrDeltaPhi")) {
            _rlm = _rlm.Define("CorrT1METJet_muonSubtrDeltaPhi_safe", [](const floats &pts) { return floats(pts.size(), 0.0f); }, {"CorrT1METJet_rawPt"});
        } else {
            _rlm = _rlm.Define("CorrT1METJet_muonSubtrDeltaPhi_safe", "CorrT1METJet_muonSubtrDeltaPhi");
        }

        if (!isDefined("CorrT1METJet_EmEF")) {
            // Per CMS_jmeType1METCorrection.pdf Section 1.4: relax EM fraction cut if branch is absent
            _rlm = _rlm.Define("CorrT1METJet_EmEF_safe", [](const floats &pts) { return floats(pts.size(), 0.0f); }, {"CorrT1METJet_rawPt"});
        } else {
            _rlm = _rlm.Define("CorrT1METJet_EmEF_safe", "CorrT1METJet_EmEF");
        }
    }

    // 4. Safe baseline for PuppiMET
    const bool hasRawPuppiMET = isDefined("RawPuppiMET_pt") && isDefined("RawPuppiMET_phi");
    if (!hasRawPuppiMET) {
        std::cerr << "[ERROR] setupJetMETCorrection: 'RawPuppiMET_pt' or 'RawPuppiMET_phi' not found! Falling back to delta propagation from PuppiMET_pt/phi." << std::endl;
        _rlm = _rlm.Define("RawPuppiMET_pt_safe",  "PuppiMET_pt")
                   .Define("RawPuppiMET_phi_safe", "PuppiMET_phi")
                   .Define("rawPuppi_flag",        []() { return 0.0f; });
    } else {
        _rlm = _rlm.Define("RawPuppiMET_pt_safe",  "RawPuppiMET_pt")
                   .Define("RawPuppiMET_phi_safe", "RawPuppiMET_phi")
                   .Define("rawPuppi_flag",        []() { return 1.0f; });
    }

    if (isDefined("RawPuppiMET_sumEt")) {
        _rlm = _rlm.Define("RawPuppiMET_sumEt_safe", "RawPuppiMET_sumEt");
    } else if (isDefined("PuppiMET_sumEt")) {
        _rlm = _rlm.Define("RawPuppiMET_sumEt_safe", "PuppiMET_sumEt");
    } else {
        _rlm = _rlm.Define("RawPuppiMET_sumEt_safe", []() { return 0.0f; });
    }

    // 5. Unified Type-1 MET calculation lambda
    auto calcType1MET = [this, _L1FastJet, _L2Relative, _L2L3Residual, dataMc](
        float rawmet_pt, float rawmet_phi, float rawmet_sumEt, float raw_flag,
        const floats &jetpts_uncorr, const floats &jetpts_full,
        const floats &jetetas, const floats &jetphis, const floats &jetareas,
        const floats &jetrawf, const floats &jetmuonf, const floats &jetmuondphi,
        const floats &chEmEF, const floats &neEmEF,
        const floats &corrT1_rawPt, const floats &corrT1_eta, const floats &corrT1_phi,
        const floats &corrT1_area, const floats &corrT1_muonf, const floats &corrT1_muondphi,
        const floats &corrT1_EmEF,
        float rho, unsigned int run) -> floats {

        if (jetpts_uncorr.size() != jetpts_full.size() || jetpts_uncorr.size() != jetetas.size() ||
            jetpts_uncorr.size() != jetphis.size() || jetpts_uncorr.size() != jetareas.size() ||
            jetpts_uncorr.size() != jetrawf.size() || jetpts_uncorr.size() != jetmuonf.size() ||
            jetpts_uncorr.size() != jetmuondphi.size() || jetpts_uncorr.size() != chEmEF.size() ||
            jetpts_uncorr.size() != neEmEF.size()) {
            std::cerr << "[ERROR] setupJetMETCorrection: Jet vector size mismatch in calcType1MET!" << std::endl;
            throw std::runtime_error("setupJetMETCorrection: Jet vector size mismatch");
        }

        const bool isRaw = (raw_flag > 0.5f);
        double metx = static_cast<double>(rawmet_pt) * std::cos(static_cast<double>(rawmet_phi));
        double mety = static_cast<double>(rawmet_pt) * std::sin(static_cast<double>(rawmet_phi));
        double sumEt_shift = 0.0;

        // 5a. Loop over standard AK4 jets
        for (size_t i = 0; i < jetpts_uncorr.size(); i++) {
            float rawFact = (1.0f - jetrawf[i]);
            float rawjetpt = jetpts_uncorr[i] * rawFact;
            float rawjetpt_musub = rawjetpt * (1.0f - jetmuonf[i]);
            float L1corr = _L1FastJet->evaluate({jetareas[i], jetetas[i], rawjetpt_musub, rho});
            float pT_L1 = rawjetpt_musub * L1corr;
            float pT_full = jetpts_full[i];

            float emFraction = chEmEF[i] + neEmEF[i];
            // Selection per CMS_jmeType1METCorrection.pdf Section 1.4: pT_full > 15 GeV, |eta| < 5.2, EM fraction < 0.9
            if (pT_full > 15.0f && std::abs(jetetas[i]) < 5.2f && emFraction < 0.9f) {
                float dpt = isRaw ? (pT_full - pT_L1) : (pT_full - jetpts_uncorr[i]);
                float phi = jetphis[i] + jetmuondphi[i];
                metx -= static_cast<double>(dpt) * std::cos(static_cast<double>(phi));
                mety -= static_cast<double>(dpt) * std::sin(static_cast<double>(phi));
                sumEt_shift += static_cast<double>(dpt);
            }
        }

        // 5b. Loop over low-pT MET-specific jets (CorrT1METJet), if present
        for (size_t j = 0; j < corrT1_rawPt.size(); j++) {
            float rawPt = corrT1_rawPt[j];
            float rawPt_musub = rawPt * (1.0f - corrT1_muonf[j]);
            float L1corr = _L1FastJet->evaluate({corrT1_area[j], corrT1_eta[j], rawPt_musub, rho});
            float pT_L1 = rawPt_musub * L1corr;

            float L2corr = 1.0f;
            if (_isRun24 || _isRun23BPix)
                L2corr = _L2Relative->evaluate({corrT1_eta[j], corrT1_phi[j], pT_L1});
            else
                L2corr = _L2Relative->evaluate({corrT1_eta[j], pT_L1});
            float pT_L2 = pT_L1 * L2corr;

            float L2L3corr = 1.0f;
            if (dataMc) {
                L2L3corr = _L2L3Residual->evaluate({float(run), corrT1_eta[j], pT_L2});
            }
            float pT_full = pT_L2 * L2L3corr;

            float emFraction = corrT1_EmEF[j];
            if (pT_full > 15.0f && std::abs(corrT1_eta[j]) < 5.2f && emFraction < 0.9f) {
                float dpt = pT_full - pT_L1;
                float phi = corrT1_phi[j] + corrT1_muondphi[j];
                metx -= static_cast<double>(dpt) * std::cos(static_cast<double>(phi));
                mety -= static_cast<double>(dpt) * std::sin(static_cast<double>(phi));
                sumEt_shift += static_cast<double>(dpt);
            }
        }

        float final_pt    = static_cast<float>(std::sqrt(metx * metx + mety * mety));
        float final_phi   = static_cast<float>(std::atan2(mety, metx));
        float final_sumEt = static_cast<float>(std::max(0.0, static_cast<double>(rawmet_sumEt) + sumEt_shift));

        return floats{final_pt, final_phi, final_sumEt};
    };

    // 6. Redefine PuppiMET branches
    _rlm = _rlm.Define("Type1_PuppiMET_vec", calcType1MET,
                       {"RawPuppiMET_pt_safe", "RawPuppiMET_phi_safe", "RawPuppiMET_sumEt_safe", "rawPuppi_flag",
                        "Jet_pt_uncorr", "Jet_pt",
                        "Jet_eta", "Jet_phi", "Jet_area",
                        "Jet_rawFactor", "Jet_muonSubtrFactor_safe", "Jet_muonSubtrDeltaPhi_safe",
                        "Jet_chEmEF", "Jet_neEmEF",
                        "CorrT1METJet_rawPt_safe", "CorrT1METJet_eta_safe", "CorrT1METJet_phi_safe",
                        "CorrT1METJet_area_safe", "CorrT1METJet_muonSubtrFactor_safe", "CorrT1METJet_muonSubtrDeltaPhi_safe",
                        "CorrT1METJet_EmEF_safe",
                        "Rho_fixedGridRhoFastjetAll", "run"})
               .Define("PuppiMET_pt_t1",    [](const floats &v) { return v[0]; }, {"Type1_PuppiMET_vec"})
               .Define("PuppiMET_phi_t1",   [](const floats &v) { return v[1]; }, {"Type1_PuppiMET_vec"})
               .Redefine("PuppiMET_pt",     "PuppiMET_pt_t1")
               .Redefine("PuppiMET_phi",    "PuppiMET_phi_t1");

    if (isDefined("PuppiMET_sumEt")) {
        _rlm = _rlm.Define("PuppiMET_sumEt_t1", [](const floats &v) { return v[2]; }, {"Type1_PuppiMET_vec"})
                   .Redefine("PuppiMET_sumEt",  "PuppiMET_sumEt_t1");
    }

    // 7. Also update PFMET (MET_pt, MET_phi, MET_sumEt) if defined
    if (isDefined("MET_pt") && isDefined("MET_phi")) {
        const bool hasRawPFMET = isDefined("RawMET_pt") && isDefined("RawMET_phi");
        if (hasRawPFMET) {
            _rlm = _rlm.Define("RawMET_pt_safe",  "RawMET_pt")
                       .Define("RawMET_phi_safe", "RawMET_phi")
                       .Define("rawPF_flag",      []() { return 1.0f; });
        } else {
            std::cerr << "[ERROR] setupJetMETCorrection: 'RawMET_pt' or 'RawMET_phi' not found! Falling back to delta propagation from MET_pt/phi." << std::endl;
            _rlm = _rlm.Define("RawMET_pt_safe",  "MET_pt")
                       .Define("RawMET_phi_safe", "MET_phi")
                       .Define("rawPF_flag",      []() { return 0.0f; });
        }
        if (isDefined("RawMET_sumEt")) {
            _rlm = _rlm.Define("RawMET_sumEt_safe", "RawMET_sumEt");
        } else if (isDefined("MET_sumEt")) {
            _rlm = _rlm.Define("RawMET_sumEt_safe", "MET_sumEt");
        } else {
            _rlm = _rlm.Define("RawMET_sumEt_safe", []() { return 0.0f; });
        }

        _rlm = _rlm.Define("Type1_PFMET_vec", calcType1MET,
                           {"RawMET_pt_safe", "RawMET_phi_safe", "RawMET_sumEt_safe", "rawPF_flag",
                            "Jet_pt_uncorr", "Jet_pt",
                            "Jet_eta", "Jet_phi", "Jet_area",
                            "Jet_rawFactor", "Jet_muonSubtrFactor_safe", "Jet_muonSubtrDeltaPhi_safe",
                            "Jet_chEmEF", "Jet_neEmEF",
                            "CorrT1METJet_rawPt_safe", "CorrT1METJet_eta_safe", "CorrT1METJet_phi_safe",
                            "CorrT1METJet_area_safe", "CorrT1METJet_muonSubtrFactor_safe", "CorrT1METJet_muonSubtrDeltaPhi_safe",
                            "CorrT1METJet_EmEF_safe",
                            "Rho_fixedGridRhoFastjetAll", "run"})
                   .Define("MET_pt_t1",    [](const floats &v) { return v[0]; }, {"Type1_PFMET_vec"})
                   .Define("MET_phi_t1",   [](const floats &v) { return v[1]; }, {"Type1_PFMET_vec"})
                   .Redefine("MET_pt",     "MET_pt_t1")
                   .Redefine("MET_phi",    "MET_phi_t1");

        if (isDefined("MET_sumEt")) {
            _rlm = _rlm.Define("MET_sumEt_t1", [](const floats &v) { return v[2]; }, {"Type1_PFMET_vec"})
                       .Redefine("MET_sumEt",  "MET_sumEt_t1");
        }
    }
}

void NanoAODAnalyzerrdframe::JetVetoMap(string jetFile, string map) {
    cout << "apply JetVetoMap" << endl;
    

    // -------------------------------------------------------------------------
    // [JME-04] Jet ID Implementation (CMS_jmeJetID.pdf)
    // - v15 (2024, _isRun24): Jet_jetId branch was removed. Evaluate AK4PUPPI_TightLeptonVeto
    //   from official jetid.json.gz using input constituent variables (CMS_jmeJetID.pdf Page 1 & 2).
    // - v12 (2022-2023, !_isRun24): Jet_jetId bit 1 (tightBit) repaired for forward regions
    //   and lepton veto per CMS_jmeJetID.pdf Page 1 & 3 official C++ recipe.
    // -------------------------------------------------------------------------
    if (_isRun24) {
        const std::vector<std::string> reqV15 = {
            "Jet_eta", "Jet_chHEF", "Jet_neHEF", "Jet_chEmEF", "Jet_neEmEF", "Jet_muEF",
            "Jet_chMultiplicity", "Jet_neMultiplicity"
        };
        for (const auto &b : reqV15) {
            if (!isDefined(b)) {
                std::cerr << "[ERROR] JetVetoMap: Mandatory NanoAODv15 branch '" << b
                          << "' not found! Cannot evaluate AK4PUPPI_TightLeptonVeto." << std::endl;
                throw std::runtime_error("JetVetoMap: missing branch " + b);
            }
        }

        const std::string jsonPath = "data/JME/" + jetFile + "/jetid.json.gz";
        std::shared_ptr<correction::CorrectionSet> jetIDreader;
        try {
            jetIDreader = loadCorrectionSet(jsonPath);
        } catch (const std::exception &e) {
            std::cerr << "[ERROR] JetVetoMap: Failed to load jetid JSON '" << jsonPath << "': " << e.what() << std::endl;
            throw;
        }

        correction::Correction::Ref _jetid;
        try {
            _jetid = jetIDreader->at("AK4PUPPI_TightLeptonVeto");
        } catch (const std::exception &e) {
            std::cerr << "[ERROR] JetVetoMap: Key 'AK4PUPPI_TightLeptonVeto' not found in '" << jsonPath << "': " << e.what() << std::endl;
            throw;
        }

        // NanoAODv15: Jet_chMultiplicity and Jet_neMultiplicity are UChar_t
        auto getJetId1 = [this, _jetid](floats &etas, floats &chHEF, floats &neHEF, floats &chEmEF, floats &neEmEF, floats &muEF, uchars &Nch, uchars &Nne) -> floats {
            floats jetIds;
            jetIds.reserve(etas.size());
            for (size_t i = 0; i < etas.size(); i++) {
                float jetId = 0.0f;
                try {
                    jetId = _jetid->evaluate({float(etas[i]), float(chHEF[i]), float(neHEF[i]), float(chEmEF[i]), float(neEmEF[i]), float(muEF[i]), int(Nch[i]), int(Nne[i]), int(Nch[i]) + int(Nne[i])});
                } catch (const std::exception &e) {
                    std::cerr << "[ERROR] JetVetoMap: Correctionlib evaluation failed for AK4PUPPI_TightLeptonVeto (jet " << i << "): " << e.what() << std::endl;
                    throw;
                }
                jetIds.emplace_back(jetId);
            }
            return jetIds;
        };
        _rlm = _rlm.Define("Jet_passJetIdTightLepVeto", getJetId1,
                           {"Jet_eta", "Jet_chHEF", "Jet_neHEF", "Jet_chEmEF", "Jet_neEmEF", "Jet_muEF", "Jet_chMultiplicity", "Jet_neMultiplicity"});
    } else {
        // NanoAODv12: Official CMS recipe per CMS_jmeJetID.pdf Page 1 & 3
        const std::vector<std::string> reqV12 = {
            "Jet_eta", "Jet_jetId", "Jet_muEF", "Jet_chEmEF", "Jet_neHEF", "Jet_neEmEF"
        };
        for (const auto &b : reqV12) {
            if (!isDefined(b)) {
                std::cerr << "[ERROR] JetVetoMap: Mandatory NanoAODv12 branch '" << b
                          << "' not found! Cannot evaluate TightLeptonVeto Jet ID." << std::endl;
                throw std::runtime_error("JetVetoMap: missing branch " + b);
            }
        }

        // NanoAODv12: Jet_jetId is UChar_t (bit-field stored as unsigned char)
        auto getJetId2 = [this](floats &etas, uchars &jetIds, floats &muEF, floats &chEmEF, floats &neHEF, floats &neEmEF) -> floats {
            floats idVec;
            idVec.reserve(etas.size());
            for (size_t i = 0; i < etas.size(); i++) {
                float id = 0.0f;
                if (std::abs(etas[i]) <= 2.7) id = (jetIds[i] & (1 << 1)) && (muEF[i] < 0.8) && (chEmEF[i] < 0.8);
                else if (std::abs(etas[i]) > 2.7 && std::abs(etas[i]) <= 3.0) id = (jetIds[i] & (1 << 1)) && (neHEF[i] < 0.99);
                else if (std::abs(etas[i]) > 3.0) id = (jetIds[i] & (1 << 1)) && (neEmEF[i] < 0.4);
                idVec.emplace_back(id);
            }
            return idVec;
        };
        _rlm = _rlm.Define("Jet_passJetIdTightLepVeto", getJetId2,
                           {"Jet_eta", "Jet_jetId", "Jet_muEF", "Jet_chEmEF", "Jet_neHEF", "Jet_neEmEF"});
    }

    _rlm = _rlm.Define("loosejetcuts", "Jet_pt>15 && (Jet_passJetIdTightLepVeto==1.0) && (Jet_neEmEF+Jet_chEmEF)<0.9");

    _rlm = _rlm.Define("Jet_pt_loosejet", "Jet_pt[loosejetcuts]")
               .Define("Jet_phi_loosejet", "Jet_phi[loosejetcuts]")
               .Define("Jet_eta_loosejet", "Jet_eta[loosejetcuts]")
               .Define("Jet_mass_loosejet", "Jet_mass[loosejetcuts]")
               .Define("njetloose", "int(Jet_pt_loosejet.size())")
               .Define("loosejet4vecs", ::gen4vec, {"Jet_pt_loosejet", "Jet_eta_loosejet", "Jet_phi_loosejet", "Jet_mass_loosejet"});

    const std::string vetoPath = "data/JME/" + jetFile + "/jetvetomaps.json.gz";
    std::shared_ptr<correction::CorrectionSet> vetoMapreader;
    try {
        vetoMapreader = loadCorrectionSet(vetoPath);
    } catch (const std::exception &e) {
        std::cerr << "[ERROR] JetVetoMap: Failed to load jetvetomaps JSON '" << vetoPath << "': " << e.what() << std::endl;
        throw;
    }

    correction::Correction::Ref _vetomap;
    try {
        _vetomap = vetoMapreader->at(map);
    } catch (const std::exception &e) {
        std::cerr << "[ERROR] JetVetoMap: Key '" << map << "' not found in '" << vetoPath << "': " << e.what() << std::endl;
        throw;
    }

    auto vetomap = [_vetomap](const std::string &key) {
        return [_vetomap, key](floats &eta, floats &phi) -> floats {
            floats xout;
            xout.reserve(eta.size());
            for (size_t i = 0; i < eta.size(); i++) {
                try {
                    xout.emplace_back(_vetomap->evaluate({key, eta[i], phi[i]}));
                } catch (const std::exception &e) {
                    std::cerr << "[ERROR] JetVetoMap: evaluate failed (key=" << key << " eta=" << eta[i] << " phi=" << phi[i] << "): " << e.what() << std::endl;
                    throw;
                }
            }
            return xout;
        };
    };
    
    // [JME-05] CMS JME POG recommendation (CMS_jmeFAQ.pdf p.6):
    // "Which correctionlib key should be evaluated for a current Run 3 jet veto map? Evaluate jetvetomap."
    // In Run 3, Jet_passJetIdTightLepVeto (with Jet_muEF < 0.8 && Jet_chEmEF < 0.8) inherently excludes
    // prompt leptons, replacing the legacy Run 2 UL manual DeltaR(jet, PF_muon) > 0.2 check.
    // Applied uniformly across all 5 Run 3 eras using "jetvetomap" alone (matching 2023BPix).
    _rlm = _rlm.Define("Jet_isVeto_loose", vetomap("jetvetomap"), {"Jet_eta_loosejet","Jet_phi_loosejet"})
               .Define("events_isVeto","Sum(Jet_isVeto_loose)");
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

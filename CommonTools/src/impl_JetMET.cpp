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

void NanoAODAnalyzerrdframe::setupJetMETCorrection(string jecFile, string jecYear, string jerMap, string jecVersion, bool dataMc){
    auto jercReader = correction::CorrectionSet::from_file("data/JME/"+jecFile+"/jet_jerc.json.gz");
    string datamcflag = "";
    if (dataMc) datamcflag = "DATA";
    else datamcflag = "MC";

    string tmp = jecYear + "_V2_" + datamcflag + "_L1FastJet_AK4PFPuppi";
    std::cout << tmp << std::endl;
    auto _L1FastJet = jercReader->at(jecYear + jecVersion + datamcflag + "_L1FastJet_AK4PFPuppi");
    auto _L2Relative = jercReader->at(jecYear + jecVersion + datamcflag + "_L2Relative_AK4PFPuppi");
    auto _L2L3Residual = jercReader->at(jecYear + jecVersion + datamcflag + "_L2L3Residual_AK4PFPuppi");

    auto applyJes = [this, _L1FastJet, _L2Relative, _L2L3Residual, dataMc](floats jetpts, floats jetetas, floats jetphis, floats jetAreas, floats jetrawf, float rho, unsigned int run, floats toCorr)->floats {
        floats corrfactors;
        corrfactors.reserve(jetpts.size());

        for (size_t i=0; i<jetpts.size(); i++) {
            float rawFact = (1.0-jetrawf[i]);
            float rawjetpt = jetpts[i] * rawFact;
            float L1corr = _L1FastJet->evaluate({jetAreas[i], jetetas[i], rawjetpt, rho});
            float L1pt = rawjetpt * L1corr;
            float L2corr = 1; float L2L3corr = 1;
            L2corr = _L2Relative->evaluate({jetetas[i], jetphis[i], L1pt});
            //if (_isRun24) L2corr = _L2Relative->evaluate({jetetas[i], jetphis[i], L1pt});
            //else  L2corr = _L2Relative->evaluate({jetetas[i], L1pt});
            float L2pt = L1pt * L2corr;
            if (dataMc) L2L3corr = _L2L3Residual->evaluate({float(run), jetetas[i], L2pt});
            float corrfactor = toCorr[i] * rawFact * L1corr * L2corr * L2L3corr; 
            corrfactors.emplace_back(corrfactor);
        }
        return corrfactors;
    };
    _rlm = _rlm.Define("Jet_pt_uncorr", "Jet_pt");
    _rlm = _rlm.Define("Jet_pt_corr", applyJes, {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_area", "Jet_rawFactor", "Rho_fixedGridRhoFastjetAll", "run", "Jet_pt"})
               .Redefine("Jet_mass", applyJes, {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_area", "Jet_rawFactor", "Rho_fixedGridRhoFastjetAll", "run", "Jet_mass"});
    if (!_isData) {
      //_rlm = _rlm.Define("Jet_pt_unc", jesUnc, {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_area", "Jet_rawFactor", "Rho_fixedGridRhoFastjetAll", "Jet_partonFlavour"});
    }
    _rlm = _rlm.Redefine("Jet_pt", "Jet_pt_corr");

    auto jerReader = correction::CorrectionSet::from_file("data/JME/jer_smear.json.gz");
    auto _jerReso = jercReader->at(jerMap + "_JRV1_MC_PtResolution_AK4PFPuppi");
    auto _jerSF = jercReader->at(jerMap + "_JRV1_MC_ScaleFactor_AK4PFPuppi");
    auto _jerSmear = jerReader->at("JERSmear");

    auto applyJer = [this, _jerReso, _jerSF, _jerSmear](floats jetpts, floats jetetas, floats jetphis, shorts genidxs, floats genpts, floats genetas, floats genphis, float rho, ULong64_t event, floats toCorr)->floats {
        floats corrfactors;
        corrfactors.reserve(jetpts.size());

        for (size_t i=0; i<jetpts.size(); i++){
            float reso = _jerReso->evaluate({jetetas[i], jetpts[i], rho});
            float sf = _jerSF->evaluate({jetetas[i], jetpts[i], "nom"});
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
    if (!_isData) {
        _rlm = _rlm.Redefine("Jet_pt_corr", applyJer, {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_genJetIdx", "GenJet_pt", "GenJet_eta", "GenJet_phi", "Rho_fixedGridRhoFastjetAll", "event", "Jet_pt"})
                   .Redefine("Jet_mass", applyJer, {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_genJetIdx", "GenJet_pt", "GenJet_eta", "GenJet_phi", "Rho_fixedGridRhoFastjetAll", "event", "Jet_mass"})
                   .Redefine("Jet_pt", "Jet_pt_corr");
    }
    //_rlm = _rlm.Define("Jet_pt_unc", jesUnc, {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_area", "Jet_rawFactor", "Rho_fixedGridRhoFastjetAll", "Jet_partonFlavour"});
}

void NanoAODAnalyzerrdframe::JetVetoMap(string jetFile, string map) {
    cout << "apply JetVetoMap" << endl;
    

    //TODO: _isRun24 is implementation for NanoAODv15
    if (_isRun24){
        auto jetIDreader = correction::CorrectionSet::from_file("data/JME/2024_Summer24/jetid.json.gz");
        auto _jetid = jetIDreader->at("AK4PUPPI_TightLeptonVeto");
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
        auto getJetId2= [this](floats &etas, uchars &jetIds, floats &muEF, floats &chEmEF, floats &neHEF, floats neEmEF)->floats{
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
    
    auto vetoMapreader = correction::CorrectionSet::from_file("data/JME/"+jetFile+"/jetvetomaps.json.gz");
    
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


/*
 * SkimEvents.cpp
 *
 *  Created on: Dec 9, 2018
 *      Author: suyong
 */

#include "SkimEvents.h"
#include "utility.h"

SkimEvents::SkimEvents(TTree *t, std::string outfilename, std::string year, std::string ch, std::string syst, std::string jsonfname, string globaltag, int nthreads)
:NanoAODAnalyzerrdframe(t, outfilename, year, ch, syst, jsonfname, globaltag, nthreads),_year(year),_ch(ch),_syst(syst)
{
  _isSkim = true;
  _isHTstitching = false;
  if (_outfilename.find("WtoLNu-4Jets") != string::npos) {
      _isHTstitching = true;
  }
}

// Define your cuts here
void SkimEvents::defineCuts()
{
  // Cuts to be applied in order
  // These will be passed to Filter method of RDF
  // check for good json event is defined earlier


  cout << "Skim cut" << endl;
  if (_ch.find("muon") != std::string::npos) {
      addCuts("(HLT_IsoMu24 || HLT_Mu50 || HLT_CascadeMu100 || HLT_HighPtTkMu100) && nmuonpass == 1 && PV_npvsGood > 0 && Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter && Flag_BadPFMuonDzFilter && Flag_hfNoisyHitsFilter && Flag_eeBadScFilter && events_isVeto==0","0"); 
  } else if (_ch.find("electron") != std::string::npos) {
      addCuts("(HLT_Ele30_WPTight_Gsf || HLT_Ele115_CaloIdVT_GsfTrkIdT || HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165 || HLT_Photon200) && nelepass == 1 && PV_npvsGood > 0 && Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter && Flag_BadPFMuonDzFilter && Flag_hfNoisyHitsFilter && Flag_eeBadScFilter && events_isVeto==0", "0");
  }

  //Prescription to fill up WJets HT = 0-100
  //if (_isHTstitching)
  //    addCuts("LHE_HT < 40","00");
}

void SkimEvents::defineMoreVars()
{
        if (_year.find("23") != std::string::npos) {
            addVar({"Flag_filter", "Flag_goodVertices"});
        }
        // define variables that you want to store
        addVartoStore("run");
        addVartoStore("luminosityBlock");
        addVartoStore("event");
        addVartoStore("puWeight.*");
        addVartoStore("unitGenWeight.*");
        addVartoStore("lhereweight.*");
        addVartoStore("muonWeight.*");
        addVartoStore("elecWeight.*");
        addVartoStore("tauWeight.*");
        addVartoStore("LHEPdfWeight");
        addVartoStore("LHEScaleWeight");
        addVartoStore("PSWeight");
        addVartoStore("GenPart_.*");
        addVartoStore("GenJet_.*");
        addVartoStore("nGenVisTau");
        addVartoStore("GenVisTau_.*");
        addVartoStore("gen.*");
        addVartoStore("PV_npvsGood");
        addVartoStore("Pileup_nPU");
        addVartoStore("Pileup_nTrueInt");
        addVartoStore("nJet");
        addVartoStore("Jet_area");
        addVartoStore("Jet_genJetIdx");
        addVartoStore("Jet_hadronFlavour");
        addVartoStore("Jet_partonFlavour");
        addVartoStore("Jet_pt");
        addVartoStore("Jet_pt_uncorr");
        addVartoStore("Jet_pt_unc");
        addVartoStore("Jet_jer");
        addVartoStore("Jet_eta");
        addVartoStore("Jet_phi");
        addVartoStore("Jet_mass");
        addVartoStore("Jet_jetId");
        addVartoStore("Jet_rawFactor");
        addVartoStore("Jet_btagPNetB");
        addVartoStore("Jet_btagUParTAK4B");
        addVartoStore("Jet_muEF");
        addVartoStore("Jet_ch.*");
        addVartoStore("Jet_ne.*");
        addVartoStore("Jet_pass.*");
        addVartoStore("btagWeight_PNet.*");
        addVartoStore("nTau");
        addVartoStore("Tau_charge");
        addVartoStore("Tau_d.*");
        addVartoStore("Tau_eta");
        addVartoStore("Tau_gen.*");
        addVartoStore("Tau_id.*");
        addVartoStore("Tau_jetIdx");
        addVartoStore("Tau_mass");
        addVartoStore("Tau_phi");
        addVartoStore("Tau_pt");
        addVartoStore("Tau_pt_unc");
        addVartoStore("Tau_pt_uncor");
        addVartoStore("Tau_puCorr");
        addVartoStore("Tau_rawDeep.*");
        addVartoStore("Tau_rawIso.*");
        addVartoStore("MET_p.*");
        addVartoStore("MET_sumEt");
        addVartoStore("MET_MetUnclust.*");
        addVartoStore("PuppiMET_p.*");
        addVartoStore("PuppiMET_sumEt");
        addVartoStore("PuppiMET_MetUnclust.*");
        addVartoStore("RawPuppiMET.*");
        addVartoStore("LHE_HT");
        addVartoStore("nElectron");
        addVartoStore("Electron_charge");
        addVartoStore("Electron_cutBased");
        addVartoStore("Electron_mvaIso.*");
        addVartoStore("Electron_deltaEtaSC");
        addVartoStore("Electron_dxy.*");
        addVartoStore("Electron_dz.*");
        addVartoStore("Electron_eta");
        addVartoStore("Electron_mass");
        addVartoStore("Electron_pf.*");
        addVartoStore("Electron_phi");
        addVartoStore("Electron_pt");
        addVartoStore("nMuon");
        addVartoStore("Muon_charge");
        addVartoStore("Muon_eta");
        addVartoStore("Muon_mass");
        addVartoStore("Muon_pfRelIso04_all");
        addVartoStore("Muon_phi");
        addVartoStore("Muon_pt.*");
        addVartoStore("Muon_tightId");
        addVartoStore("Muon_looseId");
        addVartoStore("nmuonpass");
        addVartoStore("nelepass");
        addVartoStore("nvetomuons");
        addVartoStore("nvetoelepass");
        addVartoStore("lep4vecs");
        addVartoStore("Rho_fixedGridRhoFastjetAll");
        addVartoStore("L1PreFiringWeight_.*");
        addVartoStore("LHEPart_pt");
        addVartoStore("LHEPart_pdgId");
}

void SkimEvents::bookHists()
{
	// _hist1dinfovector contains the information of histogram definitions (as TH1DModel)
	// the variable to be used for filling
	// and the minimum cutstep for which the histogram should be filled
	//
	// The braces are used to initalize the struct
	// TH1D
    add1DHist( {"hcounter", "Number of events", 2, -0.5, 1.5}, "one", "unitGenWeight", "", "", "");
}

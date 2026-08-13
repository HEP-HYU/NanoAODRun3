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
  cout << "<< Start Skim NanoAOD >>" << endl;
  _isHTstitching = false;
  if (_outfilename.find("WtoENu-4Jets")!= string::npos || _outfilename.find("WtoMuNu-4Jets")!= string::npos || _outfilename.find("WtoTauNu-4Jets")!= string::npos) {
      cout << "WJet stitching" << endl;
      _isHTstitching = true;
  }
  if (_ch.find("muon") != std::string::npos){
      cout << "muon channel" << endl;
      _isMuonCh = true;
  } else{
      cout << "electron channel" << endl;
      _isMuonCh = false;
  }

  if (_outfilename.find("Tau-LFV-") != std::string::npos){
      cout << "Signal" << endl;
      _isSignal = true;
  } else{
      _isSignal = false;
  }
}

void SkimEvents::defineObjectSelection(std::vector<std::string> jes_var){
    std::string pileFile = "";
    std::string pileMap = "";
    std::string jetFile = "";
    std::string jetMap = "";
    std::string tauYear = "";
    std::string jecFile = "";
    std::string jecYear = "";      // MC JEC key prefix (no run-era label)
    std::string jecYearData = "";  // DATA JEC key prefix (includes run-era label)
    std::string jerMap = "";
    std::string elecFile = "";
    std::string elecYear = "";
    std::string jecVersion = "_V2_";
    std::string jerVersion  = "JRV1";  // JER version tag: JRV1 for all eras, JRV2 for 2024
    if (_isRun22) {
        pileFile = "2022_Summer22";
        pileMap = "Collisions2022_355100_357900_eraBCD_GoldenJson";
        jetFile = "2022_Summer22";
        jetMap = "Summer22_23Sep2023_RunCD_V1";
        elecFile = "2022_Summer22";
        elecYear = "2022Re-recoBCD";
        tauYear = "2022_preEE";
        jecFile = "2022_Summer22";
        jecYear = "Summer22_22Sep2023";         // MC
        jecYearData = "Summer22_22Sep2023_RunCD"; // DATA
        jerMap = "Summer22_22Sep2023";
    } else if (_isRun22EE) {
        pileFile = "2022_Summer22EE";
        pileMap = "Collisions2022_359022_362760_eraEFG_GoldenJson";
        jetFile = "2022_Summer22EE";
        jetMap = "Summer22EE_23Sep2023_RunEFG_V1";
        elecFile = "2022_Summer22EE";
        elecYear = "2022Re-recoE+PromptFG";
        tauYear = "2022_postEE";
        jecFile = "2022_Summer22EE";
        jecYear = "Summer22EE_22Sep2023";  // MC
        // DATA: Run E/F/G have separate keys. Derive from globaltag set by skimonefile.py.
        // e.g. _globaltag="Summer22EE_22Sep2023_RunE" → jecYearData="Summer22EE_22Sep2023_RunE"
        // Fall back to RunF (most common) if not determinable.
        if      (_globaltag.find("RunE") != std::string::npos) jecYearData = "Summer22EE_22Sep2023_RunE";
        else if (_globaltag.find("RunG") != std::string::npos) jecYearData = "Summer22EE_22Sep2023_RunG";
        else                                                    jecYearData = "Summer22EE_22Sep2023_RunF";
        jerMap = "Summer22EE_22Sep2023";
    } else if (_isRun23) {
        pileFile = "2023_Summer23";
        pileMap = "Collisions2023_366403_369802_eraBC_GoldenJson";
        jetFile = "2023_Summer23";
        jetMap = "Summer23Prompt23_RunC_V1";
        elecFile = "2023_Summer23";
        elecYear = "2023PromptC";
        tauYear = "2023_preBPix";
        jecFile = "2023_Summer23";
        jecYear = "Summer23Prompt23";           // MC
        jecYearData = "Summer23Prompt23_RunCv4"; // DATA
        jerMap = "Summer23Prompt23_RunCv1234";
        jecVersion = "_V1_";  // 2023 jet_jerc.json only has _V1_ keys
    } else if (_isRun23BPix) {
        pileFile = "2023_Summer23BPix";
        pileMap = "Collisions2023_369803_370790_eraD_GoldenJson";
        jetFile = "2023_Summer23BPix";
        jetMap = "Summer23BPixPrompt23_RunD_V1";
        elecFile = "2023_Summer23BPix";
        elecYear = "2023PromptD";
        tauYear = "2023_postBPix";
        jecFile = "2023_Summer23BPix";
        jecYear = "Summer23BPixPrompt23";           // MC
        jecYearData = "Summer23BPixPrompt23_RunD";  // DATA
        jerMap = "Summer23BPixPrompt23_RunD";
        jecVersion = "_V1_";  // 2023BPix jet_jerc.json only has _V1_ keys
    } else if (_isRun24) {
        //TODO
        pileFile = "2024_Summer24";
        pileMap = "Collisions24_BCDEFGHI_goldenJSON";
        jetFile = "2024_Summer24";
        jetMap = "Summer24Prompt24_RunBCDEFGHI_V1";
        elecFile = "2024_Summer24";
        elecYear = "2024Prompt";
        tauYear = "2024";
        jecFile = "2024_Summer24";
        jecYear = "Summer24Prompt24";   // MC
        // DATA: Summer24Prompt24 keys have NO run-era label — same prefix as MC
        // e.g. Summer24Prompt24_V5_DATA_L1FastJet_AK4PFPuppi
        jecYearData = "Summer24Prompt24";
        jerMap = "Summer24Prompt24";
        jecVersion = "_V5_";  // jet_jerc.json.gz updated to V5 (fixed L2L3Residual bug)
        jerVersion  = "JRV2"; // jet_jerc.json.gz updated to JRV2
    }

    // MET XY correction epoch: maps era to the epoch string used inside the JSON.
    // 2022 -> "2022", 2022EE -> "2022EE", 2023 -> "2023", 2023BPix -> "2023BPix", 2024 -> "" (no file).
    std::string metEpoch = "";
    if      (_isRun22)     metEpoch = "2022";
    else if (_isRun22EE)   metEpoch = "2022EE";
    else if (_isRun23)     metEpoch = "2023";
    else if (_isRun23BPix) metEpoch = "2023BPix";
    // 2024: no met_xyCorrections file available yet, metEpoch left empty.

    // Electron cut: CMS EGM Run3 recommendation
    //   - MVA-based WP90 isolation ID (mvaIso_WP90)
    //   - Exclude ECAL gap (1.4442 < |eta_SC| < 1.566) where reconstruction degrades
    //   - Electron_eta here is the supercluster eta (Electron_eta in NanoAOD = SC eta)
    std::string muoncut  = "Muon_pt>50.0 && abs(Muon_eta)<2.4 && Muon_tightId && Muon_pfRelIso04_all<0.15";
    std::string vetomuon = "!muoncuts && Muon_pt>15.0 && abs(Muon_eta)<2.4 && Muon_looseId && Muon_pfRelIso04_all<0.25";
    // Signal electron: WP90iso + pt>50 + |eta|<2.5 + ECAL gap exclusion
    std::string eleccut  = "Electron_pt>50 && abs(Electron_eta)<2.5 && !(abs(Electron_eta)>1.4442 && abs(Electron_eta)<1.566) && Electron_mvaIso_WP90";
    // Veto electron: loose cutBased + |eta|<2.5 + ECAL gap exclusion
    std::string vetoelec = "!elecuts && Electron_pt>15.0 && abs(Electron_eta)<2.5 && !(abs(Electron_eta)>1.4442 && abs(Electron_eta)<1.566) && Electron_cutBased == 1";
    std::string skimjet = "Jet_pt>30.0 && abs(Jet_eta)<2.5 && (Jet_passJetIdTightLepVeto==1.0) && Jet_muEF<0.8 && Jet_chEmEF<0.8";

    std::string muonid = "NUM_TightID_DEN_TrackerMuons";
    std::string muoniso = "NUM_TightPFIso_DEN_TightID";
    std::string muonhlt = "NUM_IsoMu24_or_Mu50_or_CascadeMu100_or_HighPtTkMu100_DEN_CutBasedIdTight_and_PFIsoTight";

    if (_isData){
        noiseFilter();
    }
    setupJetMETCorrection(jecFile, jecYear, jerMap, jecVersion, _isData, jecYearData, jerVersion);
    // Apply MET XY-shift correction (PuppiMET + PFMET) after JEC but before object selection.
    // Skipped automatically for 2024 (no correction file available).
    if (!metEpoch.empty()) applyMETXYCorrection(jecFile, metEpoch);
    JetVetoMap(jetFile, jetMap);
    applyWeights(pileFile, pileMap);
    if (_ch.find("muon") != std::string::npos){
        if (!_isData){
            calculateTauES(tauYear, "VTight", "Tight", "VVLoose");
            calculateMuonSF(muonid, muoniso, muonhlt);
        }
        selectMuons(muoncut, vetomuon);
    } else {
        if (!_isData){
            calculateTauES(tauYear, "VTight", "VLoose", "Tight");
            calculateElectronSF(elecFile, elecYear);
        }
        selectElectrons(eleccut, vetoelec);
    }
    skimJets(skimjet);
    if (_isSignal) matchGenReco();
}

// Define your cuts here
void SkimEvents::defineCuts()
{
  // Cuts to be applied in order
  // These will be passed to Filter method of RDF
  // check for good json event is defined earlier

  cout << "Skim cut" << endl;

  if (_ch.find("muon") != std::string::npos) {
      addCuts("(HLT_IsoMu24 || HLT_Mu50 || HLT_CascadeMu100 || HLT_HighPtTkMu100) && nmuonpass == 1 && PV_npvsGood > 0 && Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter && Flag_BadPFMuonDzFilter && Flag_hfNoisyHitsFilter && Flag_eeBadScFilter && Flag_ecalBadCalibFilter && events_isVeto==0","0");
      if (_isSignal){
          addCuts("GenPart_d_SMW1_idx >= 0 && GenPart_d_SMW2_idx >= 0", "00");
      }
  } else if (_ch.find("electron") != std::string::npos) {
      addCuts("(HLT_Ele30_WPTight_Gsf || HLT_Ele115_CaloIdVT_GsfTrkIdT || HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165 || HLT_Photon200) && nelepass == 1 && PV_npvsGood > 0 && Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter && Flag_BadPFMuonDzFilter && Flag_hfNoisyHitsFilter && Flag_eeBadScFilter && Flag_ecalBadCalibFilter && events_isVeto==0", "0");
      if (_isSignal){
          addCuts("GenPart_d_SMW1_idx >= 0 && GenPart_d_SMW2_idx >= 0", "00");
      }
  }

  //Prescription to fill up WJets 0J
  if (_isHTstitching) {
      addCuts("LHE_Njets == 0","00");
  }
}

void SkimEvents::defineMoreVars()
{
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
        addVartoStore("Jet_passJetIdTightLepVeto");
        addVartoStore("Jet_rawFactor");
        addVartoStore("Jet_btagPNet.*");
        addVartoStore("Jet_btagUParTAK4.*");
        addVartoStore("Jet_muEF");
        addVartoStore("Jet_ch.*");
        addVartoStore("Jet_ne.*");
        addVartoStore("Jet_pass.*");
        addVartoStore("btagWeight.*");
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
        addVartoStore("Tau_pt_uncor");
        addVartoStore("Tau_ES_nom");
        addVartoStore("Tau_ES_up");
        addVartoStore("Tau_ES_down");

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

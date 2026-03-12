/*
 * ttHHAnalyzer.cpp
 *
 *  Created on: April 9, 2020
 *      Author: Tae Jeong Kim
 *      Refactored: Jiwon Park
 */

#include "ttHHAnalyzer.h"
#include "utility.h"

ttHHAnalyzer::ttHHAnalyzer(TTree *t, std::string outfilename, std::string year, std::string ch, std::string syst, std::string jsonfname, bool applytauFF, string globaltag, int nthreads)
:NanoAODAnalyzerrdframe(t, outfilename, year, ch, syst, jsonfname, globaltag, nthreads), _outfilename(outfilename), _syst(syst), _year(year), _ch(ch), _applytauFF(applytauFF)
{

    cout << "<< Start Process NanoAOD >>" << endl;

    if(syst.find("jes") != std::string::npos or syst.find("jer") != std::string::npos or syst.find("metUnclust") != std::string::npos or
            syst.find("tes") != std::string::npos or syst.find("hdamp") != std::string::npos or syst.find("tune") != std::string::npos or
            syst.find("muonhighscale") != std::string::npos) {
        ext_syst = true;
    }
    _syst = syst;
    cout << "syst: " << _syst << endl;

    if (_outfilename.find("ST_LFV") != std::string::npos || _outfilename.find("TT_LFV") != std::string::npos) {
        _isSignal = true;
        cout << "Input file is LFV signal" <<endl;
    } else {
        _isSignal = false;
    }
}

void ttHHAnalyzer::defineObjectSelection(std::vector<std::string> jes_var){
    std::string cut  = "onlyveto";
    std::string vetomuon = "Muon_pt>15.0 && abs(Muon_eta)<2.4 && Muon_looseId && Muon_pfRelIso04_all<0.25";
    std::string vetoelec = "Electron_pt>15.0 && abs(Electron_eta)<2.5 && Electron_cutBased == 1";
    std::string muonid = "NUM_TightID_DEN_TrackerMuons";
    std::string muoniso = "NUM_TightPFIso_DEN_TightID";
    std::string muonhlt = "NUM_IsoMu24_or_Mu50_or_CascadeMu100_or_HighPtTkMu100_DEN_CutBasedIdTight_and_PFIsoTight";
    std::string jetcut = "Jet_pt>40.0 && abs(Jet_eta)<2.4 && Jet_passJetIdTightLepVeto && Jet_muEF < 0.8 && Jet_chEmEF < 0.8";
    std::string taucut = "Tau_pt>30.0 && abs(Tau_eta)<2.5  && Tau_idDecayModeNewDMs && Tau_decayMode != 5 && Tau_decayMode != 6";

    calculateMuonSF(muonid, muoniso, muonhlt);
    selectElectrons(cut, vetoelec);
    selectJets(jes_var, jes_var_flav, jetcut);
    if (!_isData){
        topPtReweight();
    }
}

// Define your cuts here
void ttHHAnalyzer::defineCuts() {
    cout << "Apply event selections" << endl;
    addCuts("nmuonpass >= 2 && nvetoelepass == 0 && nvetomuons == 0 && PV_npvsGood > 0", "0");
    addCuts("ncleanjetspass >= 4", "00");
    addCuts("ncleanbjetspass >= 3", "000");
}

void ttHHAnalyzer::defineMoreVars() {
    defineVar("lepvec", ::select_leadingvec, {"lep4vecs"});
    defineVar("lepMET_mt", ::calculate_MT, {"lep4vecs","PuppiMET_pt","PuppiMET_phi"});
    
    addVar({"Muon1_pt", "Muon_pt[0]", ""});
    addVar({"Muon1_eta", "Muon_eta[0]", ""});
    addVar({"Muon1_mass", "Muon_mass[0]", ""});
    addVar({"Muon1_charge", "Muon_charge[0]", ""});
    
    addVar({"Jet1_pt", "(Jet_pt.size()>0) ? Jet_pt[0] : -1", ""});
    addVar({"Jet1_eta", "Jet_eta[0]", ""});
    addVar({"Jet1_mass", "Jet_mass[0]", ""});
    addVar({"Jet1_btagPNetB","Jet_btagPNetB[0]",""});

    addVar({"Jet2_pt", "Jet_pt[1]", ""});
    addVar({"Jet2_eta", "Jet_eta[1]", ""});
    addVar({"Jet2_mass", "Jet_mass[1]", ""});
    addVar({"Jet2_btagPNetB","Jet_btagPNetB[1]",""});

    addVar({"Jet3_pt", "Jet_pt[2]", ""});
    addVar({"Jet3_eta", "Jet_eta[2]", ""});
    addVar({"Jet3_mass", "Jet_mass[2]", ""});
    addVar({"Jet3_btagPNetB","Jet_btagPNetB[2]",""});

    addVar({"bJet1_pt", "bJet_pt[0]", ""});
    addVar({"bJet1_eta", "bJet_eta[0]", ""});
    addVar({"bJet1_mass", "bJet_mass[0]", ""});


    if (_syst == "data" or _isSignal) {
        addVar({"eventWeight", "1.0"});
        addVar({"eventWeight_notau", "1.0"});
        addVar({"eventWeight_notoppt", "1.0"});
        addVar({"eventWeight_nobtag", "1.0"});
        addVar({"eventWeight_notau_nobtag", "1.0"});
        addVar({"eventWeight_pu", "1"});
        addVar({"eventWeight_mu", "1"});
        addVar({"eventWeight_pumu", "1"});
        addVar({"eventWeight_elec", "1"});
        addVar({"eventWeight_puelec", "1"});
        addVar({"eventWeight_all", "1"});
    } else {
        if (_syst == "" or _syst == "nosyst") { //TODO
            addVar({"eventWeight_genpu", "1.0"});
            addVar({"eventWeight_genpumu", "1.0"});
            addVar({"eventWeight_notau_nobtag", "1.0"}); //didn't want to duplicate entry...
            addVar({"eventWeight_nobtag", "1.0"});
            addVar({"eventWeight_nopu", "1.0"});
            addVar({"eventWeight_noprefire", "1.0"});
            addVar({"eventWeight_notoppt", "1.0"});
        } else {
            addVar({"eventWeight_genpu", "unitGenWeightFF * TopPtWeight[0] * puWeight[0] * L1PreFiringWeight_Nom"});
            addVar({"eventWeight_mu", "muonWeightId[0] * muonWeightIso[0] * muonWeightTrg[0]"});
            addVar({"eventWeight_genpumu", "eventWeight_genpu * eventWeight_mu"});
            addVar({"eventWeight_notau_nobtag", "eventWeight_genpumu"}); //didn't want to duplicate entry...
        }

        if (_syst == "" or _syst == "nosyst" or ext_syst) {
            // for external syst, we only need nominal weight
            if (_syst == "" or _syst == "nosyst") {
                addVar({"eventWeight", "unitGenWeight"});
                addVar({"eventWeight_top", "TopPtWeight[0]"});
                addVar({"eventWeight_notau", "eventWeight"});
                addVar({"eventWeight_pu", "eventWeight * puWeight[0]"});
                addVar({"eventWeight_mu", "muonWeightId[0] * muonWeightIso[0] * muonWeightTrg[0]"});
                addVar({"eventWeight_pumu", "eventWeight_pu * eventWeight_mu"});
                addVar({"eventWeight_all", "eventWeight_pu * eventWeight_mu"});
            } else {
                addVar({"eventWeight", "eventWeight_nobtag * btagWeight_PNetB[0]"});
                addVar({"eventWeight_notau", "eventWeight_genpumu * btagWeight_PNetB[0]"});
            }
        } else {
            addVar({"eventWeight", "eventWeight_nobtag * btagWeight_PNetB[0]"});
            addVar({"eventWeight_notau", "eventWeight_genpu * eventWeight_mu * btagWeight_PNetB[0]"});
        }
    }
    
    // define variables that you want to store
    addVartoStore("run");
    addVartoStore("event");
    addVartoStore("luminosityBlock");
    addVartoStore("eventWeight.*");
    addVartoStore("nmuonpass");
    addVartoStore("nelepass");
    addVartoStore("ncleanjetspass");
    addVartoStore("ncleanbjetspass");
    addVartoStore("Jet_pt");
    addVartoStore("Jet1_.*");
    addVartoStore("Jet2_.*");
    addVartoStore("Jet3_.*");
    addVartoStore("bJet1_.*");
    addVartoStore("Muon1.*");
    addVartoStore("PuppiMET_pt");
    addVartoStore("PuppiMET_phi");
    addVartoStore("chi2.*");
    addVartoStore("btagWeight_PNetB");
    addVartoStore("btagWeight_PNetB_jes");
    addVartoStore("GenPart_top_pt");
    addVartoStore("TopPtWeight");
    addVartoStore("LHEPart_pt");
    addVartoStore("LHEPart_pdgId");
}

void ttHHAnalyzer::bookHists() {
    std::string lep = "e", title_tmp = "";
    if (_isMuonCh) lep = "#mu"; 
    std::vector<std::string> init_weight = {""};
    std::vector<std::string> sf_weight = {"", "_pu", "_all"};
    sf_weight = {"", "_pu", "_mu", "_pumu", "_all"};

    std::vector<std::string> sf_weight_tau = {};

    std::vector<std::string> sf_weight_FF = {};

    std::vector<std::string> theory_weight = {};

   

    //Step definition. To replace for tauFF
    //This implies histos w/o FF must be drawn, for bSF and other weights renormalization
    std::string minstep_S1 = "0";
    std::string minstep_S2 = "00";
    std::string minstep_S3 = "000";

    // TODO refine this. Too heavy? we will see.
    std::vector<std::string> syst_weight;
    if (_syst == "" or _syst == "nosyst" or _syst == "data" or ext_syst) {
        syst_weight = init_weight;
        if (_syst != "data") {
            //We anyway need this for bSF rescaling
            add1DHist({"h_nevents", ";Number of events w/o b SF;Events", 2, -0.5, 1.5}, "one", "eventWeight", "_nobtag", minstep_S1, "");
            add1DHist({"h_nevents_notausf", ";Number of events w/o b and tau SF;Events", 2, -0.5, 1.5}, "one", "eventWeight_notau", "_nobtag", minstep_S1, "00");
            add1DHist({"h_jet_ht", ";Jet HT w/o b SF (GeV);Events", 48, 40, 1000}, "Jet_HT", "eventWeight", "_nobtag", minstep_S1, "");
        }
    }
    else {
        syst_weight = sf_weight;
        if (_syst == "theory") syst_weight.insert(syst_weight.end(), theory_weight.begin(), theory_weight.end());
    }

    // S1 w/o tau SF
    maxstep = "00"; //Must be +1 step than its cut
    for (std::string weightstr : syst_weight) {
        if (weightstr.find("notau") != std::string::npos) continue;
        if (weightstr.find("notoppt") != std::string::npos) continue;

        add1DHist({"h_nevents_notausf", ";Number of events;Events", 2, -0.5, 1.5}, "one", "eventWeight_notau", weightstr, minstep_S1, maxstep);
        add1DHist({"h_nvtx_notausf", ";Number of primary vertex;Events", 70, 0.0, 70.0}, "PV_npvsGood", "eventWeight_notau", weightstr, minstep_S1, maxstep);

        add1DHist({"h_puppimet_pt_notausf", ";PuppiMET (GeV);Events", 20, 0, 400}, "PuppiMET_pt", "eventWeight_notau", weightstr, minstep_S1, maxstep);
        add1DHist({"h_puppisum_et_notausf", ";Sum ET;Events", 50, 0.0, 5000.0}, "PuppiMET_sumEt", "eventWeight_notau", weightstr, minstep_S1, maxstep);
        add1DHist({"h_puppimet_phi_notausf", ";PuppiMET #phi;Events", 20, -3.2, 3.2}, "PuppiMET_phi", "eventWeight_notau", weightstr, minstep_S1, maxstep);

        add1DHist({"h_ncleanjetspass_notausf", ";Jet multiplicity;Events", 10, 0.0, 10.0}, "ncleanjetspass", "eventWeight_notau", weightstr, minstep_S1, maxstep);
        add1DHist({"h_ncleanbjetspass_notausf", ";b-tagged jet multiplicity;Events", 5, 0.0, 5.0}, "ncleanbjetspass", "eventWeight_notau", weightstr, minstep_S1, maxstep);

        add1DHist({"h_nmuonpass_notausf", ";Muon multiplicity;Events", 5, 0.0, 5.0}, "nmuonpass", "eventWeight_notau", weightstr, minstep_S1, maxstep);
        add1DHist({"h_muon1_pt_notausf", ";Muon p_{T} (GeV);Events", 10, 0, 200}, "Muon1_pt", "eventWeight_notau", weightstr, minstep_S1, maxstep);
        add1DHist({"h_muon1_eta_notausf", ";Muon #eta;Events", 20, -2.4, 2.4}, "Muon1_eta", "eventWeight_notau", weightstr, minstep_S1, maxstep);
        title_tmp = std::string(";m_{T}(")+lep+", MET) (GeV);Events";
        add1DHist({"h_lepMET_mt_notausf", title_tmp.c_str(), 20, 0, 400}, "lepMET_mt", "eventWeight_notau", weightstr, minstep_S1, maxstep);

        add1DHist({"h_jet1_pt_notausf", ";Leading jet p_{T} (GeV);Events", 20, 0, 400}, "Jet1_pt", "eventWeight_notau", weightstr, minstep_S1, maxstep);
        add1DHist({"h_jet1_eta_notausf", ";Leading jet #eta;Events", 20, -2.4, 2.4}, "Jet1_eta", "eventWeight_notau", weightstr, minstep_S1, maxstep);
        add1DHist({"h_jet1_mass_notausf", ";Leading jet mass (GeV);Events", 20, 0, 100}, "Jet1_mass", "eventWeight_notau", weightstr, minstep_S1, maxstep);
        add1DHist({"h_jet1_btag_notausf",";PNet score of leading jet;Events", 20, 0, 1.0}, "Jet1_btagPNetB", "eventWeight_notau", weightstr, minstep_S1, maxstep);
    }

    //for all the other nominal histograms with tauSF
    if (_syst == "all" or _syst == "theory") {
        syst_weight.insert(syst_weight.end(), sf_weight_tau.begin(), sf_weight_tau.end());
        if (_applytauFF) syst_weight.insert(syst_weight.end(), sf_weight_FF.begin(), sf_weight_FF.end());
    }

    cout << "Variations to take care :";
    for (auto i : syst_weight) cout << i << " ";
    cout << endl;

    maxstep = "";

    //TODO: This part is just for test
    syst_weight.insert(syst_weight.end(), sf_weight.begin(), sf_weight.end());
    for (std::string weightstr : syst_weight) {
        if (weightstr.find("notoppt") != std::string::npos) continue;

        add1DHist({"h_nevents", ";Number of events;Events", 2, -0.5, 1.5}, "one", "eventWeight", weightstr, minstep_S1, maxstep);
        add1DHist({"h_nvtx", ";Number of primary vertex;Events", 70, 0.0, 70.0}, "PV_npvsGood", "eventWeight", weightstr, minstep_S1, maxstep);
        add1DHist({"h_puppimet_pt", ";PuppiMET (GeV);Events", 20, 0, 400}, "PuppiMET_pt", "eventWeight", weightstr, minstep_S1, maxstep);
        add1DHist({"h_puppisum_et", ";Sum ET;Events", 50, 0.0, 5000.0}, "PuppiMET_sumEt", "eventWeight", weightstr, minstep_S1, maxstep);
        add1DHist({"h_puppimet_phi", ";PuppiMET #phi;Events", 20, -3.2, 3.2}, "PuppiMET_phi", "eventWeight", weightstr, minstep_S1, maxstep);

        add1DHist({"h_ncleanjetspass", ";Jet multiplicity;Events", 10, 0.0, 10.0}, "ncleanjetspass", "eventWeight", weightstr, minstep_S1, maxstep);
        add1DHist({"h_ncleanbjetspass", ";b-tagged jet multiplicity;Events", 5, 0.0, 5.0}, "ncleanbjetspass", "eventWeight", weightstr, minstep_S1, maxstep);

        add1DHist({"h_nmuonpass", ";Muon multiplicity;Events", 5, 0.0, 5.0}, "nmuonpass", "eventWeight", weightstr, minstep_S1, maxstep);
        add1DHist({"h_muon1_pt", ";Muon p_{T} (GeV);Events", 10, 0, 200}, "Muon1_pt", "eventWeight", weightstr, minstep_S1, maxstep);
        add1DHist({"h_muon1_eta", ";Muon #eta;Events", 20, -2.4, 2.4}, "Muon1_eta", "eventWeight", weightstr, minstep_S1, maxstep);

        title_tmp =  std::string(";m_{T}(")+lep+", MET) (GeV);Events";
        add1DHist({"h_lepMET_mt", title_tmp.c_str(), 20, 0, 400}, "lepMET_mt", "eventWeight", weightstr, minstep_S1, maxstep);

        add1DHist({"h_jet1_pt", ";Leading jet p_{T} (GeV);Events", 20, 0, 400}, "Jet1_pt", "eventWeight", weightstr, minstep_S1, maxstep);
        add1DHist({"h_jet1_eta", ";Leading jet #eta;Events", 20, -2.4, 2.4}, "Jet1_eta", "eventWeight", weightstr, minstep_S1, maxstep);
        add1DHist({"h_jet1_mass", ";Leading jet mass (GeV);Events", 20, 0, 100}, "Jet1_mass", "eventWeight", weightstr, minstep_S1, maxstep);
        add1DHist({"h_jet1_btag",";PNet score of leading jet;Events", 20, 0, 1.0}, "Jet1_btagPNetB", "eventWeight", weightstr, minstep_S1, maxstep);

        add1DHist({"h_jet2_pt", ";Subleading jet p_{T} (GeV);Events", 20, 0, 400}, "Jet2_pt", "eventWeight", weightstr, minstep_S2, maxstep);
        add1DHist({"h_jet2_eta", ";Subleading jet #eta;Events", 20, -2.4, 2.4}, "Jet2_eta", "eventWeight", weightstr, minstep_S2, maxstep);
        add1DHist({"h_jet2_mass", ";Subleading jet mass (GeV);Events", 20, 0, 100}, "Jet2_mass", "eventWeight", weightstr, minstep_S2, maxstep);
        add1DHist({"h_jet2_btag",";PNet score of subleading jet;Events", 20, 0, 1.0}, "Jet2_btagPNetB", "eventWeight", weightstr, minstep_S2, maxstep);

        add1DHist({"h_jet3_pt", ";Subsubleading jet p_{T} (GeV);Events", 20, 0, 400}, "Jet3_pt", "eventWeight", weightstr, minstep_S2, maxstep);
        add1DHist({"h_jet3_eta", ";Subsubleading jet #eta;Events", 20, -2.4, 2.4}, "Jet3_eta", "eventWeight", weightstr, minstep_S2, maxstep);
        add1DHist({"h_jet3_mass", ";Subsubleading jet mass (GeV);Events", 20, 0, 100}, "Jet3_mass", "eventWeight", weightstr, minstep_S2, maxstep);
        add1DHist({"h_jet3_btag",";PNet score of subsubleading jet;Events", 20, 0, 1.0}, "Jet3_btagPNetB", "eventWeight", weightstr, minstep_S2, maxstep);

        add1DHist({"h_bjet1_pt", ";b-tagged jet p_{T} (GeV);Events", 20, 0, 400}, "bJet1_pt", "eventWeight", weightstr, minstep_S3, maxstep);
        add1DHist({"h_bjet1_eta", ";b-tagged jet #eta;Events", 20, -2.4, 2.4}, "bJet1_eta", "eventWeight", weightstr, minstep_S3, maxstep);
        add1DHist({"h_bjet1_mass", ";b-tagged jet mass (GeV);Events", 20, 0, 100}, "bJet1_mass", "eventWeight", weightstr, minstep_S3, maxstep);

        add1DHist({"h_jet_ht", ";Jet HT (GeV);Events", 48, 40, 1000}, "Jet_HT", "eventWeight", weightstr, minstep_S1, maxstep);

    }
}

double ttHHAnalyzer::tauFF(std::string year_, std::string unc_, int direction_, floats &tau_pt_, floats &tau_gen_pt_, ints &tau_dm_) {

    double val = 1.0;

    // For geniune tau, unc and SF are always 1.0
    if (tau_gen_pt_.size() > 0) return 1.0;

    //Assume exactly one tau for FF (S5)
    //if (tau_pt_.size() != 1) return 9999999.;

    if (abs(direction_) != 1 and direction_ != 0) return val;
    std::map<std::string, std::map<std::string, double>> map_ff;

    // cleaned jet, dm binned
    if (tau_dm_[0] == 0) {
        map_ff["2016pre"]["nom"]  = 0.3029;
        map_ff["2016pre"]["stat"] = 0.1235;
        map_ff["2016pre"]["syst"] = 0.3515;
        map_ff["2016post"]["nom"]  = 0.6584;
        map_ff["2016post"]["stat"] = 0.1212;
        map_ff["2016post"]["syst"] = 0.3563;
        map_ff["2017"]["nom"]  = 0.6122;
        map_ff["2017"]["stat"] = 0.06591;
        map_ff["2017"]["syst"] = 0.2855;
        map_ff["2018"]["nom"]  = 0.5042;
        map_ff["2018"]["stat"] = 0.05853;
        map_ff["2018"]["syst"] = 0.2186;
    } else if (tau_dm_[0] == 1) {
        map_ff["2016pre"]["nom"]  = 0.5605;
        map_ff["2016pre"]["stat"] = 0.1067;
        map_ff["2016pre"]["syst"] = 0.2431;
        map_ff["2016post"]["nom"]  = 0.4114;
        map_ff["2016post"]["stat"] = 0.117;
        map_ff["2016post"]["syst"] = 0.006893;
        map_ff["2017"]["nom"]  = 0.5614;
        map_ff["2017"]["stat"] = 0.05839;
        map_ff["2017"]["syst"] = 0.2171;
        map_ff["2018"]["nom"]  = 0.5407;
        map_ff["2018"]["stat"] = 0.05101;
        map_ff["2018"]["syst"] = 0.1356;
    } else if (tau_dm_[0] == 10) {
        map_ff["2016pre"]["nom"]  = 0.5634;
        map_ff["2016pre"]["stat"] = 0.147;
        map_ff["2016pre"]["syst"] = 0.3097;
        map_ff["2016post"]["nom"]  = 0.4488;
        map_ff["2016post"]["stat"] = 0.159;
        map_ff["2016post"]["syst"] = 0.2005;
        map_ff["2017"]["nom"]  = 0.7205;
        map_ff["2017"]["stat"] = 0.08742;
        map_ff["2017"]["syst"] = 0.2134;
        map_ff["2018"]["nom"]  = 0.9003;
        map_ff["2018"]["stat"] = 0.06959;
        map_ff["2018"]["syst"] = 0.224;
    } else if (tau_dm_[0] == 11) {
        map_ff["2016pre"]["nom"]  = 0.3913;
        map_ff["2016pre"]["stat"] = 0.2608;
        map_ff["2016pre"]["syst"] = 0.1003;
        map_ff["2016post"]["nom"]  = 1.384;
        map_ff["2016post"]["stat"] = 0.2204;
        map_ff["2016post"]["syst"] = 0.09176;
        map_ff["2017"]["nom"]  = 0.208;
        map_ff["2017"]["stat"] = 0.168;
        map_ff["2017"]["syst"] = 0.2655;
        map_ff["2018"]["nom"]  = 0.3594;
        map_ff["2018"]["stat"] = 0.1152;
        map_ff["2018"]["syst"] = 0.1831;
    } else {
        std::cout << "fatal: wrong decay mode detected" << std::endl;
    }

    if      (unc_ == "nom") val = map_ff[year_][unc_];
    else if (direction_ == 1) val  = 1.0 + map_ff[year_][unc_];
    else if (direction_ == -1) val = 1.0 - map_ff[year_][unc_];

    return val;
}

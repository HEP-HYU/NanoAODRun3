/*
 * TauFakeFactorAnalyzer.cpp
 *
 *  Created on: June 14, 2023
 *      Author: Jiwon Park
 */

#include "TauFakeFactorAnalyzer.h"
#include "utility.h"

TauFakeFactorAnalyzer::TauFakeFactorAnalyzer(TTree *t, std::string outfilename, std::string year, std::string ch, std::string syst, std::string jsonfname, bool applytauFF, string globaltag, int nthreads, std::string mode)
:TopLFVAnalyzer(t, outfilename, year, ch, syst, jsonfname, applytauFF, globaltag, nthreads), _mode(mode)
{
    if (_ch.find("muon") != std::string::npos){
        cout << "muon channel" << endl;
    } else{
        cout << "electron channel" << endl;
    }

}

// Define your cuts here
void TauFakeFactorAnalyzer::defineCuts() {

    std::string lepCut = _isMuonCh ?
        "nmuonpass == 1 && nvetoelepass == 0 && nvetomuons == 0" :
        "nelepass == 1 && nvetoelepass == 0 && nvetomuons == 0";
    std::string chargeNameLoose = _isMuonCh ? "mutau_charge_loose" : "etau_charge_loose";
    std::string chargeNameTight = _isMuonCh ? "mutau_charge" : "etau_charge";

    if (_mode == "lss") { // loose tau vs jet && SS
        addCuts(lepCut + " && PV_npvsGood > 0 && nloosetaupass == 1 && ncleantaupass == 0 && " + chargeNameLoose + " > 0", "0");
        addCuts("ncleanjetsloosepass >= 3", "00");
        addCuts("ncleanbjetsloosepass == 1", "000");
    } else if (_mode == "los") { // loose tau vs jet && OS
        addCuts(lepCut + " && PV_npvsGood > 0 && nloosetaupass == 1 && ncleantaupass == 0 && " + chargeNameLoose + " < 0", "0");
        addCuts("ncleanjetsloosepass >= 3", "00");
        addCuts("ncleanbjetsloosepass == 1", "000");
    } else if (_mode == "tss") { // tight tau vs jet && SS
        addCuts(lepCut + " && PV_npvsGood > 0 && ncleantaupass == 1 && " + chargeNameTight + " > 0", "0");
        addCuts("ncleanjetspass >= 3", "00");
        addCuts("ncleanbjetspass == 1", "000");
    } else if (_mode == "tos") { // tight tau vs jet && OS
        addCuts(lepCut + " && PV_npvsGood > 0 && ncleantaupass == 1 && " + chargeNameTight + " < 0", "0");
        addCuts("ncleanjetspass >= 3", "00");
        addCuts("ncleanbjetspass == 1", "000");
    } else {
        std::cout << "WRONG MODE!!!!" << std::endl;
    }
}

void TauFakeFactorAnalyzer::defineMoreVars() {

    defineBTagNormalization();

    addVar({"Tau1_pt", "Tau_pt[0]", ""});
    addVar({"Tau1_pt_gen", "(Tau_pt_gen.size()>0) ? Tau_pt_gen[0] : -1", ""});
    addVar({"Tau1_charge", "Tau_charge[0]", ""});
    addVar({"Tau1_decayMode", "Tau_decayMode[0]", ""});
    addVar({"Tau1_decayMode_gen", "(Tau_pt_gen.size()>0) ? Tau_decayMode[0] : -1", ""});

    if (_isMuonCh) {
        addVar({"Muon1_charge", "Muon_charge[0]", ""});
        addVar({"mutau_charge", "Muon1_charge * Tau1_charge", ""});
        addVartoStore("Muon1.*");
    } else {
        addVar({"Electron1_charge", "Electron_charge[0]", ""});
        addVar({"etau_charge", "Electron1_charge * Tau1_charge", ""});
        addVartoStore("Electron1.*");
    }

    if (_mode == "lss" or _mode == "los") {
        addVar({"Tau1_pt_loose", "Tau_pt_loose[0]", ""});
        addVar({"Tau1_pt_loose_gen", "(Tau_pt_loose_gen.size()>0) ? Tau_pt_loose_gen[0] : -1", ""});
        addVar({"Tau1_charge_loose", "Tau_charge_loose[0]", ""});
        addVar({"Tau1_decayMode_loose", "Tau_decayMode_loose[0]", ""});
        addVar({"Tau1_decayMode_loose_gen", "(Tau_pt_loose_gen.size()>0) ? Tau_decayMode_loose[0] : -1", ""});
        if (_isMuonCh) {
            addVar({"mutau_charge_loose", "Muon1_charge * Tau1_charge_loose", ""});
        } else {
            addVar({"etau_charge_loose", "Electron1_charge * Tau1_charge_loose", ""});
        }
    }

    // EventWeights
    if (_syst == "data") {
        addVar({"eventWeight", "1.0"});
    } else {
        addVar({"eventWeight_genpu", "unitGenWeight * TopPtWeight[0] * puWeight[0]"});
        if (_isMuonCh) {
            addVar({"eventWeight_lep", "muonWeightId[0] * muonWeightIso[0] * muonWeightTrg[0]"});
        } else {
            addVar({"eventWeight_lep", "elecWeightReco[0] * elecWeightId[0] * elecWeightTrg[0]"});
        }
        if (_mode == "lss" or _mode == "los") {
            addVar({"eventWeight_tau", "tauWeightIdVsJet_loose[0][0] * tauWeightIdVsEl_loose[0][0] * tauWeightIdVsMu_loose[0][0]"});
        } else {
            addVar({"eventWeight_tau", "tauWeightIdVsJet[0][0] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
        }
        addVar({"eventWeight_nobtag", "eventWeight_genpu * eventWeight_lep * eventWeight_tau"});
        addVar({"eventWeight", "eventWeight_nobtag * btagWeightNorm[0]"});
    }

    // define variables that you want to store
    addVartoStore("run");
    addVartoStore("event");
    addVartoStore("Tau1.*");
    addVartoStore("eventWeight.*");
}

void TauFakeFactorAnalyzer::bookHists() {

    maxstep = "";

    add1DHist({"h_nevents", ";Number of events w/o b SF;Events", 2, -0.5, 1.5}, "one", "eventWeight", "", "0", "");

    if (_syst != "data") {
        //We anyway need this for bSF rescaling
        add1DHist({"h_nevents", ";Number of events w/o b SF;Events", 2, -0.5, 1.5}, "one", "eventWeight", "_nobtag", "0", "");
    }

    if (_mode == "lss" or _mode == "los") {
        add1DHist({"h_tau1_pt", ";#tau_{h} p_{T} (GeV);Events", 20, 0, 400}, "Tau1_pt_loose", "eventWeight", "", "0", maxstep);
        add1DHist({"h_tau1_gen_pt", ";#tau_{h} p_{T} (GeV) (MC: Tau_genPartFlav == 5);Events", 20, 0, 400}, "Tau1_pt_loose_gen", "eventWeight", "", "0", maxstep);
        add1DHist({"h_tau1_decayMode", ";Decaymode of #tau_{h};Events", 15, 0, 15}, "Tau1_decayMode_loose", "eventWeight", "", "0", maxstep);
        add1DHist({"h_tau1_decayMode_gen", ";Decaymode of #tau_{h};Events", 15, 0, 15}, "Tau1_decayMode_loose_gen", "eventWeight", "", "0", maxstep);
    } else {
        add1DHist({"h_tau1_pt", ";#tau_{h} p_{T} (GeV);Events", 20, 0, 400}, "Tau1_pt", "eventWeight", "", "0", maxstep);
        add1DHist({"h_tau1_gen_pt", ";#tau_{h} p_{T} (GeV) (MC: Tau_genPartFlav == 5);Events", 20, 0, 400}, "Tau1_pt_gen", "eventWeight", "", "0", maxstep);
        add1DHist({"h_tau1_decayMode", ";Decaymode of #tau_{h};Events", 15, 0, 15}, "Tau1_decayMode", "eventWeight", "", "0", maxstep);
        add1DHist({"h_tau1_decayMode_gen", ";Decaymode of #tau_{h};Events", 15, 0, 15}, "Tau1_decayMode_gen", "eventWeight", "", "0", maxstep);
    }
    add1DHist({"h_ncleanjetspass", ";Jet multiplicity;Events", 10, 0.0, 10.0}, "ncleanjetspass", "eventWeight", "", "0", maxstep);
    add1DHist({"h_ncleanbjetspass", ";b-tagged jet multiplicity;Events", 5, 0.0, 5.0}, "ncleanbjetspass", "eventWeight", "", "0", maxstep);
    add1DHist({"h_ncleanloosejetspass", ";Jet multiplicity;Events", 10, 0.0, 10.0}, "ncleanjetsloosepass", "eventWeight", "", "0", maxstep);
    add1DHist({"h_ncleanloosebjetspass", ";b-tagged jet multiplicity;Events", 5, 0.0, 5.0}, "ncleanbjetsloosepass", "eventWeight", "", "0", maxstep);
    add2DHist({"h_ncleanjetspass_ncleanloosejetspass", ";Jet multiplicity;Jet multiplicity (tau cleaned)", 10, 0.0, 10.0, 10, 0.0, 10.0}, "ncleanjetspass", "ncleanjetsloosepass", "eventWeight", "", "0", maxstep);
}

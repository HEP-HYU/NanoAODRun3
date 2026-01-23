/*
 * TauFakeFactorAnalyzer.cpp
 *
 *  Created on: June 14, 2023
 *      Author: Jiwon Park
 */

#include "TauFakeFactorAnalyzer.h"
#include "utility.h"

TauFakeFactorAnalyzer::TauFakeFactorAnalyzer(TTree *t, std::string outfilename, std::string year, std::string ch, std::string syst, std::string jsonfname, string globaltag, int nthreads, std::string mode)
:NanoAODAnalyzerrdframe(t, outfilename, year, ch, syst, jsonfname, globaltag, nthreads), _syst(syst), _year(year), _ch(ch), _mode(mode)
{
    _syst = syst;
}

// Define your cuts here
void TauFakeFactorAnalyzer::defineCuts() {

    if (_mode == "lss") { // loose tau vs jet && SS
        addCuts("nmuonpass == 1 && nvetoelepass == 0 && nvetomuons == 0 && PV_npvsGood > 0 && nloosetaupass == 1 && ncleantaupass == 0 && mutau_charge_loose > 0", "0");
        addCuts("ncleanjetsloosepass >= 3", "00");
        addCuts("ncleanbjetsloosepass == 1", "000");
        //addCuts("ncleanjetspass >= 3", "00");
        //addCuts("ncleanbjetspass == 1", "000");
    } else if (_mode == "los") { // loose tau vs jet && OS
        addCuts("nmuonpass == 1 && nvetoelepass == 0 && nvetomuons == 0 && PV_npvsGood > 0 && nloosetaupass == 1 && ncleantaupass == 0 && mutau_charge_loose < 0", "0");
        addCuts("ncleanjetsloosepass >= 3", "00");
        addCuts("ncleanbjetsloosepass == 1", "000");
        //addCuts("ncleanjetspass >= 3", "00");
        //addCuts("ncleanbjetspass == 1", "000");
    } else if (_mode == "tss") { // tight tau vs jet && SS
        addCuts("nmuonpass == 1 && nvetoelepass == 0 && nvetomuons == 0 && PV_npvsGood > 0 && ncleantaupass == 1 && mutau_charge > 0", "0");
        addCuts("ncleanjetspass >= 3", "00");
        addCuts("ncleanbjetspass == 1", "000");
    } else if (_mode == "tos") { // tight tau vs jet && OS
        addCuts("nmuonpass == 1 && nvetoelepass == 0 && nvetomuons == 0 && PV_npvsGood > 0 && ncleantaupass == 1 && mutau_charge < 0", "0");
        addCuts("ncleanjetspass >= 3", "00");
        addCuts("ncleanbjetspass == 1", "000");
    } else {
        std::cout << "WRONG MODE!!!!" << std::endl;
    }
}

void TauFakeFactorAnalyzer::defineMoreVars() {

    addVar({"Muon1_charge", "Muon_charge[0]", ""});

    addVar({"Tau1_pt", "Tau_pt[0]", ""});
    addVar({"Tau1_pt_gen", "(Tau_pt_gen.size()>0) ? Tau_pt_gen[0] : -1", ""});
    addVar({"Tau1_charge", "Tau_charge[0]", ""});
    addVar({"Tau1_decayMode", "Tau_decayMode[0]", ""});
    addVar({"Tau1_decayMode_gen", "(Tau_pt_gen.size()>0) ? Tau_decayMode[0] : -1", ""});
    addVar({"mutau_charge", "Muon1_charge * Tau1_charge", ""});

    if (_mode == "lss" or _mode == "los") {
        addVar({"Tau1_pt_loose", "Tau_pt_loose[0]", ""});
        addVar({"Tau1_pt_loose_gen", "(Tau_pt_loose_gen.size()>0) ? Tau_pt_loose_gen[0] : -1", ""});
        addVar({"Tau1_charge_loose", "Tau_charge_loose[0]", ""});
        addVar({"Tau1_decayMode_loose", "Tau_decayMode_loose[0]", ""});
        addVar({"Tau1_decayMode_loose_gen", "(Tau_pt_loose_gen.size()>0) ? Tau_decayMode_loose[0] : -1", ""});
        addVar({"mutau_charge_loose", "Muon1_charge * Tau1_charge_loose", ""});
    }

    // EventWeights
    if (_syst == "data") {
        addVar({"eventWeight", "1.0"});
    } else {
        addVar({"eventWeight_genpu", "unitGenWeight * TopPtWeight[0] * puWeight[0] * L1PreFiringWeight_Nom"});
        addVar({"eventWeight_mu", "muonWeightId[0] * muonWeightIso[0] * muonWeightTrg[0]"});
        if (_mode == "lss" or _mode == "los") {
            addVar({"eventWeight_tau", "tauWeightIdVsJet_loose[0][0] * tauWeightIdVsEl_loose[0][0] * tauWeightIdVsMu_loose[0][0]"});
            addVar({"eventWeight_nobtag", "eventWeight_genpu * eventWeight_mu * eventWeight_tau"});
            addVar({"eventWeight", "eventWeight_nobtag * btagWeight_PNetB_loose[0]"});
            //addVar({"eventWeight", "eventWeight_nobtag * btagWeight_PNetB[0]"});
        } else {
            addVar({"eventWeight_tau", "tauWeightIdVsJet[0][0] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight_nobtag", "eventWeight_genpu * eventWeight_mu * eventWeight_tau"});
            addVar({"eventWeight", "eventWeight_nobtag * btagWeight_PNetB[0]"});
        }
    }

    // define variables that you want to store
    addVartoStore("run");
    addVartoStore("event");
    addVartoStore("Tau1.*");
    addVartoStore("Muon1.*");
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

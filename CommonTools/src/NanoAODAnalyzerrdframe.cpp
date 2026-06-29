/*
 * NanoAODAnalyzerrdframe.cpp
 *
 *  Created on: Sep 30, 2018
 *      Author: suyong
 *      Refactored for LFV analysis: jiwon park
 */

#include "NanoAODAnalyzerrdframe.h"
#include <iostream>
#include <algorithm>
#include <typeinfo>
#include <random>
#include <chrono>
#include <ctime>
#include <sstream>

#include "TCanvas.h"
#include "Math/GenVector/VectorUtil.h"
#include <vector>
#include <fstream>
#include "utility.h"
#include <regex>
#include "ROOT/RDFHelpers.hxx"
#include "correction.h"
#include "GEScaleSyst.h"

using namespace std;

NanoAODAnalyzerrdframe::NanoAODAnalyzerrdframe(TTree *atree, std::string outfilename, std::string year, std::string ch, std::string syst, std::string jsonfname, std::string globaltag, int nthreads)
:_rd(*atree), _isData(false), _jsonOK(false), _outfilename(outfilename), _year(year), _ch(ch), _syst(syst), _jsonfname(jsonfname), _globaltag(globaltag), _inrootfile(0), _outrootfile(0), _rlm(_rd), _rnt(&_rlm), currentnode(0), PDFWeights(103, 0.0), PSWeights(4, 0.0), ScaleWeights(9, 0.0) {

    // record time
    auto start = std::chrono::system_clock::now();
    std::time_t start_time = std::chrono::system_clock::to_time_t(start);

    std::cout << "Start job on: " << std::ctime(&start_time) << std::endl;

    // Year switch
    if (_year.find("2022") != std::string::npos && _year.find("EE") == std::string::npos) {
        _isRun22 = true;
        cout << "Year : Run 2022" << endl;
    }
    if (_year.find("2022") != std::string::npos && _year.find("EE") != std::string::npos) {
        _isRun22EE = true;
        cout << "Year : Run 2022EE" << endl;
    }
    if (_year.find("2023") != std::string::npos && _year.find("BPix") == std::string::npos) {
        _isRun23 = true;
        cout << "Year : Run 2023" << endl;
    }
    if (_year.find("2023") != std::string::npos && _year.find("BPix") != std::string::npos) {
        _isRun23BPix = true;
        cout << "Year : Run 2023 BPix" << endl;
    }
    if (_year.find("2024") != std::string::npos){
        _isRun24 = true;
        cout << "Year : Run 2024" << endl;
    }

    // Data/mc switch
    if (_outfilename.find("TTto") == std::string::npos && atree->GetBranch("genWeight") == nullptr) {
        _isData = true;
        cout << "Input file is data " <<endl;
    } else {
        _isData = false;
        cout << "Input file is MC " << endl;
    }
    cout << "Outfile name: " << _outfilename << endl;

    // Systematic switch
    if (_syst != "") {
        cout << "Systematics Information" << endl;
        cout << "Activated option --syst " << _syst << endl;
    } else {
        cout << "Nominal process without systematics" << endl;
    }

    TObjArray *allbranches = atree->GetListOfBranches();
    for (int i =0; i<allbranches->GetSize(); i++) {
        TBranch *abranch = dynamic_cast<TBranch *>(allbranches->At(i));
        if (abranch!= nullptr) {
            std::string brname = abranch->GetName();
            if (brname.find("HLT_") == std::string::npos and brname.find("L1_") == std::string::npos)
                cout << brname << ", ";
            _originalvars.push_back(abranch->GetName());
        }
    }
}

NanoAODAnalyzerrdframe::~NanoAODAnalyzerrdframe() {
    // TODO Auto-generated destructor stub
    // ugly...
    std::cout<<">>  Job Done  <<"<<std::endl;
}

bool NanoAODAnalyzerrdframe::isDefined(string v) {
	auto result = std::find(_originalvars.begin(), _originalvars.end(), v);
	if (result != _originalvars.end()) return true;
	else return false;
}

void NanoAODAnalyzerrdframe::setTree(TTree *t, std::string outfilename) {
	_rd = ROOT::RDataFrame(*t);
	_rlm = RNode(_rd);
	_outfilename = outfilename;
	_hist1dinfovector.clear();
	_hist2dinfovector.clear();
	_th1dhistos.clear();
	_th2dhistos.clear();
	_varstostore.clear();
	_selections.clear();

	this->setupAnalysis();
}

void NanoAODAnalyzerrdframe::setupAnalysis() {
    if (_isData) _jsonOK = readjson();

    _rlm = _rlm.Define("one", "1.0")
               .Define("zero", "0.0");

    // Object selection will be defined in sequence.
    // Selected objects will be stored in new vectors.
    std::vector<std::string> jes_var;
    defineObjectSelection(jes_var);

    defineMoreVars();
    defineCuts();
    bookHists();
    setupCuts_and_Hists();
    setupTree();
}

bool NanoAODAnalyzerrdframe::readjson() {
    cout << "Applying Golden Json" << endl;

    auto isgoodjsonevent = [this](unsigned int runnumber, unsigned int lumisection) {

        auto key = std::to_string(runnumber);
        bool goodeventflag = false;

        if (jsonroot.isMember(key)) {
            Json::Value runlumiblocks = jsonroot[key];
            for (size_t i=0; i<runlumiblocks.size() && !goodeventflag; i++) {
                auto lumirange = runlumiblocks[i];
                if (lumisection >= lumirange[0].asUInt() && lumisection <= lumirange[1].asUInt()) goodeventflag = true;
            }
            return goodeventflag;
        } else {
            //cout << "Run not in json " << runnumber << endl;
            return false;
        }

    };

    if (_jsonfname != "") {
        std::ifstream jsoninfile;
        jsoninfile.open(_jsonfname);

        if (jsoninfile.good()) {
            jsoninfile >> jsonroot;
            _rlm = _rlm.Define("goodjsonevent", isgoodjsonevent, {"run", "luminosityBlock"}).Filter("goodjsonevent");
            _jsonOK = true;
            return true;
        } else {
            cout << "Problem reading json file " << _jsonfname << endl;
            return false;
        }
    } else {
        cout << "no JSON file given" << endl;
        return true;
    }
}

void NanoAODAnalyzerrdframe::defineObjectSelection(std::vector<std::string> jes_var){
    //Override at SkimEvents.cpp and analysisAnalyzer.cpp
    std::string muoncut  = "Muon_pt>50.0 && abs(Muon_eta)<2.4 && Muon_tightId && Muon_pfRelIso04_all<0.15";
    std::string vetomuon = "!muoncuts && Muon_pt>15.0 && abs(Muon_eta)<2.4 && Muon_looseId && Muon_pfRelIso04_all<0.25";
    std::string skimjet = "Jet_pt>30.0 && abs(Jet_eta)<2.6 && (Jet_passJetIdTightLepVeto>=0.0) && Jet_muEF<0.8 && Jet_chEmEF<0.8";
    //std::string skimjet = "Jet_pt>30.0 && abs(Jet_eta)<2.6 && Jet_muEF<0.8 && Jet_chEmEF<0.8";
    selectMuons(muoncut, vetomuon);
    skimJets(skimjet);
}

void NanoAODAnalyzerrdframe::noiseFilter(){
    auto rerun_ecalBadCalibFilter = [this] (unsigned int run, float metpt, float metphi, floats jetpts, floats jetetas, floats jetphis, floats neEmEF, floats chEmEF)->bool{
        if (run < 362433u || run > 367144u) return true;
        if (metpt < 100.) return true;
        for (size_t i=0; i<jetpts.size(); i++){
            if (jetpts[i] <= 50.) continue;
            if (jetetas[i] <= -0.5 || jetetas[i] >= -0.1) continue;
            if (jetphis[i] <= -2.1 || jetphis[i] >= -1.8) continue;
            if (neEmEF[i] + chEmEF[i] <= 0.9) continue;
            float dphi_tmp = metphi-jetphis[i];
            float dphi = std::atan2(std::sin(dphi_tmp), std::cos(dphi_tmp));
            if(std::abs(dphi) < 2.9) continue;
            return false;
        }
        return true;
    };
    _rlm = _rlm.Define("rerun_ecalBadCalibFilter", rerun_ecalBadCalibFilter, {"run", "PuppiMET_pt", "PuppiMET_phi", "Jet_pt", "Jet_eta", "Jet_phi", "Jet_neEmEF", "Jet_chEmEF"})
               .Redefine("Flag_ecalBadCalibFilter", "Flag_ecalBadCalibFilter && rerun_ecalBadCalibFilter");
}


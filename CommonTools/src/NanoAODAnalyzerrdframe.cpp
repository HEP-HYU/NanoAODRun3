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
    if (_outfilename.find("TTto") == std::string::npos && atree->GetBranch("genWeight") == nullptr || _outfilename.find("TT_LFV") != std::string::npos || _outfilename.find("ST_LFV") != std::string::npos) {
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

void NanoAODAnalyzerrdframe::applyWeights(string pileFile, string map){
    // Event weight for data it's always one. For MC, it depends on the sign
    _rlm = _rlm.Define("lhereweight","one");
    _rlm = _rlm.Define("unitGenWeight", "one");

    _rlm = _rlm.Define("isData", "true");
    if(!_isData){
        _rlm = _rlm.Redefine("isData", "false");

        if (_outfilename.find("WtoLNu") != std::string::npos) {
          _rlm = _rlm.Redefine("lhereweight","LHEWeight_originalXWGTUP/abs(LHEWeight_originalXWGTUP)");
          _rlm = _rlm.Redefine("unitGenWeight","LHEWeight_originalXWGTUP/abs(LHEWeight_originalXWGTUP)");
        };


        // Store sum of weights
        auto storePDFWeights = [this](floats weights, float gen)->floats {
            for (unsigned int i=0; i<weights.size(); i++)
                PDFWeights[i] += (gen / abs(gen)) * weights[i];
            return PDFWeights;
        };
        auto storePSWeights = [this](floats weights, float gen)->floats {
            for (unsigned int i=0; i<weights.size(); i++) {
                if (i > 3) continue; //JME Nano stores all PS
                PSWeights[i] += (gen / abs(gen)) * weights[i];
            }
            return PSWeights;
        };
        auto storeScaleWeights = [this](floats weights, float gen)->floats {
            for (unsigned int i=0; i<weights.size(); i++)
                ScaleWeights[i] += (gen / abs(gen)) * weights[i];
            return ScaleWeights;
        };
        try {
            _rlm.Foreach(storePDFWeights, {"LHEPdfWeight", "genWeight"});
        } catch (exception& e) {
            cout << e.what() << endl;
            cout << "No PDF weight in this root file!" << endl;
        }
        try {
            _rlm.Foreach(storePSWeights, {"PSWeight", "genWeight"});
        } catch (exception& e) {
            cout << e.what() << endl;
            cout << "No PS weight in this root file!" << endl;
        }
        try {
            _rlm.Foreach(storeScaleWeights, {"LHEScaleWeight", "genWeight"});
        } catch (exception& e) {
            cout << e.what() << endl;
            cout << "No Scale weight in this root file!" << endl;
        }

        ////Check Normalisation issue for genWeight
        _rlm = _rlm.Redefine("unitGenWeight","genWeight != 0 ? genWeight/abs(genWeight) : 0");
        
        auto puWeightreader = correction::CorrectionSet::from_file("data/LUM/"+pileFile+"/puWeights.json.gz");
        auto _puweight = puWeightreader->at(map);

        auto PuWeight = [this, _puweight](float x) -> floats {
            floats out;
            out.emplace_back(_puweight->evaluate({x, "nominal"}));
            out.emplace_back(_puweight->evaluate({x, "up"}));
            out.emplace_back(_puweight->evaluate({x, "down"}));

            return out;
        };

        _rlm = _rlm.Define("puWeight", PuWeight, {"Pileup_nTrueInt"});
    }
}

void NanoAODAnalyzerrdframe::defineObjectSelection(std::vector<std::string> jes_var){
    //Override at SkimEvents.cpp and analysisAnalyzer.cpp
    std::string muoncut  = "Muon_pt>50.0 && abs(Muon_eta)<2.4 && Muon_tightId && Muon_pfRelIso04_all<0.15";
    std::string vetomuon = "!muoncuts && Muon_pt>15.0 && abs(Muon_eta)<2.4 && Muon_looseId && Muon_pfRelIso04_all<0.25";
    std::string skimjet = "Jet_pt>30.0 && abs(Jet_eta)<2.6 && (Jet_JetId>=0.0) && Jet_muEF<0.8 && Jet_chEmEF<0.8";
    //std::string skimjet = "Jet_pt>30.0 && abs(Jet_eta)<2.6 && Jet_muEF<0.8 && Jet_chEmEF<0.8";
    selectMuons(muoncut, vetomuon);
    skimJets(skimjet);
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
            for (unsigned int i=0; i<runlumiblocks.size() && !goodeventflag; i++) {
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

void NanoAODAnalyzerrdframe::JetVetoMap(string jetFile, string map) {
    cout << "apply JetVetoMap" << endl;
    
    auto checkoverlap = [](FourVectorVec &seljets, FourVectorVec &sellep) {
        ints mindrlepton;
        for (auto ajet: seljets) {
            auto mindr = 6.0;
            for ( auto alepton : sellep ) {
                auto dr = ROOT::Math::VectorUtil::DeltaR(ajet, alepton);
                if (dr < mindr) mindr = dr;
            }
            int out = mindr >= 0.2 ? 1 : 0;
            mindrlepton.emplace_back(out);
        }
        return mindrlepton;
    };

    //TODO: _isRun24 is implementation for NanoAODv15
    if (_isRun24){
        auto jetIDreader = correction::CorrectionSet::from_file("data/JME/2024_Summer24/jetid.json.gz");
        auto _jetid = jetIDreader->at("AK4PUPPI_TightLeptonVeto");
        auto getJetId1 = [this, _jetid](floats &etas, floats &chHEF, floats &neHEF, floats &chEmEF, floats &neEmEF, floats &muEF, uchars &Nch, uchars &Nne)->floats{
            floats jetIds;
            jetIds.reserve(etas.size());
            for (int i=0; i<etas.size(); i++){
                float jetId = _jetid->evaluate({float(etas[i]), float(chHEF[i]), float(neHEF[i]), float(chEmEF[i]), float(neEmEF[i]), float(muEF[i]), int(Nch[i]), int(Nne[i]), int(Nch[i])+int(Nne[i])});
                jetIds.emplace_back(jetId);
            }
            return jetIds;
        };
        _rlm = _rlm.Define("Jet_JetId", getJetId1, {"Jet_eta", "Jet_chHEF", "Jet_neHEF", "Jet_chEmEF", "Jet_neEmEF", "Jet_muEF", "Jet_chMultiplicity", "Jet_neMultiplicity"});
    } 
    //TODO: below is for NanoAODv12
    //If the others are going to v15, this should be removed
    else{
        auto getJetId2= [this](floats &etas, uchars &jetIds, floats &muEF, floats &chEmEF, floats &neHEF, floats neEmEF)->floats{
            floats idVec;
            idVec.reserve(etas.size());
            for (int i=0; i<etas.size(); i++){
                float id;
                if (std::abs(etas[i]) <= 2.7) id  = (jetIds[i] & (1 << 1)) && (muEF[i] < 0.8) && (chEmEF[i] < 0.8) ;
                else if (std::abs(etas[i]) > 2.7 && std::abs(etas[i]) <= 3.0) id = (jetIds[i] & (1 << 1)) && (neHEF[i] < 0.99);
                else if (std::abs(etas[i]) > 3.0) id = (jetIds[i] & (1 << 1)) && (neEmEF[i] < 0.4);
                idVec.emplace_back(id);
            }
            return idVec;
        };
        _rlm = _rlm.Define("Jet_JetId", getJetId2, {"Jet_eta", "Jet_jetId", "Jet_muEF", "Jet_chEmEF", "Jet_neHEF", "Jet_neEmEF"});
    }

    
    _rlm = _rlm.Define("loosejetcuts", "Jet_pt>15 && (Jet_JetId==1.0) && (Jet_neEmEF+Jet_chEmEF)<0.9");
    
    _rlm = _rlm.Define("Jet_pt_loosejet", "Jet_pt[loosejetcuts]")
               .Define("Jet_phi_loosejet", "Jet_phi[loosejetcuts]")
               .Define("Jet_eta_loosejet", "Jet_eta[loosejetcuts]")
               .Define("Jet_mass_loosejet", "Jet_mass[loosejetcuts]")
               .Define("njetloose", "int(Jet_pt_loosejet.size())")
               .Define("loosejet4vecs", ::gen4vec, {"Jet_pt_loosejet", "Jet_eta_loosejet", "Jet_phi_loosejet", "Jet_mass_loosejet"});
    
    auto vetoMapreader = correction::CorrectionSet::from_file("data/JME/"+jetFile+"/jetvetomaps.json.gz");
    
    auto _vetomap = vetoMapreader->at(map);
 
    auto vetomap = [this, _vetomap](floats &eta, floats &phi)->floats {
        floats xout;
        for (unsigned int i=0; i<eta.size(); i++) {
            float es = 0.0;
            es = _vetomap->evaluate({std::string("jetvetomap"),eta[i],phi[i]});
            xout.emplace_back(es);
        }
        return xout;
    };
    
    auto vetomap_fpix = [this, _vetomap](floats &eta, floats &phi)->floats {
        floats xout;
        for (unsigned int i=0; i<eta.size(); i++) {
            float es = 0.0;
            es = _vetomap->evaluate({std::string("jetvetomap_fpix"),eta[i],phi[i]});
            xout.emplace_back(es);
        }
        return xout;
    };
    
    _rlm = _rlm.Define("Jet_isVeto_loose", vetomap, {"Jet_eta_loosejet","Jet_phi_loosejet"})
               .Define("events_isVeto","Sum(Jet_isVeto_loose)");

    if (_isRun24){
        _rlm = _rlm.Define("Jet_isVeto_fpix", vetomap_fpix, {"Jet_eta_loosejet", "Jet_phi_loosejet"})
                   .Redefine("events_isVeto", "Sum(Jet_isVeto_loose)+Sum(Jet_isVeto_fpix)");
    }
}

void NanoAODAnalyzerrdframe::selectElectrons(string cut, string vetocut) {
    if (cut.find("onlyveto") != std::string::npos){
        cout<<"selectElectrons: only vetoelecuts"<<endl;
        _rlm = _rlm.Define("vetoelecuts", vetocut)
                   .Define("nvetoelepass","Sum(vetoelecuts)");
    } else{
        cout<<"selectElectrons: electroncuts, vetoelecuts"<<endl;
        _rlm = _rlm.Define("elecuts", cut)               
                   .Define("vetoelecuts", vetocut);

        _rlm = _rlm.Define("nvetoelepass","Sum(vetoelecuts)")
                   .Redefine("Electron_pt", "Electron_pt[elecuts]")
                   .Redefine("Electron_eta", "Electron_eta[elecuts]")
                   .Redefine("Electron_phi", "Electron_phi[elecuts]")
                   .Redefine("Electron_mass", "Electron_mass[elecuts]")
                   .Redefine("Electron_charge", "Electron_charge[elecuts]")
                   //.Define("Electron_Id", "Electron_looseId[elecuts]")
                   .Redefine("Electron_pfRelIso03_all", "Electron_pfRelIso03_all[elecuts]")
                   .Define("Sel_eleidx", ::good_idx, {"elecuts"})
                   .Define("nelepass", "int(Electron_pt.size())")
                   .Define("lep4vecs", ::gen4vec, {"Electron_pt", "Electron_eta", "Electron_phi", "Electron_mass"});
    }
}

void NanoAODAnalyzerrdframe::selectMuons(string cut, string vetocut) {
    if (cut.find("onlyveto") != std::string::npos){
        cout<<"selectMuons: only vetomuoncuts"<<endl;
        _rlm = _rlm.Define("vetomuoncuts", vetocut)
                   .Define("nvetomuons", "Sum(vetomuoncuts)");
    } else{
        cout<<"selectMuons: muoncuts, vetomuoncuts"<<endl;
        _rlm = _rlm.Define("muoncuts", cut);
        _rlm = _rlm.Define("vetomuoncuts", vetocut);

        _rlm = _rlm.Define("nvetomuons","Sum(vetomuoncuts)")
                   .Redefine("Muon_pt", "Muon_pt[muoncuts]")
                   .Redefine("Muon_eta", "Muon_eta[muoncuts]")
                   .Redefine("Muon_phi", "Muon_phi[muoncuts]")
                   .Redefine("Muon_mass", "Muon_mass[muoncuts]")
                   .Redefine("Muon_charge", "Muon_charge[muoncuts]")
                   .Redefine("Muon_looseId", "Muon_looseId[muoncuts]")
                   .Redefine("Muon_pfRelIso04_all", "Muon_pfRelIso04_all[muoncuts]")
                   .Define("Sel_muonidx", ::good_idx, {"muoncuts"})
                   .Define("nmuonpass", "int(Muon_pt.size())")
                   .Define("lep4vecs", ::gen4vec, {"Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass"});
    }
}

// Exactly one muon case
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

    auto muonSFreader = correction::CorrectionSet::from_file("data/MuonSF/ScaleFactors_Muon_Z_ID_ISO_"+muonFile+"_schemaV2.json");
    auto _muonid = muonSFreader->at(muonid);
    auto _muoniso = muonSFreader->at(muoniso);
    if (_isRun24) muonFile = "2023_BPix";
    auto muonHltSFreader = correction::CorrectionSet::from_file("data/MuonSF/ScaleFactors_Muon_Z_HLT_"+muonFile+"_eta_pt_schemaV2.json");
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
    
    std::string muonYear = _year;

    auto muonhighscaleup = [muonYear](floats &pts, floats &etas, floats &phis, ints &charges)->floats {
        floats out;
        out.reserve(pts.size());
        for (unsigned int i=0; i<pts.size(); i++) {
            float pt_tmp = pts[i];
            float pt_out = pt_tmp;
            if (pt_tmp > 200) {
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

    auto muonhighscaledn = [muonYear](floats &pts, floats &etas, floats &phis, ints &charges)->floats {
        floats out;
        out.reserve(pts.size());
        for (unsigned int i=0; i<pts.size(); i++) {
            float pt_tmp = pts[i];
            float pt_out = pt_tmp;
            if (pt_tmp > 200) {
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
        for (unsigned int i=0; i<pts.size(); i++) {
            out -= ptcors[i] - pts[i];
        }
        return out;
    };

    auto muonhighscalemetphi = [](floats &pts, floats &ptcors, floats &phis, float met, float metphi)->float {
        float out = 0.;
        auto metx = met * cos(metphi);
        auto mety = met * sin(metphi);
        for (unsigned int i=0; i<pts.size(); i++) {
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

    auto elecSFreader = correction::CorrectionSet::from_file("data/ElectronSF/"+elecFile+"/electron.json.gz");
    auto _elecid = elecSFreader->at("Electron-ID-SF");

    auto elecSFId = [this, _elecid, elecYear](std::string var){
        return [this, _elecid, elecYear, var](floats &pt, floats &eta, floats &phi)->floats {
            floats wVec;
            wVec.reserve(3); //cent, up, down
            std::string wp = var;
            if (pt.size() == 1){
                float sf=1.0; float sf_up=1.0; float sf_down=1.0;
                if (elecYear.find("2022")!=std::string::npos || elecYear.find("2024")!=std::string::npos){
                    if (pt[0] >= 20 && pt[0] < 75 && wp.find("Reco")!=std::string::npos) wp = "Reco20to75";
                    sf = _elecid->evaluate({elecYear, "sf", wp, eta[0], pt[0]});
                    sf_up = _elecid->evaluate({elecYear, "sfup", wp, eta[0], pt[0]});
                    sf_down = _elecid->evaluate({elecYear, "sfdown", wp, eta[0], pt[0]});
                } else if (elecYear.find("2023")!=std::string::npos){
                    if (pt[0] >= 20 && pt[0] < 75 && wp.find("Reco")!=std::string::npos) wp = "Reco20to75";
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

    auto elecHltSFreader = correction::CorrectionSet::from_file("data/ElectronSF/"+elecFile+"/electronHlt.json.gz");
    auto _elechlt = elecHltSFreader->at("Electron-HLT-SF");

    auto elecSFTrg = [this, _elechlt, elecYear](floats &pt, floats &eta)->floats {
        floats wVec;
        wVec.reserve(3); //cent, up, down
        if (pt.size() == 1){
            float sf = _elechlt->evaluate({elecYear, "sf", "HLT_SF_Ele30_MVAiso90ID", eta[0], pt[0]});
            float sf_up = _elechlt->evaluate({elecYear, "sfup", "HLT_SF_Ele30_MVAiso90ID", eta[0], pt[0]});
            float sf_down = _elechlt->evaluate({elecYear, "sfdown", "HLT_SF_Ele30_MVAiso90ID", eta[0], pt[0]});
            wVec.emplace_back(sf);
            wVec.emplace_back(sf_up);
            wVec.emplace_back(sf_down);
        } else wVec = {1.0, 1.0, 1.0};
        return wVec;
    };

    _rlm = _rlm.Define("elecWeightTrg", elecSFTrg, {"Electron_pt","Electron_eta"});

}

void NanoAODAnalyzerrdframe::setupJetMETCorrection(string jecFile, string jecYear, string jerMap, bool dataMc){
    auto jercReader = correction::CorrectionSet::from_file("data/JME/"+jecFile+"/jet_jerc.json.gz");
    string datamcflag = "";
    if (dataMc) datamcflag = "DATA";
    else datamcflag = "MC";

    string tmp = jecYear + "_V2_" + datamcflag + "_L1FastJet_AK4PFPuppi";
    std::cout << tmp << std::endl;
    auto _L1FastJet = jercReader->at(jecYear + "_V2_" + datamcflag + "_L1FastJet_AK4PFPuppi");
    auto _L2Relative = jercReader->at(jecYear + "_V2_" + datamcflag + "_L2Relative_AK4PFPuppi");
    auto _L2L3Residual = jercReader->at(jecYear + "_V2_DATA_L2L3Residual_AK4PFPuppi");

    auto applyJes = [this, _L1FastJet, _L2Relative, _L2L3Residual, dataMc](floats jetpts, floats jetetas, floats jetphis, floats jetAreas, floats jetrawf, float rho, unsigned int run, floats toCorr)->floats {
        floats corrfactors;
        corrfactors.reserve(jetpts.size());

        for (unsigned int i=0; i<jetpts.size(); i++) {
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

        for (unsigned int i=0; i<jetpts.size(); i++){
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

void NanoAODAnalyzerrdframe::skimJets(string cut) {
    // input vector: vec[pt][vars]
    // Note: do not skim with exact value of pt!
    auto skimCol = [this](floatsVec toSkim, ints cut)->floatsVec {

        floatsVec out;
        for (size_t i=0; i<toSkim.size(); i++) {
            if (cut[i] > 0) out.emplace_back(toSkim[i]);
        }
        return out;
    };

    // skim jet collection
    // https://twiki.cern.ch/twiki/bin/view/CMS/JetID13p6TeV
    // nanoAOD Flags
    _rlm = _rlm.Define("jetcuts", cut)
               .Redefine("Jet_pt", "Jet_pt[jetcuts]")
               .Redefine("Jet_eta", "Jet_eta[jetcuts]")
               .Redefine("Jet_phi", "Jet_phi[jetcuts]")
               .Redefine("Jet_mass", "Jet_mass[jetcuts]")
               .Redefine("Jet_JetId", "Jet_JetId[jetcuts]")
               .Redefine("Jet_area", "Jet_area[jetcuts]")
               .Redefine("Jet_rawFactor", "Jet_rawFactor[jetcuts]")
               .Redefine("Jet_pt_uncorr", "Jet_pt_uncorr[jetcuts]")
               .Redefine("Jet_btagPNetB", "Jet_btagPNetB[jetcuts]")
               .Redefine("nJet", "int(Jet_pt.size())");
    if (_isRun24){
        _rlm = _rlm.Redefine("Jet_btagUParTAK4B", "Jet_btagUParTAK4B[jetcuts]");
    }
    if (!_isData) {
        //_rlm = _rlm.Redefine("Jet_pt_unc", skimCol, {"Jet_pt_unc", "jetcuts"})
        //           .Redefine("Jet_jer", skimCol, {"Jet_jer", "jetcuts"})
        _rlm = _rlm.Redefine("Jet_hadronFlavour","Jet_hadronFlavour[jetcuts]")
                   .Redefine("Jet_genJetIdx","Jet_genJetIdx[jetcuts]");
    }
}

void NanoAODAnalyzerrdframe::applyBSFs(std::vector<string> jes_var, string btagYear, string btagMap) {
    cout << "Loading Btag SF: " << btagYear << endl;

    _rlm = _rlm.Define("btag_var", [](){return strings(btag_var);})
               .Define("btag_jes_var", [jes_var](){return strings(jes_var);});

    //case3 - Shape correction
    //If you are interested in using the whole b-tagging discriminant distribution in your analysis,
    //e.g. using b-tagging variables to separate signal and background, then this method is for you
    auto bSFreader = correction::CorrectionSet::from_file("data/BTV/" + btagYear + "/btagging.json.gz");
    auto _btagSF = bSFreader->at(btagMap);
    auto btagSF_shape = [this, _btagSF](floats &pts, floats &etas, uchars &hadflav, floats &btags)->floatsVec{
        floats wVec;
        wVec.reserve(17);
        floatsVec out;
        out.reserve(pts.size());

        std::vector<std::string> systs = {"central",
            "up_hf", "down_hf", "up_lf", "down_lf",
            "up_hfstats1", "down_hfstats1", "up_hfstats2", "down_hfstats2",
            "up_lfstats1", "down_lfstats1", "up_lfstats2", "down_lfstats2",
            "up_cferr1", "down_cferr1", "up_cferr2", "down_cferr2"};
        
        for (auto i=0; i<int(pts.size()); i++){
            for (auto &syst: systs){
                float sf = 1.0;
                if (pts[i] <= 40) sf = 1.0;
                else if (syst.find("cferr")!=std::string::npos && hadflav[i]!=4) sf = 1.0;
                else if ((syst.find("hf")!=std::string::npos || syst.find("lf")!=std::string::npos) && hadflav[i]==4) sf = 1.0;
                else sf = _btagSF->evaluate({syst, int(hadflav[i]), abs(etas[i]), pts[i], btags[i]}); 
                std::cout << i << " " << syst << " " << sf << std::endl;
                wVec.emplace_back(sf);
            }
            out.emplace_back(wVec);
            wVec.clear();
        }
        return out;
    };
    auto btag_evWeight = [this](floatsVec &btagWeights) ->floatsVec{
        const int vars = 17;
        floats out(vars, 1.0);

        for (auto &jet: btagWeights){
            for (int i=0; i<vars; i++){
                out[i] *= jet[i];
            }
        }
        return out;
    }
    _rlm = _rlm.Define("Jet_btagSF", btagSF_shape, {"Jet_pt", "Jet_eta", "Jet_hadronFlavour", "Jet_btagPNetB"})
               .Define("btagWeight", btag_evWeight, {"Jet_btagSF"});

    //auto btagSF24_fixWP = [this, _btagSF](floats &pts, floats &etas, uchars &hadflav, floats &btags)->floats{
    //    cout << " btag fix wp" << endl;
    //    floats wVec;
    //    wVec.reserve(15);

    //    std::vector<std::string> systs = {"central", "up", "down",
    //        "up_correlated", "down_correlated", "up_uncorrelated", "down_uncorrelated",
    //        "up_jes", "down_jes", "up_jer", "down_jer",
    //        "up_pileup", "down_pileup", "up_statistic", "down_statistic"};
    //    
    //    if (int(pts.size()) != int(etas.size())) cout << "eta size hmmmmmmmmmmmm" << endl;
    //    if (int(pts.size()) != int(hadflav.size())) cout << "hadflav size hmmmmmmmmmmmm" << endl;
    //    if (int(pts.size()) != int(btags.size())) cout << "btag size hmmmmmmmmmmmm" << endl;
    //    if (std::min({etas.size(), hadflav.size(), btags.size()}) != pts.size()) return wVec;
    //    for (auto &syst: systs){
    //        float w = 1.0;
    //        for (auto i=0; i<int(pts.size()); i++){
    //            if (btags[i] > 0.1272) { // WP M: 0.1272
    //                try{
    //                    float sf = _btagSF->evaluate({syst, "M", int(hadflav[i]), abs(etas[i]), pts[i]});
    //                    w = w*sf;
    //                } catch (const std::exception &e){
    //                    cout << "e: " << e.what() << endl;
    //                    cout << "syst: " << syst << " flavor: " << hadflav[i] << " eta: " << etas[i] << " pt: " << pts[i] << endl;
    //                }
    //            }
    //        }
    //        wVec.emplace_back(w);
    //    }
    //    return wVec;
    //};

    //if (false){
    //    _rlm = _rlm.Define("btagWeight", btagSF24_fixWP, {"Jet_pt", "Jet_eta", "Jet_hadronFlavour", "Jet_btagUParTAK4B"});
    
    //TODO: apply b-tagging SF
    //cout<<"Generate case3 b-tagging weight"<<endl;
    //std::string column_name = output_var + "case3";
    //_rlm = _rlm.Define(column_name, btagweightgenerator3, Jets_vars_names);
    //Total event weight after shape correction
    //_rlm = _rlm.Define("evWeight", "pugenWeight*btagWeight_case3");
    //std::cout<< "BJet SF column name: " << column_name << std::endl;

    //cout << "Generate b-tagging weight" << endl;
    //_rlm = _rlm.Define("btagWeight_PNetB_perJet", btagweightgenerator, {"Jet_pt", "Jet_eta", "Jet_hadronFlavour", "Jet_btagPNetB"})
    //           .Define("btagWeight_PNetB_jes_perJet", btagweightgenerator, {"Jet_pt_unc", "Jet_eta", "Jet_hadronFlavour", "Jet_btagPNetB"});
}

void NanoAODAnalyzerrdframe::selectJets(std::vector<std::string> jes_var, std::vector<std::string> jes_var_flav, std::string cut) {
    //// input vector: vec[pt][vars], for bSF
    auto skimCol = [this](doublesVec toSkim, ints cut)->doublesVec {

        doublesVec out;
        for (size_t i=0; i<toSkim.size(); i++) {
            if (cut[i] > 0) out.emplace_back(toSkim[i]);
        }
        return out;
    };

    //// input vector: vec[pt][vars]
    auto calcBSF = [this](doublesVec perJetSF, int nvar)->doubles {

        doubles out;
        out.reserve(nvar);
        for (size_t i=0; i<nvar; i++) {
            double bSF = 1.0;
            for (size_t j=0; j<perJetSF.size(); j++) {
                if (perJetSF[j].empty()) continue;
                bSF *= perJetSF[j][i];
            }
            out.emplace_back(bSF);
        }
        return out;
    };

    //if (!_isData) {

    //    auto syst_unc = _syst;

    //    auto selectJer = [syst_unc](floatsVec unc)->floats {

    //        int idx = 0;
    //        if (syst_unc == "jerup") idx = 1;
    //        else if (syst_unc == "jerdown") idx = 2;
    //        floats selected;
    //        selected.reserve(unc.size());

    //        for (size_t i=0; i<unc.size(); i++) {
    //            selected.emplace_back(unc[i][idx]);
    //        }
    //        return selected;
    //    };

    //    if (_syst.find("jes") != std::string::npos) {

    //        auto selectJes = [syst_unc, jes_var, jes_var_flav](floatsVec unc)->floats {

    //            floats selected;
    //            selected.reserve(unc.size());

    //            std::vector<std::string> jes_var_all = jes_var;
    //            jes_var_all.insert(jes_var_all.end(), jes_var_flav.begin(), jes_var_flav.end());
    //            unsigned int jesidx = -1;
    //            for (size_t i=0; i<jes_var_all.size(); i++) {
    //                if (jes_var_all[i] == syst_unc) jesidx = i;
    //            }
    //            if (int(jesidx) == -1) cerr << "Found No JES Unc Name!!" << endl;

    //            for (size_t i=0; i<unc.size(); i++) {
    //                selected.emplace_back(unc[i][jesidx]);
    //            }
    //            return selected;

    //        };

    //        _rlm = _rlm.Define("Jet_pt_unc_toapply", selectJes, {"Jet_pt_unc"})
    //                   .Define("Jet_jer_toapply", selectJer, {"Jet_jer"})
    //                   .Redefine("Jet_pt", "Jet_pt * Jet_jer_toapply * Jet_pt_unc_toapply")
    //                   .Redefine("Jet_mass", "Jet_mass * Jet_jer_toapply * Jet_pt_unc_toapply");

    //    } else {
    //        _rlm = _rlm.Define("Jet_jer_toapply", selectJer, {"Jet_jer"})
    //                   .Redefine("Jet_pt", "Jet_pt * Jet_jer_toapply")
    //                   .Redefine("Jet_mass", "Jet_mass * Jet_jer_toapply");
    //    }

    //    if (_syst.find("metUnclust") != std::string::npos) {

    //        if (_syst.find("up") != std::string::npos) {
    //            _rlm = _rlm.Redefine("MET_pt", "float (sqrt(pow(MET_pt * cos(MET_phi) + MET_MetUnclustEnUpDeltaX, 2) + pow(MET_pt * sin(MET_phi) + MET_MetUnclustEnUpDeltaY, 2)) )")
    //                       .Redefine("MET_phi", "float (atan2((MET_pt * sin(MET_phi) + MET_MetUnclustEnUpDeltaY), (MET_pt * cos(MET_phi) + MET_MetUnclustEnUpDeltaX)) )");
    //        } else if (_syst.find("down") != std::string::npos) {
    //            _rlm = _rlm.Redefine("MET_pt", "float (sqrt(pow(MET_pt * cos(MET_phi) - MET_MetUnclustEnUpDeltaX, 2) + pow(MET_pt * sin(MET_phi) - MET_MetUnclustEnUpDeltaY, 2)) )")
    //                       .Redefine("MET_phi", "float (atan2((MET_pt * sin(MET_phi) - MET_MetUnclustEnUpDeltaY), (MET_pt * cos(MET_phi) - MET_MetUnclustEnUpDeltaX)) )");
    //        }
    //    }
    //}

    //https://twiki.cern.ch/twiki/bin/view/CMS/JetID13p6TeV
    //nanoAOD Flags
    _rlm = _rlm.Define("jetcuts", cut);
    _rlm = _rlm.Redefine("Jet_pt", "Jet_pt[jetcuts]")
               .Redefine("Jet_eta", "Jet_eta[jetcuts]")
               .Redefine("Jet_phi", "Jet_phi[jetcuts]")
               .Redefine("Jet_mass", "Jet_mass[jetcuts]")
               .Redefine("Jet_btagPNetB", "Jet_btagPNetB[jetcuts]")
               .Define("jet4vecs", ::gen4vec, {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass"});
    if (_isRun24){
        _rlm = _rlm.Redefine("Jet_btagUParTAK4B", "Jet_btagUParTAK4B[jetcuts]");
    }

    //if (!_isData) {
    //    _rlm = _rlm.Redefine("btagWeight_PNetB_perJet", skimCol, {"btagWeight_PNetB_perJet", "jetcuts"})
    //               .Redefine("btagWeight_PNetB_jes_perJet", skimCol, {"btagWeight_PNetB_jes_perJet", "jetcuts"});
    //}

    // for checking overlapped jets with leptons
    auto checkoverlap = [](FourVectorVec &seljets, FourVectorVec &sellep) {

        ints mindrlepton;
        for (auto ajet: seljets) {
            auto mindr = 6.0;
            for ( auto alepton : sellep ) {
                auto dr = ROOT::Math::VectorUtil::DeltaR(ajet, alepton);
                if (dr < mindr) mindr = dr;
            }
            int out = mindr >= 0.4 ? 1 : 0;
            mindrlepton.emplace_back(out);
        }
        return mindrlepton;
    };

    // Overlap removal with muon / electron (used for btagging SF)
    if (isDefined("cleantau4vecs")){
        _rlm = _rlm.Define("lepjetoverlap", checkoverlap, {"jet4vecs","lep4vecs"});
        _rlm = _rlm.Define("taujetoverlap", checkoverlap, {"jet4vecs","cleantau4vecs"})
                   .Define("loosetaujetoverlap", checkoverlap, {"jet4vecs","cleanloosetau4vecs"})
                   .Define("jetoverlap","lepjetoverlap && taujetoverlap")
                   .Define("jetoverlaploose","lepjetoverlap && loosetaujetoverlap");

    } else{
        _rlm = _rlm.Define("jetoverlap", checkoverlap, {"jet4vecs","lep4vecs"});
        _rlm = _rlm.Define("jetoverlaploose", checkoverlap, {"jet4vecs","lep4vecs"});
    }
    _rlm = _rlm.Define("Jet_pt_loose", "Jet_pt[jetoverlaploose]")
               .Define("Jet_btagPNetB_loose", "Jet_btagPNetB[jetoverlaploose]")
               .Define("ncleanjetsloosepass", "int(Jet_pt_loose.size())");
    if (_isRun24){
        _rlm = _rlm.Define("Jet_btagUParTAK4B_loose", "Jet_btagUParTAK4B[jetoverlaploose]");
    }

    _rlm = _rlm.Redefine("Jet_pt", "Jet_pt[jetoverlap]")
               .Redefine("Jet_eta", "Jet_eta[jetoverlap]")
               .Redefine("Jet_phi", "Jet_phi[jetoverlap]")
               .Redefine("Jet_mass", "Jet_mass[jetoverlap]")
               .Redefine("Jet_btagPNetB", "Jet_btagPNetB[jetoverlap]")
               .Define("ncleanjetspass", "int(Jet_pt.size())")
               .Define("cleanjet4vecs", ::gen4vec, {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass"})
               .Define("Jet_HT", "Sum(Jet_pt)");
    if (_isRun24){
        _rlm = _rlm.Redefine("Jet_btagUParTAK4B", "Jet_btagUParTAK4B[jetoverlap]");
    }

    if (!_isData) {
        int nbsf_var = btag_var.size();
        int njes_var = jes_var.size();
        _rlm = _rlm.Define("nbsf_var", [nbsf_var](){return int(nbsf_var);})
                   .Define("njes_var", [njes_var](){return int(njes_var);})
                   //.Define("btagWeight_PNetB_perJet_loose", skimCol, {"btagWeight_PNetB_perJet", "jetoverlaploose"})
                   //.Define("btagWeight_PNetB_loose", calcBSF, {"btagWeight_PNetB_perJet_loose", "nbsf_var"})
                   //.Redefine("btagWeight_PNetB_perJet", skimCol, {"btagWeight_PNetB_perJet", "jetoverlap"})
                   //.Redefine("btagWeight_PNetB_jes_perJet", skimCol, {"btagWeight_PNetB_jes_perJet", "jetoverlap"})
                   .Define("btagWeight", calcBSF, {"btagWeight", "nbsf_var"})
                   //.Define("btagWeight_PNetB_jes", calcBSF, {"btagWeight_PNetB_jes_perJet", "njes_var"});
    }


    // b-tagging
    // https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer22/#ak4-b-tagging
    if (_isRun22) { 
        _rlm = _rlm.Define("btagcuts", "Jet_btagPNetB>0.245") //l: 0.047, m: 0.245, t: 0.6734
                   .Define("btagcuts_loose", "Jet_btagPNetB>0.047");
    } else if (_isRun22EE) { 
        _rlm = _rlm.Define("btagcuts", "Jet_btagPNetB>0.2605") //l: 0.0499, m: 0.2605, t: 0.6915
                   .Define("btagcuts_loose", "Jet_btagPNetB>0.0499");
    } else if (_isRun23) { 
        _rlm = _rlm.Define("btagcuts", "Jet_btagPNetB>0.1917") //l: 0.0358, m: 0.1917, t: 0.6172
                   .Define("btagcuts_loose", "Jet_btagPNetB>0.0358");
    } else if (_isRun23BPix) { 
        _rlm = _rlm.Define("btagcuts", "Jet_btagPNetB>0.1919") //l: 0.0359, m: 0.1919, t: 0.6133
                   .Define("btagcuts_loose", "Jet_btagPNetB>0.0359");
    } else if (_isRun24){
        // TODO: copy from 23BPix
        _rlm = _rlm.Define("btagcuts", "Jet_btagUParTAK4B>0.1272") //l: 0.0246, m: 0.1272, t: 0.4648
                   .Define("btagcuts_loose", "Jet_btagUParTAK4B>0.0246");
    }

    _rlm = _rlm.Define("bJet_pt", "Jet_pt[btagcuts]")
               .Define("bJet_eta", "Jet_eta[btagcuts]")
               .Define("bJet_phi", "Jet_phi[btagcuts]")
               .Define("bJet_mass", "Jet_mass[btagcuts]")
               .Define("bJet_btagPNetB", "Jet_btagPNetB[btagcuts]")
               .Define("ncleanbjetspass", "int(bJet_pt.size())")
               .Define("bJet_HT", "Sum(bJet_pt)")
               .Define("cleanbjet4vecs", ::gen4vec, {"bJet_pt", "bJet_eta", "bJet_phi", "bJet_mass"})
               .Define("bJet_pt_loose", "Jet_pt_loose[btagcuts_loose]")
               .Define("ncleanbjetsloosepass", "int(bJet_pt_loose.size())");
    if (_isRun24){
        _rlm = _rlm.Define("bJet_btagUParTAK4B", "Jet_btagUParTAK4B[btagcuts]");
    }
}

void NanoAODAnalyzerrdframe::selectTaus(string cut, string tauYear) {
    auto syst_unc = _syst;

    //TES var.
    if (!_isData) {

        auto selectTES = [syst_unc](floatsVec unc)->floats {

            int idx = -1;
            if (syst_unc.find("tesup") != std::string::npos) idx = 0;
            else if (syst_unc.find("tesdown") != std::string::npos) idx = 1;
            floats selected;
            selected.reserve(unc.size());

            for (size_t i=0; i<unc.size(); i++) {
                if (idx < 0) selected.emplace_back(1.0f);
                selected.emplace_back(unc[i][idx]);
                std::cout << idx << " " << unc[i][idx]  << endl;
            }
            return selected;
        };

        if (_syst.find("tes") != std::string::npos) {
          _rlm = _rlm.Define("Tau_pt_unc_toapply", selectTES, {"Tau_pt_unc"})
                     .Redefine("Tau_pt", "Tau_pt * Tau_pt_unc_toapply")
                     .Redefine("Tau_mass", "Tau_mass * Tau_pt_unc_toapply");
        }
    }

    auto overlap_removal_leptau = [](FourVectorVec &lep4vecs, FourVectorVec &tau4vecs) {
        ints out;
        for (auto tau: tau4vecs) {
            int check = 0;
            for (auto lep: lep4vecs) {
                auto dR = ROOT::Math::VectorUtil::DeltaR(lep, tau);
                if( dR >= 0.4 ) check = 1;
            }
            out.emplace_back(check);
        }
        return out;
    };

    _rlm = _rlm.Define("tau4vecs", ::gen4vec, {"Tau_pt", "Tau_eta", "Tau_phi", "Tau_mass"});
    _rlm = _rlm.Define("leptauoverlap", overlap_removal_leptau, {"lep4vecs","tau4vecs"});

    // input vector: vec[pt][vars]
    auto skimCol = [this](floatsVec toSkim, ints cut)->floatsVec {

        floatsVec out;
        for (size_t i=0; i<toSkim.size(); i++) {
            if (cut[i] > 0) out.emplace_back(toSkim[i]);
        }
        if (out.size() > 0) return out;
        else return {{1.0f, 1.0f, 1.0f}};
    };

    // Fake factor study - loose but not tight
    // Hadronic Tau Object Selections
    _rlm = _rlm.Define("taucuts", cut)
               .Define("deeptauidcuts","Tau_idDeepTau2018v2p5VSmu >= 4 && Tau_idDeepTau2018v2p5VSe >= 6 && Tau_idDeepTau2018v2p5VSjet >= 7")
               .Define("deeptauidcuts_loose","Tau_idDeepTau2018v2p5VSmu >= 4 && Tau_idDeepTau2018v2p5VSe >= 2 && Tau_idDeepTau2018v2p5VSjet >= 4");

    // Hadronic Tau Selection
    _rlm = _rlm.Define("seltaucuts_loose","taucuts && deeptauidcuts_loose && leptauoverlap")
               .Define("Tau_pt_loose", "Tau_pt[seltaucuts_loose]")
               .Define("Tau_pt_loose_gen", "Tau_pt[seltaucuts_loose]")
               .Define("Tau_eta_loose", "Tau_eta[seltaucuts_loose]")
               .Define("Tau_phi_loose", "Tau_phi[seltaucuts_loose]")
               .Define("Tau_mass_loose", "Tau_mass[seltaucuts_loose]")
               .Define("Tau_charge_loose", "Tau_charge[seltaucuts_loose]")
               .Define("Tau_decayMode_loose", "Tau_decayMode[seltaucuts_loose]")
               .Define("nloosetaupass", "int(Tau_pt_loose.size())")
               .Define("cleanloosetau4vecs", ::gen4vec, {"Tau_pt_loose", "Tau_eta_loose", "Tau_phi_loose", "Tau_mass_loose"});


    _rlm = _rlm.Define("seltaucuts","taucuts && deeptauidcuts && leptauoverlap")
               .Define("Tau_pt_gen", "Tau_pt[seltaucuts]")
               .Redefine("Tau_pt", "Tau_pt[seltaucuts]")
               .Redefine("Tau_eta", "Tau_eta[seltaucuts]")
               .Redefine("Tau_phi", "Tau_phi[seltaucuts]")
               .Redefine("Tau_mass", "Tau_mass[seltaucuts]")
               .Redefine("Tau_charge", "Tau_charge[seltaucuts]")
               .Redefine("Tau_jetIdx", "Tau_jetIdx[seltaucuts]")
               .Redefine("Tau_decayMode", "Tau_decayMode[seltaucuts]")
               .Define("ncleantaupass", "int(Tau_pt.size())")
               .Define("cleantau4vecs", ::gen4vec, {"Tau_pt", "Tau_eta", "Tau_phi", "Tau_mass"});

    if (!_isData) {
        // Tau SF
        cout << "Loading Tau SF" << endl;
        std::string tauid_vsjet = "VTight";
        std::string tauid_vse   = "Tight";
        std::string tauid_vsmu  = "Tight";

        cout << "Tau ID WP vsJet : " << tauid_vsjet << endl;
        cout << "Tau ID WP vsMuon : " << tauid_vsmu << endl;
        cout << "Tau ID WP vsElectron : " << tauid_vse << endl;

        auto tauSFreader = correction::CorrectionSet::from_file("data/TauIDSFs/tau_DeepTau2018v2p5_" + tauYear + ".json.gz");
        auto _tauidSFjet = tauSFreader->at("DeepTau2018v2p5VSjet");
        auto _tauidSFele = tauSFreader->at("DeepTau2018v2p5VSe");
        auto _tauidSFmu  = tauSFreader->at("DeepTau2018v2p5VSmu");

        auto tauSFIdVsJet = [this, _tauidSFjet, tauid_vsjet, tauid_vse](floats &pt, uchars &dm, uchars &genmatch)->floatsVec {
            //std::cout << "tauSFIdVsJet vsjet_wp: "<< tauid_vsjet << std::endl;

            floats uncSources;
            uncSources.reserve(3);
            floatsVec wVec;
            wVec.reserve(pt.size());

            if (pt.size() > 0) {
                for (unsigned int i=0; i<pt.size(); i++) {
                    uncSources.emplace_back(_tauidSFjet->evaluate({pt[i], int(dm[i]), int(genmatch[i]), tauid_vsjet, tauid_vse, "nom", "dm"}));
                    uncSources.emplace_back(_tauidSFjet->evaluate({pt[i], int(dm[i]), int(genmatch[i]), tauid_vsjet, tauid_vse, "up", "dm"}));
                    uncSources.emplace_back(_tauidSFjet->evaluate({pt[i], int(dm[i]), int(genmatch[i]), tauid_vsjet, tauid_vse, "down", "dm"}));
                    wVec.emplace_back(uncSources);
                    uncSources.clear();
                }
            }
            return wVec;
        };


        auto tauSFIdVsEl = [this, _tauidSFele, tauid_vse](floats &eta, uchars &dm, uchars &genmatch)->floatsVec {

            //std::cout << "tauSFIdVsEl" << std::endl;
            floats uncSources;
            uncSources.reserve(3);
            floatsVec wVec;
            wVec.reserve(eta.size());

            if (eta.size() > 0) {
                for (unsigned int i=0; i<eta.size(); i++) {
                    uncSources.emplace_back(_tauidSFele->evaluate({eta[i], int(dm[i]), int(genmatch[i]), tauid_vse, "nom"}));
                    uncSources.emplace_back(_tauidSFele->evaluate({eta[i], int(dm[i]), int(genmatch[i]), tauid_vse, "up"}));
                    uncSources.emplace_back(_tauidSFele->evaluate({eta[i], int(dm[i]), int(genmatch[i]), tauid_vse, "down"}));
                    wVec.emplace_back(uncSources);
                    uncSources.clear();
                }
            }
            return wVec;
        };

        auto tauSFIdVsMu = [this, _tauidSFmu, tauid_vsmu](floats &eta, uchars &genmatch)->floatsVec {

            floats uncSources;
            uncSources.reserve(3);
            floatsVec wVec;
            wVec.reserve(eta.size());

            if (eta.size() > 0) {
                for (unsigned int i=0; i<eta.size(); i++) {
                    uncSources.emplace_back(_tauidSFmu->evaluate({eta[i], int(genmatch[i]), tauid_vsmu, "nom"}));
                    uncSources.emplace_back(_tauidSFmu->evaluate({eta[i], int(genmatch[i]), tauid_vsmu, "up"}));
                    uncSources.emplace_back(_tauidSFmu->evaluate({eta[i], int(genmatch[i]), tauid_vsmu, "down"}));
                    wVec.emplace_back(uncSources);
                    uncSources.clear();
                }
            }
            return wVec;
        };

        _rlm = _rlm.Define("tauWeightIdVsJet", tauSFIdVsJet, {"Tau_pt","Tau_decayMode","Tau_genPartFlav"})
                   .Define("tauWeightIdVsEl", tauSFIdVsEl, {"Tau_eta","Tau_decayMode", "Tau_genPartFlav"})
                   .Define("tauWeightIdVsMu", tauSFIdVsMu, {"Tau_eta","Tau_genPartFlav"});

        tauid_vsjet = "Loose";

        _rlm = _rlm.Define("tauWeightIdVsJet_loose", tauSFIdVsJet, {"Tau_pt","Tau_decayMode","Tau_genPartFlav"});
        _rlm = _rlm.Define("Tau_genPartFlav_loose","Tau_genPartFlav[seltaucuts_loose]")
                   .Define("taugencut_loose","Tau_genPartFlav_loose == 5")
                   .Redefine("Tau_pt_loose_gen", "Tau_pt_loose_gen[taugencut_loose]")
                   .Redefine("tauWeightIdVsJet_loose", skimCol, {"tauWeightIdVsJet_loose", "seltaucuts_loose"})
                   .Define("tauWeightIdVsEl_loose", skimCol, {"tauWeightIdVsEl", "seltaucuts_loose"})
                   .Define("tauWeightIdVsMu_loose", skimCol, {"tauWeightIdVsMu", "seltaucuts_loose"});

        _rlm = _rlm.Redefine("Tau_genPartFlav","Tau_genPartFlav[seltaucuts]")
                   .Define("taugencut","Tau_genPartFlav == 5")
                   .Redefine("Tau_pt_gen", "Tau_pt[taugencut]")
                   .Redefine("tauWeightIdVsJet", skimCol, {"tauWeightIdVsJet", "seltaucuts"})
                   .Redefine("tauWeightIdVsEl", skimCol, {"tauWeightIdVsEl", "seltaucuts"})
                   .Redefine("tauWeightIdVsMu", skimCol, {"tauWeightIdVsMu", "seltaucuts"});
    }
}

void NanoAODAnalyzerrdframe::matchGenReco() {
    if (_isMuonCh){
        _rlm = _rlm.Define("FinalGenPart_idx", ::FinalGenPart_idx, {"GenPart_pdgId", "GenPart_genPartIdxMother"});
    } else {
        _rlm = _rlm.Define("FinalGenPart_idx", ::FinalGenPart_idx_elec, {"GenPart_pdgId", "GenPart_genPartIdxMother"});
    }
    _rlm = _rlm.Define("GenPart_SMb_idx", "FinalGenPart_idx[3]")
               .Define("GenPart_SMW1_idx", "FinalGenPart_idx[4]")
               .Define("GenPart_SMW2_idx", "FinalGenPart_idx[5]")
               .Define("GenPart_SMtop_idx", "FinalGenPart_idx[7]")
               .Define("GenPart_d_SMb_idx", "FinalGenPart_idx[8]")
               .Define("GenPart_d_SMW_idx", "FinalGenPart_idx[9]")
               .Define("GenPart_d_SMW1_idx", "FinalGenPart_idx[10]")
               .Define("GenPart_d_SMW2_idx", "FinalGenPart_idx[11]")
               .Define("GenPart_d_SMtop_idx", "FinalGenPart_idx[12]");

    _rlm = _rlm.Define("drmax1", "float(0.15)")
               .Define("drmax2", "float(0.4)")
               .Define("GenJet_SMb_idx",::dRmatching,{"GenPart_SMb_idx","drmax2","GenPart_pt","GenPart_eta","GenPart_phi","GenPart_mass","GenJet_pt","GenJet_eta","GenJet_phi","GenJet_mass"})
               .Define("GenJet_SMW1_idx",::dRmatching,{"GenPart_SMW1_idx","drmax2","GenPart_pt","GenPart_eta","GenPart_phi","GenPart_mass","GenJet_pt","GenJet_eta","GenJet_phi","GenJet_mass"})
               .Define("GenJet_SMW2_idx",::dRmatching,{"GenPart_SMW2_idx","drmax2", "GenPart_pt","GenPart_eta","GenPart_phi","GenPart_mass","GenJet_pt","GenJet_eta","GenJet_phi","GenJet_mass"})
               .Define("Jet_SMb_idx",::dRmatching,{"GenJet_SMb_idx","drmax2","GenJet_pt","GenJet_eta","GenJet_phi","GenJet_mass","Jet_pt","Jet_eta","Jet_phi","Jet_mass"})
               .Define("Jet_SMW1_idx",::dRmatching,{"GenJet_SMW1_idx","drmax2","GenJet_pt","GenJet_eta","GenJet_phi","GenJet_mass","Jet_pt","Jet_eta","Jet_phi","Jet_mass"})
               .Define("Jet_SMW2_idx",::dRmatching,{"GenJet_SMW2_idx","drmax2", "GenJet_pt","GenJet_eta","GenJet_phi","GenJet_mass","Jet_pt","Jet_eta","Jet_phi","Jet_mass"});

    _rlm = _rlm.Define("GenPart_SMbmatched_4vecs", ::gen4vec_withidx, {"GenPart_pt","GenPart_eta","GenPart_phi","GenPart_mass", "GenPart_SMb_idx"})
               .Define("GenPart_SMW1matched_4vecs", ::gen4vec_withidx, {"GenPart_pt","GenPart_eta","GenPart_phi","GenPart_mass", "GenPart_SMW1_idx"})
               .Define("GenPart_SMW2matched_4vecs", ::gen4vec_withidx, {"GenPart_pt","GenPart_eta","GenPart_phi","GenPart_mass", "GenPart_SMW2_idx"});

    _rlm = _rlm.Define("GenPart_d_SMbmatched_4vecs", ::gen4vec_withidx, {"GenPart_pt","GenPart_eta","GenPart_phi","GenPart_mass", "GenPart_d_SMb_idx"})
               .Define("GenPart_d_SMWmatched_4vecs", ::gen4vec_withidx, {"GenPart_pt","GenPart_eta","GenPart_phi","GenPart_mass", "GenPart_d_SMW_idx"})
               .Define("GenPart_d_SMW1matched_4vecs", ::gen4vec_withidx, {"GenPart_pt","GenPart_eta","GenPart_phi","GenPart_mass", "GenPart_d_SMW1_idx"})
               .Define("GenPart_d_SMW2matched_4vecs", ::gen4vec_withidx, {"GenPart_pt","GenPart_eta","GenPart_phi","GenPart_mass", "GenPart_d_SMW2_idx"})
               .Define("GenPart_d_SMtopmatched_4vecs", ::gen4vec_withidx, {"GenPart_pt","GenPart_eta","GenPart_phi","GenPart_mass", "GenPart_d_SMtop_idx"})

               .Define("GenJet_SMbmatched_4vecs", ::gen4vec_withidx, {"GenJet_pt","GenJet_eta","GenJet_phi","GenJet_mass", "GenJet_SMb_idx"})
               .Define("GenJet_SMW1matched_4vecs", ::gen4vec_withidx, {"GenJet_pt","GenJet_eta","GenJet_phi","GenJet_mass", "GenJet_SMW1_idx"})
               .Define("GenJet_SMW2matched_4vecs", ::gen4vec_withidx, {"GenJet_pt","GenJet_eta","GenJet_phi","GenJet_mass", "GenJet_SMW2_idx"})

               .Define("Jet_SMbmatched_4vecs", ::gen4vec_withidx, {"Jet_pt","Jet_eta","Jet_phi","Jet_mass", "Jet_SMb_idx"})
               .Define("Jet_SMW1matched_4vecs", ::gen4vec_withidx, {"Jet_pt","Jet_eta","Jet_phi","Jet_mass", "Jet_SMW1_idx"})
               .Define("Jet_SMW2matched_4vecs", ::gen4vec_withidx, {"Jet_pt","Jet_eta","Jet_phi","Jet_mass", "Jet_SMW2_idx"});

    cout << "match Gen Reco ended" << endl;
}

void NanoAODAnalyzerrdframe::selectFatJets() {
    _rlm = _rlm.Define("fatjetcuts", "FatJet_pt>400.0 && abs(FatJet_eta)<2.4 && FatJet_tau1>0.0 && FatJet_tau2>0.0 && FatJet_tau3>0.0 && FatJet_tau3/FatJet_tau2<0.5")
               .Define("Sel_fatjetpt", "FatJet_pt[fatjetcuts]")
               .Define("Sel_fatjeteta", "FatJet_eta[fatjetcuts]")
               .Define("Sel_fatjetphi", "FatJet_phi[fatjetcuts]")
               .Define("Sel_fatjetmass", "FatJet_mass[fatjetcuts]")
               .Define("nfatjetspass", "int(Sel_fatjetpt.size())")
               //.Define("Sel_fatjetweight", "std::vector<double>(nfatjetspass, evWeight)")
               .Define("Sel_fatjet4vecs", ::gen4vec, {"Sel_fatjetpt", "Sel_fatjeteta", "Sel_fatjetphi", "Sel_fatjetmass"});
}

void NanoAODAnalyzerrdframe::topPtReweight() {
    _rlm = _rlm.Define("gentopcut", "abs(GenPart_pdgId) == 6 && GenPart_statusFlags & (1 << 13)")
               .Define("GenPart_top_pt", "GenPart_pt[gentopcut]");

    auto topPtLOtoNLO = [](floats toppt)->float {

        float out = 1.0;
        std::vector<float> xbins = {0,50,100,150,200,250,300,350,400,450,500,550,600,800,1000,2000};
        std::vector<float> sfs = {1.1536, 1.0578, 0.993, 0.9373, 0.8881, 0.8456, 0.8087, 0.7809,
                                  0.7559, 0.7388, 0.7247, 0.7245, 0.7124, 0.7284, 0.7317};

        if (toppt.size() != 2) out = 1.0;
        else {
            float pt1 = toppt[0];
            float pt2 = toppt[1];
            if (pt1 >= 2000) pt1 = 1999;
            if (pt2 >= 2000) pt2 = 1999;
            int xbin1 = (std::upper_bound(xbins.begin(), xbins.end(), pt1)-1) - xbins.begin();
            int xbin2 = (std::upper_bound(xbins.begin(), xbins.end(), pt2)-1) - xbins.begin();
            out = std::sqrt(sfs.at(xbin1) * sfs.at(xbin2));
        }
        return out;
    };

    auto topPtNLOtoNNLO = [](floats toppt)->floats {

        floats out;
        out.reserve(3);
        if (toppt.size() != 2) {
            out.emplace_back(1.0); out.emplace_back(1.0); out.emplace_back(1.0);
        } else {
            float pt1 = toppt[0];
            float pt2 = toppt[1];
            float nom_w1 = (0.103 * std::exp(-0.0118 * pt1) - 0.000134 * pt1 + 0.973);
            float nom_w2 = (0.103 * std::exp(-0.0118 * pt2) - 0.000134 * pt2 + 0.973);
            float nom_weight = std::sqrt(nom_w1 * nom_w2);
            out.emplace_back(nom_weight);

            // blessed by TOP-22-003 team (ANv16) this goes to nom weight fct
            std::vector<float> xbins = {0,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,
                                        1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500};
            std::vector<float> unc_up = {0, 0.2634291, 0.125821, 0.1204998, 0.1185625, 0.1225926, 0.1286782,
                                         0.1353967, 0.1415252, 0.1471104, 0.1555594, 0.1580664, 0.1604668,
                                         0.163553, 0.1690409, 0.1725839, 0.1755679, 0.1779287, 0.1808263,
                                         0.1841962, 0.1866275, 0.1877652, 0.1897599, 0.1920097, 0.20095 };
            std::vector<float> unc_dn = {0, 0.4420274, 0.1457186, 0.1366773, 0.1278782, 0.1296307, 0.1372049,
                                         0.1473315, 0.1561993, 0.1645779, 0.168799, 0.1779081, 0.1835472,
                                         0.1903936, 0.1969376, 0.2013376, 0.2056985, 0.2082445, 0.214616,
                                         0.2186682, 0.2240442, 0.2240283, 0.2287682, 0.2321543, 0.249874 };

            if (pt1 >= 2500) pt1 = 2499;
            if (pt2 >= 2500) pt2 = 2499;
            int xbin1 = (std::upper_bound(xbins.begin(), xbins.end(), pt1)-1) - xbins.begin();
            int xbin2 = (std::upper_bound(xbins.begin(), xbins.end(), pt2)-1) - xbins.begin();

            float up = nom_weight + std::sqrt(std::pow((nom_w2 * unc_up.at(xbin1)), 2) * std::pow((nom_w1 * unc_up.at(xbin2)), 2));
            float dn = nom_weight - std::sqrt(std::pow((nom_w2 * unc_dn.at(xbin1)), 2) * std::pow((nom_w1 * unc_dn.at(xbin2)), 2));
            out.emplace_back(up); out.emplace_back(dn);
        }

        return out;
    };

    // NLO ttbar: NLO to theory weight
    // https://twiki.cern.ch/twiki/bin/view/CMS/TopPtReweighting#TOP_PAG_corrections_based_on_the
    if (_outfilename.find("TTto") != std::string::npos) {
        _rlm = _rlm.Define("TopPtWeight", topPtNLOtoNNLO, {"GenPart_top_pt"});
    } else if (_outfilename.find("TT_LFV") != std::string::npos) {
        _rlm = _rlm.Define("TopPtWeight_LO", topPtLOtoNLO, {"GenPart_top_pt"})
                   .Define("TopPtWeight_NLO", topPtNLOtoNNLO, {"GenPart_top_pt"})
                   //.Define("TopPtWeight", "TopPtWeight_LO * TopPtWeight_NLO");
                   .Define("TopPtWeight", "floats v{TopPtWeight_LO * TopPtWeight_NLO[0],\
                                                    TopPtWeight_NLO[1] + (TopPtWeight_LO - 1) * TopPtWeight_NLO[0],\
                                                    TopPtWeight_NLO[2] + (TopPtWeight_LO - 1) * TopPtWeight_NLO[0]};\
                                                    return v;"); //keep magnitude of error but change nominal value
    } else {
        _rlm = _rlm.Define("TopPtWeight", "floats v{1.0, 1.0, 1.0}; return v;");
    }
}

void NanoAODAnalyzerrdframe::calculateTauES(string tauYear, string tauid_vsjet, string tauid_vsmu, string tauid_vse) {
    // Tau SF
    cout << "Loading Tau SF" << endl;

    cout << "Tau ID WP vsJet : " << tauid_vsjet << endl;
    cout << "Tau ID WP vsMuon : " << tauid_vsmu << endl;
    cout << "Tau ID WP vsElectron : " << tauid_vse << endl;


    auto tauSFreader = correction::CorrectionSet::from_file("data/TauIDSFs/tau_DeepTau2018v2p5_" + tauYear + ".json.gz");
    auto _testool  = tauSFreader->at("tau_energy_scale");
    auto _tauVsJet = tauSFreader->at("DeepTau2018v2p5VSjet");
    auto _tauVsEle = tauSFreader->at("DeepTau2018v2p5VSe");
    auto _tauVsMu  = tauSFreader->at("DeepTau2018v2p5VSmu");

    // Tau ID SF
    cout << "Applying Tau ID SF" << endl;
    auto tauIdSF = [this, _tauVsJet, _tauVsEle, _tauVsMu, tauid_vsjet, tauid_vse, tauid_vsmu](floats &pts, floats &etas, uchars &dms, uchars &genids)->floats{
        floats xout;
        xout.reserve(pts.size());
        for (unsigned int i=0; i<pts.size(); i++){
            float sf = 1.0;
            if (int(dms[i]) != 5 && int(dms[i]) != 6) {
                if (int(genids[i]) == 5) { //genuine tau
                    string flag = "dm";
                    if (pts[i] > 140.) flag = "pt";
                    sf = _tauVsJet->evaluate({pts[i], int(dms[i]), int(genids[i]), tauid_vsjet, tauid_vse, "nom", flag});
                } else if (int(genids[i]) == 1 || int(genids[i]) == 3){ //genuine electron
                    sf = _tauVsEle->evaluate({etas[i], int(dms[i]), int(genids[i]), tauid_vse, "nom"});
                } else if (int(genids[i]) == 2 || int(genids[i]) == 4) { //genuine muon
                    sf = _tauVsMu->evaluate({etas[i], int(genids[i]), tauid_vsmu, tauid_vse, tauid_vsjet, "nom"});
                }
            }
            xout.emplace_back(sf);
        }
        return xout;
    };

    // Tau ES
    cout<<"Applying TauES on Genuine taus"<<endl;
    auto tauES = [this, _testool, tauid_vsjet, tauid_vse](std::string var) {
        return [this, _testool, tauid_vsjet, tauid_vse, var](floats &pts, floats &etas, uchars &dms, uchars &genids)->floats {
            floats xout;
            xout.reserve(pts.size());
            for (unsigned int i=0; i<pts.size(); i++) {
                float es = 1.0;
                int dm_tmp = int (dms[i]);
                try{
                    if (dm_tmp!=5 && dm_tmp!=6 && int(genids[i])!=2 && int(genids[i])!=4 && int(genids[i])!=6){
                        if (dm_tmp == 2) dm_tmp = 1;
                        es = _testool->evaluate({pts[i], etas[i], dm_tmp, int(genids[i]), "DeepTau2018v2p5", tauid_vsjet, tauid_vse, var});
                    }
                } catch (const exception &e) {
                    std::cout << "e: " << e.what() <<std::endl;
                    std::cout << "pt : " << pts[i] << " eta " << etas[i] << " dm " << int(dms[i]) << " genid " << int(genids[i]) << std::endl;
                }
                xout.emplace_back(es);
            }
            return xout;
        };
    };

    _rlm = _rlm.Define("Tau_pt_uncor", "Tau_pt")
           .Define("Tau_ES_nom", tauES("nom"), {"Tau_pt_uncor", "Tau_eta", "Tau_decayMode", "Tau_genPartFlav"})
           .Define("Tau_ES_up", tauES("up"), {"Tau_pt_uncor", "Tau_eta", "Tau_decayMode", "Tau_genPartFlav"})
           .Define("Tau_ES_down", tauES("down"), {"Tau_pt_uncor", "Tau_eta", "Tau_decayMode", "Tau_genPartFlav"})
           .Redefine("Tau_pt", "Tau_pt_uncor*Tau_ES_nom")
           .Redefine("Tau_mass", "Tau_mass*Tau_ES_nom");
}


void NanoAODAnalyzerrdframe::helper_1DHistCreator(std::string hname, std::string title, const int nbins, const double xlow, const double xhi, std::string rdfvar, std::string evWeight, RNode *anode) {

	RDF1DHist histojets = anode->Histo1D({hname.c_str(), title.c_str(), nbins, xlow, xhi}, rdfvar, evWeight); // Fill with weight given by evWeight
	_th1dhistos[hname] = histojets;
}

void NanoAODAnalyzerrdframe::helper_2DHistCreator(std::string hname, std::string title, const int nbinsx, const double xlow, const double xhi, const int nbinsy, const double ylow, const double yhi, std::string rdfvarX, std::string rdfvarY, std::string evWeight, RNode *anode) {

	RDF2DHist histojets = anode->Histo2D({hname.c_str(), title.c_str(), nbinsx, xlow, xhi, nbinsy, ylow, yhi}, rdfvarX, rdfvarY, evWeight); // Fill with weight given by evWeight
	_th2dhistos[hname] = histojets;
}

// Automatically loop to create
void NanoAODAnalyzerrdframe::setupCuts_and_Hists() {

    cout << "setting up definitions, cuts, and histograms" <<endl;

    for ( auto &c : _varinfovector) {
        if (c.mincutstep.length()==0) _rlm = _rlm.Define(c.varname, c.vardefinition);
    }

    for (auto &x : _hist1dinfovector) {
        std::string hpost = "";

        if (x.mincutstep.length()==0) {
            helper_1DHistCreator(std::string(x.hmodel.fName.Data())+hpost+x.systname,  std::string(x.hmodel.fTitle.Data()), x.hmodel.fNbinsX, x.hmodel.fXLow, x.hmodel.fXUp, x.varname, x.weightname+x.systname, &_rlm);
        }
    }

    for (auto &x : _hist2dinfovector) {
        std::string hpost = "";

        if (x.mincutstep.length()==0) {
            helper_2DHistCreator(std::string(x.hmodel.fName.Data())+hpost+x.systname,  std::string(x.hmodel.fTitle.Data()), x.hmodel.fNbinsX, x.hmodel.fXLow, x.hmodel.fXUp, x.hmodel.fNbinsY, x.hmodel.fYLow, x.hmodel.fYUp, x.varnameX, x.varnameY, x.weightname+x.systname, &_rlm);
        }
    }

    _rnt.setRNode(&_rlm);

    for (auto acut : _cutinfovector) {
        //std::string cutname = "S" + to_string(acut.idx.length()-1);
        std::string cutname = "S" + to_string(acut.idx.length());
        std::string hpost = "_"+cutname;
        RNode *r = _rnt.getParent(acut.idx)->getRNode();
        auto rnext = new RNode(r->Define(cutname, acut.cutdefinition));
        *rnext = rnext->Filter(cutname);

        for ( auto &c : _varinfovector) {
            if (acut.idx.compare(c.mincutstep)==0) *rnext = rnext->Define(c.varname, c.vardefinition);
        }
        for (auto &x : _hist1dinfovector) {
            if (acut.idx.compare(0, x.mincutstep.length(), x.mincutstep)==0) {
                bool reachedMax = false;
                if (x.maxcutstep.length() > 0 and acut.idx.compare(0, x.maxcutstep.length(), x.maxcutstep)>=0) reachedMax = true;
                if (!reachedMax) {
                    helper_1DHistCreator(std::string(x.hmodel.fName.Data())+hpost+x.systname, std::string(x.hmodel.fTitle.Data()), x.hmodel.fNbinsX, x.hmodel.fXLow, x.hmodel.fXUp, x.varname, x.weightname+x.systname, rnext);
                }
            }
        }

        for (auto &x : _hist2dinfovector) {
            if (acut.idx.compare(0, x.mincutstep.length(), x.mincutstep)==0) {
                bool reachedMax = false;
                if (x.maxcutstep.length() > 0 and acut.idx.compare(0, x.maxcutstep.length(), x.maxcutstep)>=0) reachedMax = true;
                if (!reachedMax) {
                    helper_2DHistCreator(std::string(x.hmodel.fName.Data())+hpost+x.systname, std::string(x.hmodel.fTitle.Data()), x.hmodel.fNbinsX, x.hmodel.fXLow, x.hmodel.fXUp, x.hmodel.fNbinsY, x.hmodel.fYLow, x.hmodel.fYUp, x.varnameX, x.varnameY, x.weightname+x.systname, rnext);
                }
            }
        }
        _rnt.addDaughter(rnext, acut.idx);

        /*
        _rlm = _rlm.Define(cutname, acut.cutdefinition);
        _rlm = _rlm.Filter(cutname);

        for ( auto &c : _varinfovector) {
            if (acut.idx.compare(c.mincutstep)==0) _rlm = _rlm.Define(c.varname, c.vardefinition);
        }
        for (auto &x : _hist1dinfovector) {
            if (acut.idx.compare(0, x.mincutstep.length(), x.mincutstep)==0) {
                helper_1DHistCreator(std::string(x.hmodel.fName)+hpost,  std::string(x.hmodel.fTitle)+hpost, x.hmodel.fNbinsX, x.hmodel.fXLow, x.hmodel.fXUp, x.varname, x.weightname);
            }
        }
        _rnt.addDaughter(&_rlm, acut.idx);
        */
    }
}

void NanoAODAnalyzerrdframe::add1DHist(TH1DModel histdef, std::string variable, std::string weight, string syst, string mincutstep, string maxcutstep) {

	_hist1dinfovector.push_back({histdef, variable, weight, syst, mincutstep, maxcutstep});
}


void NanoAODAnalyzerrdframe::add2DHist(TH2DModel histdef, std::string variableX, std::string variableY, std::string weight, string syst, string mincutstep, string maxcutstep) {

	_hist2dinfovector.push_back({histdef, variableX, variableY, weight, syst, mincutstep, maxcutstep});
}

void NanoAODAnalyzerrdframe::drawHists(RNode t) {

	cout << "processing" <<endl;
	t.Count();
}

void NanoAODAnalyzerrdframe::addVar(varinfo v) {

	_varinfovector.push_back(v);
}

void NanoAODAnalyzerrdframe::addVartoStore(string varname) {

    // varname is assumed to be a regular expression.
    // e.g. if varname is "Muon_eta" then "Muon_eta" will be stored
    // if varname=="Muon_.*", then any branch name that starts with "Muon_" string will
    // be saved
    _varstostore.push_back(varname);
    /*
    std::regex b(varname);
    bool foundmatch = false;
    for (auto a: _rlm.GetColumnNames()) {
        if (std::regex_match(a, b)) {
            _varstostore.push_back(a);
            foundmatch = true;
        }
    }
    */

}

void NanoAODAnalyzerrdframe::setupTree() {
    cout << "setting up Tree" <<endl;

    vector<RNodeTree *> rntends;
    _rnt.getRNodeLeafs(rntends);
    for (auto arnt: rntends) {
        RNode *arnode = arnt->getRNode();
        string nodename = arnt->getIndex();
        vector<string> varforthistree;
        std::map<string, int> varused;

        for (auto varname: _varstostore) {
            bool foundmatch = false;
            std::regex b(varname);
            for (auto a: arnode->GetColumnNames()) {
                if (std::regex_match(a, b) && varused[a]==0) {
                    varforthistree.push_back(a);
                    varused[a]++;
                    foundmatch = true;
                }
            }
            if (!foundmatch) {
                cout << varname << " not found at "<< nodename << endl;
            }
        }
        _varstostorepertree[nodename]  = varforthistree;
    }
}

void NanoAODAnalyzerrdframe::addCuts(string cut, string idx) {

	_cutinfovector.push_back({cut, idx});
}


void NanoAODAnalyzerrdframe::run(bool saveAll, string outtreename) {

    /*
    if (saveAll) {
        _rlm.Snapshot(outtreename, _outfilename);
    }
    else {
        // use the following if you want to store only a few variables
        _rlm.Snapshot(outtreename, _outfilename, _varstostore);
    }
    */

    vector<RNodeTree *> rntends;
    _rnt.getRNodeLeafs(rntends);
    _rnt.Print();
    cout << rntends.size() << endl;
    // on master, regex_replace doesn't work somehow
    //std::regex rootextension("\\.root");

    for (auto arnt: rntends) {
        string nodename = arnt->getIndex();
        //string outname = std::regex_replace(_outfilename, rootextension, "_"+nodename+".root");
        string outname = _outfilename;

        // if producing many root files due to branched selection criteria,  each root file will get a different name
        if (rntends.size()>1) outname.replace(outname.find(".root"), 5, "_"+nodename+".root");
        _outrootfilenames.push_back(outname);
        RNode *arnode = arnt->getRNode();
        cout << arnt->getIndex();
        //cout << ROOT::RDF::SaveGraph(_rlm) << endl;

        if (saveAll) {
            std::cout << "save all branches" << std::endl;
            arnode->Snapshot(outtreename, outname);
        } else {
            // use the following if you want to store only a few variables
            //arnode->Snapshot(outtreename, outname, _varstostore);
            cout << " writing branches" << endl;
            for (auto bname: _varstostorepertree[nodename]) {
                cout << bname << ", ";
                //try{
                //    arnode->Snapshot(outtreename, outname, {bname});
                //} catch (const std::exception &e){
                //    std::cout << "exception: " << e.what() << std::endl;
                //    std::cout << "BAD: " << bname << std::endl;
                //}
            }
            arnode->Snapshot(outtreename, outname, _varstostorepertree[nodename]);
            cout<<endl<<"Snapshot is done"<<endl;
        }
        _outrootfile = new TFile(outname.c_str(),"UPDATE");
        for (auto &h : _th1dhistos) {
            if (h.second.GetPtr() != nullptr) {
                h.second.GetPtr()->Print();
                h.second.GetPtr()->Write();
                //std::cout<<h.second->GetName()<<std::endl;
                //h.second->Write();
                //std::cout<<"Histogram is written"<<std::endl;
            }
        }
        for (auto &h : _th2dhistos) {
            if (h.second.GetPtr() != nullptr) {
                h.second.GetPtr()->Print();
                h.second.GetPtr()->Write();
                //std::cout<<h.second->GetName()<<std::endl;
                //h.second->Write();
                //std::cout<<"Histogram is written"<<std::endl;
            }
        }

        if (_isSkim == true) {
            TH1F* hPDFWeights = new TH1F("LHEPdfWeightSum", "LHEPdfWeightSum", 103, 0, 103);
            for (size_t i=0; i<PDFWeights.size(); i++)
                hPDFWeights->SetBinContent(i+1, PDFWeights[i]);

            TH1F* hPSWeights = new TH1F("PSWeightSum", "PSWeightSum", 4, 0, 4);
            for (size_t i=0; i<PSWeights.size(); i++)
                hPSWeights->SetBinContent(i+1, PSWeights[i]);

            TH1F* hScaleWeights = new TH1F("ScaleWeightSum", "ScaleWeightSum", 9, 0, 9);
            for (size_t i=0; i<ScaleWeights.size(); i++)
                hScaleWeights->SetBinContent(i+1, ScaleWeights[i]);
        }

        _outrootfile->Write(0, TObject::kOverwrite);
        _outrootfile->Close();
    }
}

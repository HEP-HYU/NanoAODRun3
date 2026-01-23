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

    // Skim switch
    if (_isSkim == true) {
        cout << "<< Start Skim NanoAOD >>" << endl;
    } else {
        cout << "<< Start Process NanoAOD >>" << endl;
    }

    // Channel(electron/muon) switch  // Default is electron channel
    if (_ch.find("muon") != std::string::npos) {
        cout << _ch << " channel" << endl;
        _isMuonCh = true;
    } else {
        cout << _ch << " channel" << endl;
        _isMuonCh = false;
    }

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
        cout << "Input file is data: " << _outfilename << endl;
    } else {
        _isData = false;
        cout << "Input file is MC: " << _outfilename << endl;
    }

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

    // Event weight for data it's always one. For MC, it depends on the sign
    if(_isSkim){
        _rlm = _rlm.Define("isData", "true");
        _rlm = _rlm.Define("lhereweight","one");
        _rlm = _rlm.Define("unitGenWeight", "one");

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
            
            std::string pileFile = "";
            std::string map = "";
            if (_isRun22) {
                pileFile = "2022_Summer22";
                map = "Collisions2022_355100_357900_eraBCD_GoldenJson";
            } else if (_isRun22EE) {
                pileFile = "2022_Summer22EE";
                map = "Collisions2022_359022_362760_eraEFG_GoldenJson";
            } else if (_isRun23) {
                pileFile = "2023_Summer23";
                map = "Collisions2023_366403_369802_eraBC_GoldenJson";
            } else if (_isRun23BPix) {
                pileFile = "2023_Summer23BPix";
                map = "Collisions2023_369803_370790_eraD_GoldenJson";
            } else if (_isRun24) {
                //TODO
                pileFile = "2023_Summer23BPix";
                map = "Collisions2023_369803_370790_eraD_GoldenJson";
            }
            auto puWeightreader = correction::CorrectionSet::from_file("data/LUM/"+pileFile+"/puWeights.json.gz");
            auto _puweight = puWeightreader->at(map);

            auto PuWeight = [this, _puweight](float x) -> floats {
                floats out;
                out.emplace_back(_puweight->evaluate({x, "nominal"}));
                out.emplace_back(_puweight->evaluate({x, "up"}));
                out.emplace_back(_puweight->evaluate({x, "down"}));

                return out;
            };

            if (_isRun24){
                _rlm = _rlm.Define("puWeight", "std::vector<double>{1.0, 1.0, 1.0}");
            }
            else {
                _rlm = _rlm.Define("puWeight", PuWeight, {"Pileup_nTrueInt"});
            }
        }    
    }

    std::vector<std::string> jes_var;

    // Object selection will be defined in sequence.
    // Selected objects will be stored in new vectors.
    if (_isSkim) {
        JetVetoMap();
        if (_isMuonCh){
            selectMuons();
        } else {
            selectElectrons();
        }
        setupJetMETCorrection(_globaltag, jes_var, jes_var_flav, "AK4PFPuppi", _isData);
        skimJets();
        if (!_isData){
            calculateEvWeight();
        //    applyBSFs(jes_var);
        }
    }
    else {
        calculateSF();
        if (_isMuonCh){
            selectElectrons();
        } else {
            selectMuons();
        }
        selectTaus();
        selectJets(jes_var, jes_var_flav);
        if (!_isData){
            topPtReweight();
        }
    }
    if (_isRun24){
        _rlm = _rlm.Define("puWeight", "std::vector<double>{1.0, 1.0, 1.0}");
    }
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

void NanoAODAnalyzerrdframe::JetVetoMap() {
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
        auto getJetId = [this, _jetid](floats &etas, floats &chHEF, floats &neHEF, floats &chEmEF, floats &neEmEF, floats &muEF, uchars &Nch, uchars &Nne)->floats{
            floats jetIds;
            jetIds.reserve(etas.size());
            for (int i=0; i<etas.size(); i++){
                float jetId = _jetid->evaluate({etas[i], chHEF[i], neHEF[i], chEmEF[i], neEmEF[i], muEF[i], float(Nch[i]), float(Nne[i]), float(Nch[i]+Nne[i])});
                jetIds.emplace_back(jetId);
            }
            return jetIds;
        };
        _rlm = _rlm.Define("Jet_passJetIdTightLepVeto", getJetId, {"Jet_eta", "Jet_chHEF", "Jet_neHEF", "Jet_chEmEF", "Jet_neEmEF", "Jet_muEF", "Jet_chMultiplicity", "Jet_neMultiplicity"});
    } 
    //TODO: below is for NanoAODv12
    //If the others are going to v15, this should be removed
    else{
        auto isJetIdTightLepVeto = [this](floats &etas, uchars &jetIds, floats &muEF, floats &chEmEF, floats &neHEF, floats neEmEF)->floats{
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
        _rlm = _rlm.Define("Jet_passJetIdTightLepVeto", isJetIdTightLepVeto, {"Jet_eta", "Jet_jetId", "Jet_muEF", "Jet_chEmEF", "Jet_neHEF", "Jet_neEmEF"});
    }

    
    _rlm = _rlm.Define("loosejetcuts", "Jet_pt>15 && Jet_passJetIdTightLepVeto && (Jet_neEmEF+Jet_chEmEF)<0.9");
    
    _rlm = _rlm.Define("Jet_pt_loosejet", "Jet_pt[loosejetcuts]")
               .Define("Jet_phi_loosejet", "Jet_phi[loosejetcuts]")
               .Define("Jet_eta_loosejet", "Jet_eta[loosejetcuts]")
               .Define("Jet_mass_loosejet", "Jet_mass[loosejetcuts]")
               .Define("njetloose", "int(Jet_pt_loosejet.size())")
               .Define("loosejet4vecs", ::gen4vec, {"Jet_pt_loosejet", "Jet_eta_loosejet", "Jet_phi_loosejet", "Jet_mass_loosejet"});
    
    _rlm = _rlm.Define("vetomuoncut","Muon_pt>15.0 && abs(Muon_eta)<2.4 && Muon_looseId && Muon_pfRelIso04_all<0.25");

    
    _rlm = _rlm.Define("Muon_pt_veto", "Muon_pt[vetomuoncut]")
               .Define("Muon_phi_veto", "Muon_phi[vetomuoncut]")
               .Define("Muon_eta_veto", "Muon_eta[vetomuoncut]")
               .Define("Muon_mass_veto", "Muon_mass[vetomuoncut]")
               .Define("nmuonveto", "int(Muon_pt_veto.size())")
               .Define("vetomuon4vecs", ::gen4vec, {"Muon_pt_veto", "Muon_eta_veto", "Muon_phi_veto", "Muon_mass_veto"});
    
    _rlm = _rlm.Define("vetomuonjetoverlap", checkoverlap, {"loosejet4vecs","vetomuon4vecs"});
    
    _rlm = _rlm.Redefine("Jet_pt_loosejet", "Jet_pt_loosejet[vetomuonjetoverlap]")
               .Redefine("Jet_phi_loosejet", "Jet_phi_loosejet[vetomuonjetoverlap]")
               .Redefine("Jet_eta_loosejet", "Jet_eta_loosejet[vetomuonjetoverlap]")
               .Redefine("Jet_mass_loosejet", "Jet_mass_loosejet[vetomuonjetoverlap]")
               .Redefine("njetloose", "int(Jet_pt_loosejet.size())")
               .Redefine("loosejet4vecs", ::gen4vec, {"Jet_pt_loosejet", "Jet_eta_loosejet", "Jet_phi_loosejet", "Jet_mass_loosejet"});
    
    _rlm = _rlm.Redefine("Jet_phi_loosejet", [](floats phis) {
                    floats corrected(phis.size());
                    std::transform(phis.begin(), phis.end(), corrected.begin(), [](float phi) {
                        return (std::abs(phi) > 3.141592653589790f) ? ((phi > 0 ? 1.0f : -1.0f) * 3.141592653589790f) : phi;
                    });
                    return corrected;
               }, {"Jet_phi_loosejet"})
               .Redefine("Jet_eta_loosejet", [](floats etas) {
                    floats corrected(etas.size());
                    std::transform(etas.begin(), etas.end(), corrected.begin(), [](float eta) {
                        return (std::abs(eta) > 5.191f) ? ((eta > 0 ? 1.0f : -1.0f) * 5.190f) : eta;
                    });
	                return corrected;
               }, {"Jet_eta_loosejet"});

    //Check if one of this is in the veto map if so skip the event
    // .Define("badjets",int())

    std::string jetFile = "";
    std::string map = "";
    
    if (_isRun22) {
        jetFile = "2022_Summer22";
        map = "Summer22_23Sep2023_RunCD_V1";
    } else if (_isRun22EE) {
        jetFile = "2022_Summer22EE";
        map = "Summer22EE_23Sep2023_RunEFG_V1";
    } else if (_isRun23) {
        jetFile = "2023_Summer23";
        map = "Summer23Prompt23_RunC_V1";
    } else if (_isRun23BPix) {
        jetFile = "2023_Summer23BPix";
        map = "Summer23BPixPrompt23_RunD_V1";
    } else if (_isRun24) {
        //TODO
        jetFile = "2023_Summer23BPix";
        map = "Summer23BPixPrompt23_RunD_V1";
    }
    
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

    //if (_isRun24){
    //    _rlm = _rlm.Define("Jet_isVeto_fpix", vetomap_fpix, {"Jet_eta_loosejet", "Jet_phi_loosejet"})
    //               .Redefine("events_isVeto", "Sum(Jet_isVeto_loose)+Sum(Jet_isVeto_fpix)");
    //}
}

void NanoAODAnalyzerrdframe::selectElectrons() {
    if (_isMuonCh) {
        cout<<"selectElectrons muon channel vetoelecuts"<<endl;
        _rlm = _rlm.Define("vetoelecuts", "Electron_pt>15.0 && abs(Electron_eta)<2.5 && Electron_cutBased == 1")
                   .Define("nvetoelepass","Sum(vetoelecuts)");
    } else  {
        cout<<"selectElectrons electron channel vetoelecuts"<<endl;
        _rlm = _rlm.Define("elecuts", "Electron_pt>50 && abs(Electron_eta)<2.5 && Electron_mvaIso_WP90")
                   .Define("vetoelecuts", "!elecuts && Electron_pt>15.0 && abs(Electron_eta)<2.5 && Electron_cutBased == 1");

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

void NanoAODAnalyzerrdframe::selectMuons() {
    if (_isMuonCh) {
        cout<<"selectMuons muon channel vetomuoncuts"<<endl;
        _rlm = _rlm.Define("muoncuts", "Muon_pt>50.0 && abs(Muon_eta)<2.4 && Muon_tightId && Muon_pfRelIso04_all<0.15");
        _rlm = _rlm.Define("vetomuoncuts", "!muoncuts && Muon_pt>15.0 && abs(Muon_eta)<2.4 && Muon_looseId && Muon_pfRelIso04_all<0.25");

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
    } else {
        cout<<"selectMuons electron channel vetomuoncuts"<<endl;
        _rlm = _rlm.Define("vetomuoncuts", "Muon_pt>15.0 && abs(Muon_eta)<2.4 && Muon_looseId && Muon_pfRelIso04_all<0.25")
                   .Define("nvetomuons","Sum(vetomuoncuts)");
    }
}

void NanoAODAnalyzerrdframe::calculateSF(){
    if (_isMuonCh) {
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
        cout<<"Loading Muon SF L527"<<endl;

        auto muonSFreader = correction::CorrectionSet::from_file("data/MuonSF/ScaleFactors_Muon_Z_ID_ISO_"+muonFile+"_schemaV2.json");
        cout<<"Loading Muon SF L528"<<endl;
        auto _muonid = muonSFreader->at("NUM_TightID_DEN_TrackerMuons");
        cout<<"Loading Muon SF L529"<<endl;
        auto _muoniso = muonSFreader->at("NUM_TightPFIso_DEN_TightID");
        cout<<"Loading Muon SF L530"<<endl;
        if (_isRun24) muonFile = "2023_BPix";
        auto muonHltSFreader = correction::CorrectionSet::from_file("data/MuonSF/ScaleFactors_Muon_Z_HLT_"+muonFile+"_eta_pt_schemaV2.json");
        cout<<"Loading Muon SF L531"<<endl;
        auto _muontrg = muonHltSFreader->at("NUM_IsoMu24_or_Mu50_or_CascadeMu100_or_HighPtTkMu100_DEN_CutBasedIdTight_and_PFIsoTight");

        cout<<"Loading Muon SF L534"<<endl;
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

        cout<<"Loading Muon SF L587"<<endl;
        _rlm = _rlm.Define("muonWeightId", muonSFId, {"Muon_pt","Muon_eta"})
                   .Define("muonWeightIso", muonSFIso, {"Muon_pt","Muon_eta"});
        if (_isRun24){
            _rlm = _rlm.Define("muonWeightTrg", "std::vector<double>{1.0, 1.0, 1.0}");
        }
        else{
            _rlm = _rlm.Define("muonWeightTrg", muonSFTrg, {"Muon_pt","Muon_eta"});
        }
        
        std::string muonYear = "";

        muonYear = _year + "_UL";
        cout<<"Loading Muon SF L595"<<endl;

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

        cout<<"Loading Muon SF L655"<<endl;

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
     
    } else {
        //Electron SF
        cout << "Loading Electron SF" << endl;
        std::string elecFile = _year;
        std::string elecYear = "";
        if (_isRun22) {
            elecFile = "2022_Summer22";
            elecYear = "2022Re-recoBCD";
        } else if (_isRun22EE) {
            elecFile = "2022_Summer22EE";
            elecYear = "2022Re-recoE+PromptFG";
        } else if (_isRun23) {
            elecFile = "2023_Summer23";
            elecYear = "2023PromptC";
        } else if (_isRun23BPix) {
            elecFile = "2023_Summer23BPix";
            elecYear = "2023PromptD";
        } else if (_isRun24){
            // TODO
            elecFile = "2024_Summer24";
            elecYear = "2024Prompt";
        }

        auto elecSFreader = correction::CorrectionSet::from_file("data/ElectronSF/"+elecFile+"/electron.json.gz");
        auto _elecreco = elecSFreader->at("Electron-ID-SF");
        auto _elecid = elecSFreader->at("Electron-ID-SF");
        if (_isRun24){
            auto elecIdSFreader = correction::CorrectionSet::from_file("data/ElectronSF/"+elecFile+"/electronID.json.gz");
            _elecid = elecIdSFreader->at("Electron-ID-SF");
        }

        auto elecSFreco = [this, _elecreco, elecYear](floats &pt, floats &eta, floats &phi)->floats {
            floats wVec;
            wVec.reserve(3); //cent, up, down
            if (pt.size() == 1){
                float sf=1.0; float sf_up=1.0; float sf_down=1.0;
                if (elecYear.find("2022")!=std::string::npos || elecYear.find("2024")!=std::string::npos){
                    if (pt[0] >= 20 && pt[0] < 75){
                        sf = _elecreco->evaluate({elecYear, "sf", "Reco20to75", eta[0], pt[0]});
                        sf_up = _elecreco->evaluate({elecYear, "sfup", "Reco20to75", eta[0], pt[0]});
                        sf_down = _elecreco->evaluate({elecYear, "sfdown", "Reco20to75", eta[0], pt[0]});
                    } else if (pt[0] >= 75){
                        sf = _elecreco->evaluate({elecYear, "sf", "RecoAbove75", eta[0], pt[0]});
                        sf_up = _elecreco->evaluate({elecYear, "sfup", "RecoAbove75", eta[0], pt[0]});
                        sf_down = _elecreco->evaluate({elecYear, "sfdown", "RecoAbove75", eta[0], pt[0]});
                    }
                } else if (elecYear.find("2023")!=std::string::npos){
                    if (pt[0] >= 20 && pt[0] < 75){
                        sf = _elecreco->evaluate({elecYear, "sf", "Reco20to75", eta[0], pt[0], phi[0]});
                        sf_up = _elecreco->evaluate({elecYear, "sfup", "Reco20to75", eta[0], pt[0], phi[0]});
                        sf_down = _elecreco->evaluate({elecYear, "sfdown", "Reco20to75", eta[0], pt[0], phi[0]});
                    } else if (pt[0] >= 75){
                        sf = _elecreco->evaluate({elecYear, "sf", "RecoAbove75", eta[0], pt[0], phi[0]});
                        sf_up = _elecreco->evaluate({elecYear, "sfup", "RecoAbove75", eta[0], pt[0], phi[0]});
                        sf_down = _elecreco->evaluate({elecYear, "sfdown", "RecoAbove75", eta[0], pt[0], phi[0]});
                    }
                }
                wVec.emplace_back(sf);
                wVec.emplace_back(sf_up);
                wVec.emplace_back(sf_down);
            } else wVec = {1.0, 1.0, 1.0};
            return wVec;
        };

        auto elecSFId = [this, _elecid, elecYear](floats &pt, floats &eta, floats &phi)->floats {
            floats wVec;
            wVec.reserve(3); //cent, up, down
            if (pt.size() == 1){
                float sf=1.0; float sf_up=1.0; float sf_down=1.0;
                if (elecYear.find("2022")!=std::string::npos || elecYear.find("2024")!=std::string::npos){
                    sf = _elecid->evaluate({elecYear, "sf", "wp90iso", eta[0], pt[0]});
                    sf_up = _elecid->evaluate({elecYear, "sfup", "wp90iso", eta[0], pt[0]});
                    sf_down = _elecid->evaluate({elecYear, "sfdown", "wp90iso", eta[0], pt[0]});
                } else if (elecYear.find("2023")!=std::string::npos){
                    sf = _elecid->evaluate({elecYear, "sf", "wp90iso", eta[0], pt[0], phi[0]});
                    sf_up = _elecid->evaluate({elecYear, "sfup", "wp90iso", eta[0], pt[0], phi[0]});
                    sf_down = _elecid->evaluate({elecYear, "sfdown", "wp90iso", eta[0], pt[0], phi[0]});
                }
                wVec.emplace_back(sf);
                wVec.emplace_back(sf_up);
                wVec.emplace_back(sf_down);
            } else wVec = {1.0, 1.0, 1.0};
            return wVec;
        };


        _rlm = _rlm.Define("elecWeightReco", elecSFId, {"Electron_pt","Electron_eta","Electron_phi"})
                   .Define("elecWeightId", elecSFId, {"Electron_pt","Electron_eta","Electron_phi"});
        if (_isRun24) _rlm = _rlm.Define("elecWeightTrg", "one");
        else {
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
    }
}

/*
void NanoAODAnalyzerrdframe::selectMET()
{
    _rlm = _rlm.Define("met4vec", ::genmet4vec, {"MET_pt","MET_phi"});
}*/

void NanoAODAnalyzerrdframe::setupJetMETCorrection(string globaltag, std::vector<std::string> jes_var, std::vector<std::string> jes_var_flav, std::string jetalgo, bool dataMc) {
    std::vector<JetCorrectionUncertainty*> regroupedUnc;
    std::vector<JetCorrectionUncertainty*> flavPureUnc;
    FactorizedJetCorrector* _jetCorrector;

    if (_globaltag != "") {
        cout << "Applying new JetMET corrections. GT: " + _globaltag + " on jetAlgo: AK4PFPuppi" << endl;
        string basedirectory = "data/jes/";

        string datamcflag = "";
        if (dataMc) datamcflag = "DATA";
        else datamcflag = "MC";

        // set file names that contain the parameters for corrections
        string dbfilenamel1 = basedirectory + _globaltag + "_" + datamcflag + "_L1FastJet_" + jetalgo + ".txt";
        string dbfilenamel2 = basedirectory + _globaltag + "_" + datamcflag + "_L2Relative_" + jetalgo + ".txt";
        string dbfilenamel3 = basedirectory + _globaltag + "_" + datamcflag + "_L3Absolute_" + jetalgo + ".txt";
        string dbfilenamel2l3 = basedirectory + _globaltag + "_" + datamcflag + "_L2L3Residual_" + jetalgo + ".txt";


        JetCorrectorParameters *L1JetCorrPar = new JetCorrectorParameters(dbfilenamel1);
        if (!L1JetCorrPar->isValid()) {
          std::cerr << "L1FastJet correction parameters not read" << std::endl;
          exit(1);
        }

        JetCorrectorParameters *L2JetCorrPar = new JetCorrectorParameters(dbfilenamel2);
        if (!L2JetCorrPar->isValid()) {
            std::cerr << "L2Relative correction parameters not read" << std::endl;
            exit(1);
        }

        JetCorrectorParameters *L3JetCorrPar = new JetCorrectorParameters(dbfilenamel3);
        if (!L3JetCorrPar->isValid()) {
            std::cerr << "L3Absolute correction parameters not read" << std::endl;
            exit(1);
        }

        JetCorrectorParameters *L2L3JetCorrPar = new JetCorrectorParameters(dbfilenamel2l3);
        if (!L2L3JetCorrPar->isValid()) {
            std::cerr << "L2L3Residual correction parameters not read" << std::endl;
            exit(1);
        }

        // to apply all the corrections, first collect them into a vector
        std::vector<JetCorrectorParameters> jetc;
        jetc.push_back(*L1JetCorrPar);
        jetc.push_back(*L2JetCorrPar);
        jetc.push_back(*L3JetCorrPar);
        jetc.push_back(*L2L3JetCorrPar);

        // apply the various corrections
        _jetCorrector = new FactorizedJetCorrector(jetc);

        // object to calculate uncertainty
        if (!dataMc) {
            cout<<"Applying JEC Uncertainty"<<endl;
            for (std::string src : jes_var) {
                if (src.find("up") != std::string::npos) {
                    if (src.find("HEM") != std::string::npos) continue;
                    auto uncsource = src.substr(3, src.size()-2-3);
                    cout << "JEC Uncertainty Source : " + uncsource << endl;
                    string dbfilenameunc = basedirectory + "RegroupedV2_" + _globaltag + "_MC_UncertaintySources_AK4PFPuppi.txt";
                    JetCorrectorParameters* uncCorrPar = new JetCorrectorParameters(dbfilenameunc, uncsource);
                    JetCorrectionUncertainty* _jetCorrectionUncertainty = new JetCorrectionUncertainty(*uncCorrPar);
                    regroupedUnc.emplace_back(_jetCorrectionUncertainty);
                } else {
                    continue; //We only need var name, no up/down
                }
            }
            //cout << "Treating breakdown of FlavorPure*" << endl;
            for (std::string src : jes_var_flav) {
                if (src.find("up") != std::string::npos) {
                    auto uncsource = src.substr(3, src.size()-2-3);
                    cout << "JEC Uncertainty Source : " + uncsource << endl;
                    string dbfilenameunc = basedirectory + _globaltag + "_MC_UncertaintySources_AK4PFPuppi.txt";
                    JetCorrectorParameters* uncCorrPar = new JetCorrectorParameters(dbfilenameunc, uncsource);
                    JetCorrectionUncertainty* _jetCorrectionUncertainty = new JetCorrectionUncertainty(*uncCorrPar);
                    flavPureUnc.emplace_back(_jetCorrectionUncertainty);
                } else {
                    continue; //We only need var name, no up/down
                }
            }
        }
    }

    auto applyJes = [this, _jetCorrector](floats jetpts, floats jetetas, floats jetphis, floats jetAreas, floats jetrawf, float rho, floats tocorrect)->floats {

        floats corrfactors;
        corrfactors.reserve(jetpts.size());

        for (unsigned int i=0; i<jetpts.size(); i++) {
            float rawfrac = 1.0-jetrawf[i];
            float rawjetpt = jetpts[i] * (rawfrac);
            _jetCorrector->setJetPt(rawjetpt);
            _jetCorrector->setJetEta(jetetas[i]);
            _jetCorrector->setJetPhi(jetphis[i]);
            _jetCorrector->setJetA(jetAreas[i]);
            _jetCorrector->setRho(rho);
            float corrfactor = _jetCorrector->getCorrection();
            if (abs(corrfactor) > 100.) corrfactor = 1.0;
            corrfactors.emplace_back(tocorrect[i] * rawfrac * corrfactor);
        }
        return corrfactors;
    };

    // structure: jes[jetIdx][varIdx]
    auto jesUnc = [this, jes_var, jes_var_flav, regroupedUnc, flavPureUnc](floats jetpts, floats jetetas, floats jetphis, floats jetAreas, floats jetrawf, float rho, shorts &partflav)->floatsVec {

        floats uncSources;
        uncSources.reserve(2 * regroupedUnc.size() + flavPureUnc.size());
        floatsVec uncertainties;
        uncertainties.reserve(jetpts.size());

        for (unsigned int i=0; i<jetpts.size(); i++) {
            for (size_t j=0; j<regroupedUnc.size(); j++) {
                auto corrector = regroupedUnc[j];

                corrector->setJetPt(jetpts[i]);
                corrector->setJetEta(jetetas[i]);
                float unc = corrector->getUncertainty(true);
                if (abs(unc) > 100.) unc = 0.;
                uncSources.emplace_back(1.0f + unc);
                uncSources.emplace_back(1.0f - unc);
            }
            // FlavorPure *
            for (size_t j=0; j<flavPureUnc.size(); j++) {
                auto corrector = flavPureUnc[j];

                bool isCorrectFlav = false;
                // The order of jes_var_flav is important: Gluon / Quark / Charm / Bottom
                if      (j == 0 and (abs(partflav[i]) == 21 or abs(partflav[i]) == 0)) isCorrectFlav = true;
                else if (j == 1 and (abs(partflav[i]) == 1 or abs(partflav[i]) == 2 or abs(partflav[i]) == 3)) isCorrectFlav = true;
                else if (j == 2 and abs(partflav[i]) == 4) isCorrectFlav = true;
                else if (j == 3 and abs(partflav[i]) == 5) isCorrectFlav = true;

                if (isCorrectFlav) {
                    corrector->setJetPt(jetpts[i]);
                    corrector->setJetEta(jetetas[i]);
                    float unc = corrector->getUncertainty(true);
                    if (abs(unc) > 100.) unc = 0.;
                    uncSources.emplace_back(1.0f + unc);
                    uncSources.emplace_back(1.0f - unc);
                } else {
                    uncSources.emplace_back(1.0f);
                    uncSources.emplace_back(1.0f);
                }
            }
            uncertainties.emplace_back(uncSources);
            uncSources.clear();
        }
        return uncertainties;
    };
    /*
    auto metCorr = [=](float met, float metphi, floats jetptsbefore, floats jetptsafter, floats jetphis, int npv, unsigned int _runnb)->float {

        int runnb = int(_runnb);

        auto metx = met * cos(metphi);
        auto mety = met * sin(metphi);

        for (unsigned int i=0; i<jetphis.size(); i++) {
            if (jetptsafter[i] > 15.0) {
                metx -= (jetptsafter[i] - jetptsbefore[i])*cos(jetphis[i]);
                mety -= (jetptsafter[i] - jetptsbefore[i])*sin(jetphis[i]);
            }
        }

        //starting phi modulation correction on met
        // https://lathomas.web.cern.ch/lathomas/METStuff/XYCorrections/XYMETCorrection_withUL17andUL18andUL16.h

        auto uncormet = float(sqrt(metx*metx + mety*mety));
        auto uncormet_phi = float(atan2(mety, metx));

        if(npv>100) npv=100;
        auto METxcorr(0.),METycorr(0.);
        if (_isData) {
            //UL2018
            if(runnb >= 315252 && runnb <= 316995 ) {METxcorr = -(0.263733 * npv + -1.91115); METycorr = -(0.0431304 * npv + -0.112043);}
            if(runnb >= 316998 && runnb <= 319312 ) {METxcorr = -(0.400466 * npv + -3.05914); METycorr = -(0.146125 * npv + -0.533233);}
            if(runnb >= 319313 && runnb <= 320393 ) {METxcorr = -(0.430911 * npv + -1.42865); METycorr = -(0.0620083 * npv + -1.46021);}
            if(runnb >= 320394 && runnb <= 325273 ) {METxcorr = -(0.457327 * npv + -1.56856); METycorr = -(0.0684071 * npv + -0.928372);}
            //UL2017
            if(runnb >= 297020 && runnb <= 299329 ) {METxcorr = -(-0.211161 * npv + 0.419333); METycorr = -(0.251789 * npv + -1.28089);}
            if(runnb >= 299337 && runnb <= 302029 ) {METxcorr = -(-0.185184 * npv + -0.164009); METycorr = -(0.200941 * npv + -0.56853);}
            if(runnb >= 302030 && runnb <= 303434 ) {METxcorr = -(-0.201606 * npv + 0.426502); METycorr = -(0.188208 * npv + -0.58313);}
            if(runnb >= 303435 && runnb <= 304826 ) {METxcorr = -(-0.162472 * npv + 0.176329); METycorr = -(0.138076 * npv + -0.250239);}
            if(runnb >= 304911 && runnb <= 306462 ) {METxcorr = -(-0.210639 * npv + 0.72934); METycorr = -(0.198626 * npv + 1.028);}
            //UL2016
            if(runnb >= 272007 && runnb <= 275376 ) {METxcorr = -(-0.0214894 * npv + -0.188255); METycorr = -(0.0876624 * npv + 0.812885);}
            if(runnb >= 275657 && runnb <= 276283 ) {METxcorr = -(-0.032209 * npv + 0.067288); METycorr = -(0.113917 * npv + 0.743906);}
            if(runnb >= 276315 && runnb <= 276811 ) {METxcorr = -(-0.0293663 * npv + 0.21106); METycorr = -(0.11331 * npv + 0.815787);}
            if(runnb >= 276831 && runnb <= 277420 ) {METxcorr = -(-0.0132046 * npv + 0.20073); METycorr = -(0.134809 * npv + 0.679068);}
            if(((runnb >= 277772 && runnb <= 278768) || runnb==278770)) {METxcorr = -(-0.0543566 * npv + 0.816597); METycorr = -(0.114225 * npv + 1.17266);};
            if(((runnb >= 278801 && runnb <= 278808) || runnb==278769)) {METxcorr = -(0.134616 * npv + -0.89965); METycorr = -(0.0397736 * npv + 1.0385);};
            if(runnb >= 278820 && runnb <= 280385 ) {METxcorr = -(0.121809 * npv + -0.584893); METycorr = -(0.0558974 * npv + 0.891234);};
            if(runnb >= 280919 && runnb <= 284044 ) {METxcorr = -(0.0868828 * npv + -0.703489); METycorr = -(0.0888774 * npv + 0.902632);};
        }
        auto CorrectedMET_x = uncormet * cos( uncormet_phi ) + METxcorr;
        auto CorrectedMET_y = uncormet * sin( uncormet_phi ) + METycorr;

        float CorrectedMET = sqrt(CorrectedMET_x * CorrectedMET_x + CorrectedMET_y * CorrectedMET_y);

        return CorrectedMET;
    };

    auto metUnc = [](float met, float metphi, floats jetptsbefore, floatsVec jetptscorr, floats jetphis)->floats {

        floats corrfactors;
        corrfactors.reserve(static_cast<int>(jetptscorr.size()));

        for (size_t j=0; j<jetptscorr[0].size(); j++) {
            auto metx = met * cos(metphi);
            auto mety = met * sin(metphi);

            for (unsigned int i=0; i<jetphis.size(); i++) {
                if (jetptsbefore[i] * jetptscorr[i][j] > 15.0) {
                    metx -= (jetptscorr[i][j] - 1.0) * jetptsbefore[i] * cos(jetphis[i]);
                    mety -= (jetptscorr[i][j] - 1.0) * jetptsbefore[i] * sin(jetphis[i]);
                }
            }
            corrfactors.emplace_back(float(sqrt(metx*metx + mety*mety)));
        }
        return corrfactors;
    };

    auto metPhiCorr = [=](float met, float metphi, floats jetptsbefore, floats jetptsafter, floats jetphis, int npv, unsigned int _runnb)->float {

        int runnb = int(_runnb);

        auto metx = met * cos(metphi);
        auto mety = met * sin(metphi);

        for (unsigned int i=0; i<jetphis.size(); i++) {
            if (jetptsafter[i] > 15.0) {
                metx -= (jetptsafter[i] - jetptsbefore[i])*cos(jetphis[i]);
                mety -= (jetptsafter[i] - jetptsbefore[i])*sin(jetphis[i]);
            }
        }

        //starting phi modulation correction on met phi
        // https://lathomas.web.cern.ch/lathomas/METStuff/XYCorrections/XYMETCorrection_withUL17andUL18andUL16.h

        auto uncormet = float(sqrt(metx*metx + mety*mety));
        auto uncormet_phi = float(atan2(mety, metx));

        if(npv>100) npv=100;
        auto METxcorr(0.),METycorr(0.);
        if (_isData) {
            //UL2018
            if (runnb >= 315252 && runnb <= 316995 ) {METxcorr = -(0.263733 * npv + -1.91115); METycorr = -(0.0431304 * npv + -0.112043);}
            if (runnb >= 316998 && runnb <= 319312 ) {METxcorr = -(0.400466 * npv + -3.05914); METycorr = -(0.146125 * npv + -0.533233);}
            if (runnb >= 319313 && runnb <= 320393 ) {METxcorr = -(0.430911 * npv + -1.42865); METycorr = -(0.0620083 * npv + -1.46021);}
            if (runnb >= 320394 && runnb <= 325273 ) {METxcorr = -(0.457327 * npv + -1.56856); METycorr = -(0.0684071 * npv + -0.928372);}
            //UL2017
            if (runnb >= 297020 && runnb <= 299329 ) {METxcorr = -(-0.211161 * npv + 0.419333); METycorr = -(0.251789 * npv + -1.28089);}
            if (runnb >= 299337 && runnb <= 302029 ) {METxcorr = -(-0.185184 * npv + -0.164009); METycorr = -(0.200941 * npv + -0.56853);}
            if (runnb >= 302030 && runnb <= 303434 ) {METxcorr = -(-0.201606 * npv + 0.426502); METycorr = -(0.188208 * npv + -0.58313);}
            if (runnb >= 303435 && runnb <= 304826 ) {METxcorr = -(-0.162472 * npv + 0.176329); METycorr = -(0.138076 * npv + -0.250239);}
            if (runnb >= 304911 && runnb <= 306462 ) {METxcorr = -(-0.210639 * npv + 0.72934); METycorr = -(0.198626 * npv + 1.028);}
            //UL2016
            if (runnb >= 272007 && runnb <= 275376 ) {METxcorr = -(-0.0214894 * npv + -0.188255); METycorr = -(0.0876624 * npv + 0.812885);}
            if (runnb >= 275657 && runnb <= 276283 ) {METxcorr = -(-0.032209 * npv + 0.067288); METycorr = -(0.113917 * npv + 0.743906);}
            if (runnb >= 276315 && runnb <= 276811 ) {METxcorr = -(-0.0293663 * npv + 0.21106); METycorr = -(0.11331 * npv + 0.815787);}
            if (runnb >= 276831 && runnb <= 277420 ) {METxcorr = -(-0.0132046 * npv + 0.20073); METycorr = -(0.134809 * npv + 0.679068);}
            if (((runnb >= 277772 && runnb <= 278768) || runnb==278770)) {METxcorr = -(-0.0543566 * npv + 0.816597); METycorr = -(0.114225 * npv + 1.17266);};
            if (((runnb >= 278801 && runnb <= 278808) || runnb==278769)) {METxcorr = -(0.134616 * npv + -0.89965); METycorr = -(0.0397736 * npv + 1.0385);};
            if (runnb >= 278820 && runnb <= 280385 ) {METxcorr = -(0.121809 * npv + -0.584893); METycorr = -(0.0558974 * npv + 0.891234);};
            if (runnb >= 280919 && runnb <= 284044 ) {METxcorr = -(0.0868828 * npv + -0.703489); METycorr = -(0.0888774 * npv + 0.902632);};
        }

        auto CorrectedMET_x = uncormet * cos( uncormet_phi ) + METxcorr;
        auto CorrectedMET_y = uncormet * sin( uncormet_phi ) + METycorr;

        float CorrectedMETPhi;

        if      (CorrectedMET_x==0 && CorrectedMET_y>0) CorrectedMETPhi = TMath::Pi();
        else if (CorrectedMET_x==0 && CorrectedMET_y<0 ) CorrectedMETPhi = -TMath::Pi();
        else if (CorrectedMET_x >0) CorrectedMETPhi = TMath::ATan(CorrectedMET_y / CorrectedMET_x);
        else if (CorrectedMET_x <0 && CorrectedMET_y>0) CorrectedMETPhi = TMath::ATan(CorrectedMET_y / CorrectedMET_x) + TMath::Pi();
        else if (CorrectedMET_x <0 && CorrectedMET_y<0) CorrectedMETPhi = TMath::ATan(CorrectedMET_y / CorrectedMET_x) - TMath::Pi();
        else    CorrectedMETPhi =0;

        return CorrectedMETPhi;
    };

    auto metPhiUnc = [](float met, float metphi, floats jetptsbefore, floatsVec jetptscorr, floats jetphis)->floats {

        floats corrfactors;
        corrfactors.reserve(static_cast<int>(jetptscorr.size()));

        for (size_t j=0; j<jetptscorr[0].size(); j++) {
            auto metx = met * cos(metphi);
            auto mety = met * sin(metphi);

            for (unsigned int i=0; i<jetphis.size(); i++) {
                if (jetptsbefore[i] * jetptscorr[i][j] > 15.0) {
                    metx -= (jetptscorr[i][j] - 1.0) * jetptsbefore[i] * cos(jetphis[i]);
                    mety -= (jetptscorr[i][j] - 1.0) * jetptsbefore[i] * sin(jetphis[i]);
                }
            }
            corrfactors.emplace_back(float(atan2(mety, metx)));
        }
        return corrfactors;
    };
    */
    //FIXME: should correct jet mass. but can we do it at once?
    
    if (_jetCorrector != 0) {
        _rlm = _rlm.Define("Jet_pt_uncorr", "Jet_pt");
        _rlm = _rlm.Define("Jet_pt_corr", applyJes, {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_area", "Jet_rawFactor", "Rho_fixedGridRhoFastjetAll", "Jet_pt"})
	               .Redefine("Jet_mass", applyJes, {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_area", "Jet_rawFactor", "Rho_fixedGridRhoFastjetAll", "Jet_mass"});
	  //.Redefine("MET_pt", metCorr, {"MET_pt", "MET_phi", "Jet_pt", "Jet_pt_corr", "Jet_phi","PV_npvsGood", "run"})
	  //.Redefine("MET_phi", metPhiCorr, {"MET_pt", "MET_phi", "Jet_pt", "Jet_pt_corr", "Jet_phi", "PV_npvsGood", "run"});
        if (!_isData) {
	        _rlm = _rlm.Define("Jet_pt_unc", jesUnc, {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_area", "Jet_rawFactor", "Rho_fixedGridRhoFastjetAll", "Jet_partonFlavour"});
	                   //.Define("MET_pt_unc", metUnc, {"MET_pt", "MET_phi", "Jet_pt", "Jet_pt_unc", "Jet_phi"})
                       //.Define("MET_phi_unc", metPhiUnc, {"MET_pt", "MET_phi", "Jet_pt", "Jet_pt_unc", "Jet_phi"});
        }
        _rlm = _rlm.Redefine("Jet_pt", "Jet_pt_corr");
    }


    // JER
    std::string jetResFilePath_ = "data/jer/";
    std::string jetResSFFilePath_ = "data/jer/";

    if (_isRun22) {
        jetResFilePath_ += "Summer22_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi.txt";
        jetResSFFilePath_ += "Summer22_22Sep2023_JRV1_MC_SF_AK4PFPuppi.txt";
    } else if (_isRun22EE) {
        jetResFilePath_ += "Summer22EE_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi.txt";
        jetResSFFilePath_ += "Summer22EE_22Sep2023_JRV1_MC_SF_AK4PFPuppi.txt";
    } else if (_isRun23) {
        jetResFilePath_ += "Summer23Prompt23_RunCv1234_JRV1_MC_PtResolution_AK4PFPuppi.txt";
        jetResSFFilePath_ += "Summer23Prompt23_RunCv1234_JRV1_MC_SF_AK4PFPuppi.txt";
    } else if (_isRun23BPix) {
        jetResFilePath_ += "Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi.txt";
        jetResSFFilePath_ += "Summer23BPixPrompt23_RunD_JRV1_MC_SF_AK4PFPuppi.txt";
    } else if (_isRun24) {
        //TODO
        jetResFilePath_ += "Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi.txt";
        jetResSFFilePath_ += "Summer23BPixPrompt23_RunD_JRV1_MC_SF_AK4PFPuppi.txt";
    }

    
    JME::JetResolution jetResObj;
    JME::JetResolutionScaleFactor jetResSFObj;
    if (!_isData) {//JME::JetResolution
        jetResObj = JME::JetResolution(jetResFilePath_);
        jetResSFObj = JME::JetResolutionScaleFactor(jetResSFFilePath_);
    }

    // Compute the JER and Unc ( v[pt][unc], unc = nom, up, down)
    // cattool + PhysicsTools/PatUtils/interface/SmearedJetProducerT.h
    auto applyJer = [this, jetResObj, jetResSFObj](floats jetpts, floats jetetas, floats jetphis, floats jetms,
                    floats genjetpts, floats genjetetas, floats genjetphis, floats genjetms, shorts genidx, float rho, unsigned long long event)
                    ->floatsVec {

        floatsVec out;
        out.reserve(jetpts.size());
        floats var;
        var.reserve(3);

        if (jetpts.size() > 0) {
            for (size_t i=0; i<jetpts.size(); i++) {

                JME::JetParameters jetPars = {{JME::Binning::JetPt, jetpts[i]},
                                              {JME::Binning::JetEta, jetetas[i]},
                                              {JME::Binning::Rho, rho}};

                const float jetRes = static_cast<float>(jetResObj.getResolution(jetPars)); // Note: this is relative resolution.
                const float cJER   = static_cast<float>(jetResSFObj.getScaleFactor(jetPars));
                const float cJERUp = static_cast<float>(jetResSFObj.getScaleFactor(jetPars, Variation::UP));
                const float cJERDn = static_cast<float>(jetResSFObj.getScaleFactor(jetPars, Variation::DOWN));

                bool _isGenMatch = false;
                ROOT::Math::PtEtaPhiMVector jetv(jetpts[i], jetetas[i], jetphis[i], jetms[i]);
                ROOT::Math::PtEtaPhiMVector genv(genjetpts[i], genjetetas[i], genjetphis[i], genjetms[i]);
                if (genidx[i] >= 0 && ROOT::Math::VectorUtil::DeltaR(jetv, genv) < 0.2 &&
                    std::abs(genjetpts[i] - jetpts[i]) < jetRes * 3 * jetpts[i])
                    _isGenMatch = true;


                // JER - apply scaling method if matched genJet is found,
                //       apply gaussian smearing method if unmatched
                float jetpt = jetpts[i];
                float genjetpt = genjetpts[i];
                if (_isGenMatch) {
                    float dPt = jetpt - genjetpt;
                    float jersf = 1 + (dPt * (cJER - 1))/jetpt;
                    float jersfup = 1 + (dPt * (cJERUp - 1))/jetpt;
                    float jersfdn = 1 + (dPt * (cJERDn - 1))/jetpt;

                    floats tmpjers = {jersf, jersfup, jersfdn};
                    for (size_t j=0; j<tmpjers.size(); j++) {
                        float tmpjer = tmpjers[j];
                        if (std::isnan(tmpjer) or std::isinf(tmpjer) or tmpjer<0 ) var.emplace_back(1.0f);
                        else var.emplace_back(std::max(0.0f, tmpjer));
                    }

                } else if (cJER > 1){
                    std::uint32_t seed = uint32_t(jetetas[0] * 100) + static_cast<unsigned int>(event);
                    std::mt19937 m_random_generator = std::mt19937(seed);

                    float sigma = jetRes * std::sqrt(cJER * cJER - 1);
                    float sigmaUp = jetRes * std::sqrt(cJERUp * cJERUp - 1);
                    float sigmaDn = jetRes * std::sqrt(cJERDn * cJERDn - 1);
                    std::normal_distribution<float> d(0, sigma);
                    std::normal_distribution<float> dup(0, sigmaUp);
                    std::normal_distribution<float> ddn(0, sigmaDn);

                    floats tmpjers = {1.0f + d(m_random_generator), 1.0f + dup(m_random_generator), 1.0f + ddn(m_random_generator)};
                    for (size_t j=0; j<tmpjers.size(); j++) {
                        float tmpjer = tmpjers[j];
                        if (std::isnan(tmpjer) or std::isinf(tmpjer) or tmpjer<0 ) var.emplace_back(1.0f);
                        else var.emplace_back(std::max(0.0f, tmpjer));
                    }
                } else {
                    var = {1.0f, 1.0f, 1.0f};
                }
                out.emplace_back(var);
                var.clear();
            }
        }
        return out;
    };

    if (!dataMc) {
        _rlm = _rlm.Define("Jet_jer", applyJer, {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass", "GenJet_pt", "GenJet_eta", "GenJet_phi", "GenJet_mass", "Jet_genJetIdx", "Rho_fixedGridRhoFastjetAll", "event"});
    }
}

void NanoAODAnalyzerrdframe::skimJets() {
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
    _rlm = _rlm.Define("jetcuts", "Jet_pt>30.0 && abs(Jet_eta)<2.6 && Jet_passJetIdTightLepVeto && Jet_muEF<0.8 && Jet_chEmEF<0.8")
               .Redefine("Jet_pt", "Jet_pt[jetcuts]")
               .Redefine("Jet_eta", "Jet_eta[jetcuts]")
               .Redefine("Jet_phi", "Jet_phi[jetcuts]")
               .Redefine("Jet_mass", "Jet_mass[jetcuts]")
               .Redefine("Jet_passJetIdTightLepVeto", "Jet_passJetIdTightLepVeto[jetcuts]")
               .Redefine("Jet_area", "Jet_area[jetcuts]")
               .Redefine("Jet_rawFactor", "Jet_rawFactor[jetcuts]")
               //.Redefine("Jet_pt_uncorr", "Jet_pt_uncorr[jetcuts]")
               //.Redefine("Jet_rawFactor", "Jet_rawFactor[jetcuts]")
               .Redefine("Jet_btagPNetB", "Jet_btagPNetB[jetcuts]")
               .Redefine("nJet", "int(Jet_pt.size())");
    //if (!_isData) {
    //    _rlm = _rlm.Redefine("Jet_pt_unc", skimCol, {"Jet_pt_unc", "jetcuts"})
    //               .Redefine("Jet_jer", skimCol, {"Jet_jer", "jetcuts"})
    //               .Redefine("Jet_hadronFlavour","Jet_hadronFlavour[jetcuts]")
    //               .Redefine("Jet_genJetIdx","Jet_genJetIdx[jetcuts]");
    //}
}

void NanoAODAnalyzerrdframe::applyBSFs(std::vector<string> jes_var) {
    cout << "Loading Btag SF" << endl;
    std::string btagYear = "";
    if (_isRun22) {
        btagYear = "2022_Summer22";
    } else if (_isRun22EE) {
        btagYear = "2022_Summer22EE";
    } else if (_isRun23) {
        btagYear = "2023_Summer23";
    } else if (_isRun23BPix) {
        btagYear = "2023_Summer23BPix";
    } else if (_isRun24){
        //TODO
        btagYear = "2023_Summer23BPix";
    }

    _rlm = _rlm.Define("btag_var", [](){return strings(btag_var);})
               .Define("btag_jes_var", [jes_var](){return strings(jes_var);});

    //case3 - Shape correction
    //If you are interested in using the whole b-tagging discriminant distribution in your analysis,
    //e.g. using b-tagging variables to separate signal and background, then this method is for you
    auto bSFreader = correction::CorrectionSet::from_file("data/BTV/" + btagYear + "/btagging.json.gz");
    auto _PNetSF = bSFreader->at("particleNet_shape");
    auto btagweightgenerator= [this, _PNetSF](floats &pts, floats &etas, ints &hadflav, floats &btags)->floats{
        floats uncSources;
        uncSources.reserve(3);
        floats wVec;
        wVec.reserve(pts.size());
    
        for (auto i=0; i<int(pts.size()); i++){
            if(std::abs(etas[i])>2.5 || pts[i]<30) continue;
            uncSources.emplace_back(_PNetSF->evaluate({"central", hadflav[i], abs(etas[i]), pts[i], btags[i]}));
            if (hadflav[i] != 4) {
                uncSources.emplace_back(_PNetSF->evaluate({"up_hfstats2", hadflav[i], abs(etas[i]), pts[i], btags[i]}));
                uncSources.emplace_back(_PNetSF->evaluate({"down_hfstats2", hadflav[i], abs(etas[i]), pts[i], btags[i]}));
            }
            else {
                uncSources.emplace_back(_PNetSF->evaluate({"up_cferr1", hadflav[i], abs(etas[i]), pts[i], btags[i]}));
                uncSources.emplace_back(_PNetSF->evaluate({"down_cferr1", hadflav[i], abs(etas[i]), pts[i], btags[i]}));
            }
            //wVec.emplace_back(uncSources);
            uncSources.clear();
        }
        return wVec;
    };
    
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

void NanoAODAnalyzerrdframe::selectJets(std::vector<std::string> jes_var, std::vector<std::string> jes_var_flav) {
    //// input vector: vec[pt][vars], for bSF
    //auto skimCol = [this](doublesVec toSkim, ints cut)->doublesVec {

    //    doublesVec out;
    //    for (size_t i=0; i<toSkim.size(); i++) {
    //        if (cut[i] > 0) out.emplace_back(toSkim[i]);
    //    }
    //    return out;
    //};

    //// input vector: vec[pt][vars]
    //auto calcBSF = [this](doublesVec perJetSF, int nvar)->doubles {

    //    doubles out;
    //    out.reserve(nvar);
    //    for (size_t i=0; i<nvar; i++) {
    //        double bSF = 1.0;
    //        for (size_t j=0; j<perJetSF.size(); j++) {
    //            if (perJetSF[j].empty()) continue;
    //            bSF *= perJetSF[j][i];
    //        }
    //        out.emplace_back(bSF);
    //    }
    //    return out;
    //};

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
    _rlm = _rlm.Define("jetcuts", "Jet_pt>40.0 && abs(Jet_eta)<2.4 && Jet_passJetIdTightLepVeto && Jet_muEF < 0.8 && Jet_chEmEF < 0.8");
    _rlm = _rlm.Redefine("Jet_pt", "Jet_pt[jetcuts]")
               .Redefine("Jet_eta", "Jet_eta[jetcuts]")
               .Redefine("Jet_phi", "Jet_phi[jetcuts]")
               .Redefine("Jet_mass", "Jet_mass[jetcuts]")
               .Redefine("Jet_btagPNetB", "Jet_btagPNetB[jetcuts]")
               .Define("jet4vecs", ::gen4vec, {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass"});

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
    _rlm = _rlm.Define("lepjetoverlap", checkoverlap, {"jet4vecs","lep4vecs"});

    _rlm = _rlm.Define("taujetoverlap", checkoverlap, {"jet4vecs","cleantau4vecs"})
               .Define("loosetaujetoverlap", checkoverlap, {"jet4vecs","cleanloosetau4vecs"})
               .Define("jetoverlap","lepjetoverlap && taujetoverlap")
               .Define("jetoverlaploose","lepjetoverlap && loosetaujetoverlap");

    _rlm = _rlm.Define("Jet_pt_loose", "Jet_pt[jetoverlaploose]")
               .Define("Jet_btagPNetB_loose", "Jet_btagPNetB[jetoverlaploose]")
               .Define("ncleanjetsloosepass", "int(Jet_pt_loose.size())");

    _rlm = _rlm.Redefine("Jet_pt", "Jet_pt[jetoverlap]")
               .Redefine("Jet_eta", "Jet_eta[jetoverlap]")
               .Redefine("Jet_phi", "Jet_phi[jetoverlap]")
               .Redefine("Jet_mass", "Jet_mass[jetoverlap]")
               .Redefine("Jet_btagPNetB", "Jet_btagPNetB[jetoverlap]")
               .Define("ncleanjetspass", "int(Jet_pt.size())")
               .Define("cleanjet4vecs", ::gen4vec, {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass"})
               .Define("Jet_HT", "Sum(Jet_pt)");

    //if (!_isData) {
    //    int nbsf_var = btag_var.size();
    //    int njes_var = jes_var.size();
    //    _rlm = _rlm.Define("nbsf_var", [nbsf_var](){return int(nbsf_var);})
    //               .Define("njes_var", [njes_var](){return int(njes_var);})
    //               .Define("btagWeight_PNetB_perJet_loose", skimCol, {"btagWeight_PNetB_perJet", "jetoverlaploose"})
    //               .Define("btagWeight_PNetB_loose", calcBSF, {"btagWeight_PNetB_perJet_loose", "nbsf_var"})
    //               .Redefine("btagWeight_PNetB_perJet", skimCol, {"btagWeight_PNetB_perJet", "jetoverlap"})
    //               .Redefine("btagWeight_PNetB_jes_perJet", skimCol, {"btagWeight_PNetB_jes_perJet", "jetoverlap"})
    //               .Define("btagWeight_PNetB", calcBSF, {"btagWeight_PNetB_perJet", "nbsf_var"})
    //               .Define("btagWeight_PNetB_jes", calcBSF, {"btagWeight_PNetB_jes_perJet", "njes_var"});
    //}


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
        _rlm = _rlm.Define("btagcuts", "Jet_btagPNetB>0.1990") //l: 0.0365, m: 0.1990, t: 0.6373
                   .Define("btagcuts_loose", "Jet_btagPNetB>0.0365");
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
}

void NanoAODAnalyzerrdframe::selectTaus() {

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
    _rlm = _rlm.Define("taucuts", "Tau_pt>30.0 && abs(Tau_eta)<2.5  && Tau_idDecayModeNewDMs && Tau_decayMode != 5 && Tau_decayMode != 6")
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

        std::string tauYear = "";
        if (_isRun22) {
            tauYear = "2022_preEE";
        } else if (_isRun22EE) {
            tauYear = "2022_postEE";
        } else if (_isRun23) {
            tauYear = "2023_preBPix";
        } else if (_isRun23BPix) {
            tauYear = "2023_postBPix";
        } else if (_isRun24){
            //TODO
            tauYear = "2023_postBPix";
        }

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

            //std::cout << "size eta: " << eta.size() << " dm: " << dm.size() << " genmatch: " << genmatch.size() << std::endl;
            if (eta.size() > 0) {
                for (unsigned int i=0; i<eta.size(); i++) {
                    //std::cout << "eta: " << eta[i] << " dm: " << int(dm[i]) << " genmatch: " << int(genmatch[i]) << std::endl;
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
    _rlm = _rlm.Define("GenPart_LFVup_idx", "FinalGenPart_idx[0]")
               .Define("GenPart_LFVmuon_idx", "FinalGenPart_idx[1]")
               .Define("GenPart_LFVtau_idx", "FinalGenPart_idx[2]")
               .Define("GenPart_SMb_idx", "FinalGenPart_idx[3]")
               .Define("GenPart_SMW1_idx", "FinalGenPart_idx[4]")
               .Define("GenPart_SMW2_idx", "FinalGenPart_idx[5]")
               .Define("GenPart_LFVtop_idx", "FinalGenPart_idx[6]")
               .Define("GenPart_SMtop_idx", "FinalGenPart_idx[7]");

    _rlm = _rlm.Define("drmax1", "float(0.15)")
               .Define("drmax2", "float(0.4)");
    if (_isMuonCh){
        _rlm = _rlm.Define("Muon_matched", ::dRmatching_binary,{"GenPart_LFVmuon_idx","drmax1","GenPart_pt","GenPart_eta","GenPart_phi","GenPart_mass","Muon_pt","Muon_eta","Muon_phi","Muon_mass"});
    } else {
        _rlm = _rlm.Define("Muon_matched", ::dRmatching_binary,{"GenPart_LFVmuon_idx","drmax1","GenPart_pt","GenPart_eta","GenPart_phi","GenPart_mass","Electron_pt","Electron_eta","Electron_phi","Electron_mass"});
    }
    _rlm = _rlm.Define("Tau_matched",::dRmatching_binary,{"GenPart_LFVtau_idx","drmax2", "GenPart_pt","GenPart_eta","GenPart_phi","GenPart_mass","Tau_pt","Tau_eta","Tau_phi","Tau_mass"})
               .Define("Jet_LFVup_matched",::dRmatching_binary,{"GenPart_LFVup_idx","drmax2","GenPart_pt","GenPart_eta","GenPart_phi","GenPart_mass","Jet_pt","Jet_eta","Jet_phi","Jet_mass"})
               .Define("Jet_SMb_matched",::dRmatching_binary,{"GenPart_SMb_idx","drmax2","GenPart_pt","GenPart_eta","GenPart_phi","GenPart_mass","Jet_pt","Jet_eta","Jet_phi","Jet_mass"})
               .Define("Jet_SMW1_matched",::dRmatching_binary,{"GenPart_SMW1_idx","drmax2","GenPart_pt","GenPart_eta","GenPart_phi","GenPart_mass","Jet_pt","Jet_eta","Jet_phi","Jet_mass"})
               .Define("Jet_SMW2_matched",::dRmatching_binary,{"GenPart_SMW2_idx","drmax2", "GenPart_pt","GenPart_eta","GenPart_phi","GenPart_mass","Jet_pt","Jet_eta","Jet_phi","Jet_mass"});
    if (_isMuonCh){
        _rlm = _rlm.Define("Sel_muon_matched","Muon_matched[muoncuts]");
    } else {
        _rlm = _rlm.Define("Sel_muon_matched", "Muon_matched[elecuts]");
    }
    _rlm = _rlm.Define("Sel_tau_matched","Tau_matched[seltaucuts]")
               .Define("Sel2_LFVupjet_matched","Jet_LFVup_matched[jetcuts][lepjetoverlap][taujetoverlap]")
               .Define("Sel2_SMbjet_matched","Jet_SMb_matched[jetcuts][lepjetoverlap][taujetoverlap]")
               .Define("Sel2_SMW1jet_matched","Jet_SMW1_matched[jetcuts][lepjetoverlap][taujetoverlap]")
               .Define("Sel2_SMW2jet_matched","Jet_SMW2_matched[jetcuts][lepjetoverlap][taujetoverlap]")
               .Define("Sel_muonmatched_idx","ArgMax(Sel_muon_matched)")
               .Define("Sel_taumatched_idx","ArgMax(Sel_tau_matched)")
               .Define("Sel2_LFVupjet_matched_idx","ArgMax(Sel2_LFVupjet_matched)")
               .Define("Sel2_SMbjet_matched_idx","ArgMax(Sel2_SMbjet_matched)")
               .Define("Sel2_SMW1jet_matched_idx","ArgMax(Sel2_SMW1jet_matched)")
               .Define("Sel2_SMW2jet_matched_idx","ArgMax(Sel2_SMW2jet_matched)")
               .Define("nmuonmatched", "Sum(Muon_matched)")
               .Define("ntaumatched", "Sum(Tau_matched)")
               .Define("nJet_LFVup_matched", "Sum(Jet_LFVup_matched)")
               .Define("nJet_SMb_matched", "Sum(Jet_SMb_matched)")
               .Define("nJet_SMW1_matched", "Sum(Jet_SMW1_matched)")
               .Define("nJet_SMW2_matched", "Sum(Jet_SMW2_matched)");
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

void NanoAODAnalyzerrdframe::calculateEvWeight() {

    // Put SFs that need to be calculated in skim:
    // e.g. tau selection is done at processing but we want to have SFs already

    // Tau SF
    cout << "Loading Tau SF" << endl;
    std::string tauid_vsjet = "VTight";
    std::string tauid_vse   = "Tight";
    std::string tauid_vsmu  = "Tight";

    cout << "Tau ID WP vsJet : " << tauid_vsjet << endl;
    cout << "Tau ID WP vsMuon : " << tauid_vsmu << endl;
    cout << "Tau ID WP vsElectron : " << tauid_vse << endl;

    std::string tauYear = "";
    if (_isRun22) {
        tauYear = "2022_preEE";
    } else if (_isRun22EE) {
        tauYear = "2022_postEE";
    } else if (_isRun23) {
        tauYear = "2023_preBPix";
    } else if (_isRun23BPix) {
        tauYear = "2023_postBPix";
    } else if (_isRun24) {
        //TODO
        tauYear = "2023_postBPix";
    }

    auto tauSFreader = correction::CorrectionSet::from_file("data/TauIDSFs/tau_DeepTau2018v2p5_" + tauYear + ".json.gz");
    auto _testool    = tauSFreader->at("tau_energy_scale");

    // Tau ES
    cout<<"Applying TauES on Genuine taus"<<endl;
    auto tauES = [this, _testool, tauid_vsjet, tauid_vse](floats &pt, floats &eta, uchars &dm, uchars &genid, floats &x)->floats {

        floats xout;

        for (unsigned int i=0; i<pt.size(); i++) {
            float es = 1.0;
            if (int(genid[i])==1 || int(genid[i])==3 || int(genid[i])==5) {
                if (dm[i]!=5 and dm[i]!=6)
                    es = _testool->evaluate({pt[i], eta[i], dm[i], int(genid[i]), "DeepTau2018v2p5", tauid_vsjet, tauid_vse, "nom"});
            }
            xout.emplace_back(x[i]*es);
        }
        return xout;
    };

    auto tauESUnc = [this, _testool, tauid_vsjet, tauid_vse](floats &pt, floats &eta, uchars &dm, uchars &genid, floats &x)->floatsVec {

        floats uncSources;
        uncSources.reserve(2);
        floatsVec xout;
        xout.reserve(pt.size());

        for (unsigned int i=0; i<pt.size(); i++) {
          if (int(genid[i])==1 || int(genid[i])==3 || int(genid[i])==5) {
                if (dm[i]!=5 and dm[i]!=6) {
                    uncSources.emplace_back(_testool->evaluate({pt[i], eta[i], dm[i], int(genid[i]), "DeepTau2018v2p5", tauid_vsjet, tauid_vse, "up"}));
                    uncSources.emplace_back(_testool->evaluate({pt[i], eta[i], dm[i], int(genid[i]), "DeepTau2018v2p5", tauid_vsjet, tauid_vse, "down"}));
                }
            }
            else uncSources = {1.0f, 1.0f};
            xout.emplace_back(uncSources);
            uncSources.clear();
        }
        return xout;
    };

    _rlm = _rlm.Define("Tau_pt_uncor", "Tau_pt")
           .Redefine("Tau_pt", tauES, {"Tau_pt_uncor", "Tau_eta", "Tau_decayMode", "Tau_genPartFlav", "Tau_pt_uncor"})
           .Redefine("Tau_mass", tauES, {"Tau_pt_uncor", "Tau_eta", "Tau_decayMode", "Tau_genPartFlav", "Tau_mass"})
           .Define("Tau_pt_unc", tauESUnc, {"Tau_pt_uncor", "Tau_eta", "Tau_decayMode", "Tau_genPartFlav", "Tau_pt_uncor"});


    // ID SFs
    /*
    auto tauSFIdVsJet = [this, _tauidSFjet, tauYear, tauid_vsjet, tauid_vse](floats &pt, floats &eta, uchars &genid, ints &dm)->floatsVec {

        floats uncSources;
        uncSources.reserve(11);
        floatsVec wVec;
        wVec.reserve(pt.size());

        auto year_ = tauYear.substr(2);
        //std::vector<std::string> uncerts = {"uncert0", "uncert1", "syst_alleras", "syst_" + year_, "syst_dmX_" + year_};

        if (pt.size() > 0) {
            for (unsigned int i=0; i<pt.size(); i++) {
                // TauSFTool will take care of pt > 140 SF by setting pT = 140
                uncSources.emplace_back(_tauidSFjet->evaluate({pt[i], dm[i], int(genid[i]), tauid_vsjet, tauid_vse, "nom", "pt"}));

                //for (auto unc : uncerts) {
                //    size_t pos = unc.find("dmX");
                //    if (pos != std::string::npos) unc.replace(pos, 3, "dm"+std::to_string(dm[i]));
                //    if (pt[i] <= 140) {
                //        uncSources.emplace_back(_tauidSFjet->evaluate({pt[i], dm[i], int(genid[i]), tauid_vsjet, tauid_vse, unc + "_up", "pt"}));
                //        uncSources.emplace_back(_tauidSFjet->evaluate({pt[i], dm[i], int(genid[i]), tauid_vsjet, tauid_vse, unc + "_down", "pt"}));
                //    } else {
                //        uncSources.emplace_back(1.0 + (min(pt[i], static_cast<float>(500.)) - 40.) * (pt[i] > 40.) * 0.00018);
                //        uncSources.emplace_back(1.0 - (min(pt[i], static_cast<float>(500.)) - 40.) * (pt[i] > 40.) * 0.00018);
                //    }
                //}
                wVec.emplace_back(uncSources);
                uncSources.clear();
            }
        }
        return wVec;
    };
    
    TauIDSFTool* _tauidSFjet;
    TauIDSFTool* _tauidSFjetHighPt;
    _tauidSFjet = new TauIDSFTool(tauYear, "DeepTau2017v2p1VSjet", tauid_vsjet, tauid_vse, false, true, false, false);
    _tauidSFjetHighPt = new TauIDSFTool(tauYear, "DeepTau2017v2p1VSjet", tauid_vsjet, tauid_vse, false, false, false, true);

    auto tauSFIdVsJet = [this, _tauidSFjet, _tauidSFjetHighPt, tauYear](floats &pt, floats &eta, uchars &genid, ints dm)->floatsVec {

        floats uncSources;
        uncSources.reserve(27); //nom + syst 40 + highPT 2
        floatsVec wVec;
        wVec.reserve(pt.size());

        auto year_ = tauYear.substr(2);
        std::vector<std::string> uncerts = {"uncert0", "uncert1", "syst_alleras", "syst_" + year_, "syst_dmX_" + year_};
        std::vector<std::string> uncertsHighPt = {"stat","stat_bin1","stat_bin2","syst","extrap"};

        if (pt.size() > 0) {
            for (unsigned int i=0; i<pt.size(); i++) {
                // TauSFTool will take care of pt > 140 SF by setting pT = 140
                float nomsf = _tauidSFjet->getSFvsDMandPT(pt[i], dm[i], int(genid[i]));
                float nomsf_highpt = _tauidSFjetHighPt->getHighPTSFvsPT(pt[i], int(genid[i]));

                if (pt[i] <= 140) {
                    uncSources.emplace_back(nomsf);
                    for (auto unc : uncerts) {
                        size_t pos = unc.find("dmX");
                        if (pos != std::string::npos) { //indices 9-16
                            unc.replace(pos, 3, "dm"+std::to_string(dm[i]));
                            float upsf = _tauidSFjet->getSFvsDMandPT(pt[i], dm[i], int(genid[i]), unc + "_up");
                            float dnsf = _tauidSFjet->getSFvsDMandPT(pt[i], dm[i], int(genid[i]), unc + "_down");
                            std::vector<float> dmsf;
                            if      (dm[i] == 0)  dmsf = {upsf, dnsf, nomsf, nomsf, nomsf, nomsf, nomsf, nomsf};
                            else if (dm[i] == 1)  dmsf = {nomsf, nomsf, upsf, dnsf, nomsf, nomsf, nomsf, nomsf};
                            else if (dm[i] == 10) dmsf = {nomsf, nomsf, nomsf, nomsf, upsf, dnsf, nomsf, nomsf};
                            else if (dm[i] == 11) dmsf = {nomsf, nomsf, nomsf, nomsf, nomsf, nomsf, upsf, dnsf};
                            else                  dmsf = {nomsf, nomsf, nomsf, nomsf, nomsf, nomsf, nomsf, nomsf}; //placeholder
                            uncSources.insert(uncSources.end(), dmsf.begin(), dmsf.end());
                        } else { //indices 1-8
                            uncSources.emplace_back(_tauidSFjet->getSFvsDMandPT(pt[i], dm[i], int(genid[i]), unc + "_up"));
                            uncSources.emplace_back(_tauidSFjet->getSFvsDMandPT(pt[i], dm[i], int(genid[i]), unc + "_down"));
                        }
                    }
                    uncSources.insert(uncSources.end(), 10, nomsf);
                } else {
                    uncSources.emplace_back(nomsf_highpt);
                    uncSources.insert(uncSources.end(), 16, nomsf_highpt);
                    for (auto unc : uncertsHighPt) {
                        uncSources.emplace_back(_tauidSFjetHighPt->getHighPTSFvsPT(pt[i], int(genid[i]), unc + "_up"));
                        uncSources.emplace_back(_tauidSFjetHighPt->getHighPTSFvsPT(pt[i], int(genid[i]), unc + "_down"));
                    }
                }
                wVec.emplace_back(uncSources);
                uncSources.clear();
            }
        }
        return wVec;
    };*/

    /*
    auto tauSFIdVsJet = [this, _tauidSFjet, tauid_vsjet, tauid_vse](floats &pt, uchars &dm, uchars &genmatch)->floatsVec {

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

        std::cout << "tauSFIdVsEl" << std::endl;
        floats uncSources;
        uncSources.reserve(3);
        floatsVec wVec;
        wVec.reserve(eta.size());

        std::cout << "size eta: " << eta.size() << " dm: " << dm.size() << " genmatch: " << genmatch.size() << std::endl;
        if (eta.size() > 0) {
            for (unsigned int i=0; i<eta.size(); i++) {
                std::cout << "eta: " << eta[i] << " dm: " << int(dm[i]) << " genmatch: " << int(genmatch[i]) << std::endl;
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

    //tauid_vsjet = "Loose";
    //_tauidSFjet = new TauIDSFTool(tauYear, "DeepTau2018v2p5VSjet", tauid_vsjet, tauid_vse, false, true, false, false);
    //_tauidSFjetHighPt = new TauIDSFTool(tauYear, "DeepTau2018v2p5VSjet", tauid_vsjet, tauid_vse, false, false, false, true);

    //_rlm = _rlm.Define("tauWeightIdVsJet_loose", tauSFIdVsJet, {"Tau_pt","Tau_decayMode","Tau_genPartFlav"});
    */
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
            arnode->Snapshot(outtreename, outname);
        } else {
            // use the following if you want to store only a few variables
            //arnode->Snapshot(outtreename, outname, _varstostore);
            cout << " writing branches" << endl;
            for (auto bname: _varstostorepertree[nodename]) {
                cout << bname << ", ";
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

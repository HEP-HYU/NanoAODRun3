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

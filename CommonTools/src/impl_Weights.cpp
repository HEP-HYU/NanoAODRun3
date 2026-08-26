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

void NanoAODAnalyzerrdframe::applyWeights(string pileFile, string map){
    // Event weight for data it's always one. For MC, it depends on the sign
    _rlm = _rlm.Define("lhereweight","one");
    _rlm = _rlm.Define("unitGenWeight", "one");

    _rlm = _rlm.Define("isData", "true");
    if(!_isData){
        _rlm = _rlm.Redefine("isData", "false");


        // Store sum of weights
        auto storePDFWeights = [this](floats weights, float gen)->floats {
            for (size_t i=0; i<weights.size(); i++)
                PDFWeights[i] += (gen / abs(gen)) * weights[i];
            return PDFWeights;
        };
        auto storePSWeights = [this](floats weights, float gen)->floats {
            for (size_t i=0; i<weights.size(); i++) {
                if (i > 3) continue; //JME Nano stores all PS
                PSWeights[i] += (gen / abs(gen)) * weights[i];
            }
            return PSWeights;
        };
        auto storeScaleWeights = [this](floats weights, float gen)->floats {
            for (size_t i=0; i<weights.size(); i++)
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
        if (_outfilename.find("WtoLNu") != std::string::npos || _outfilename.find("WtoENu") != std::string::npos
            || _outfilename.find("WtoMuNu") != std::string::npos || _outfilename.find("WtoTauNu") != std::string::npos) {
            std::cout << "WtoLNu" << std::endl;
            _rlm = _rlm.Redefine("lhereweight","LHEWeight_originalXWGTUP/abs(LHEWeight_originalXWGTUP)");
            _rlm = _rlm.Redefine("unitGenWeight","LHEWeight_originalXWGTUP/abs(LHEWeight_originalXWGTUP)");
        };

        
        auto puWeightreader = loadCorrectionSet("data/LUM/"+pileFile+"/puWeights.json.gz");
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


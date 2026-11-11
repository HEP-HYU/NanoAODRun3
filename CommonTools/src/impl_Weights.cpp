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
        bool hasPDF   = _rd.HasColumn("LHEPdfWeight");
        bool hasPS    = _rd.HasColumn("PSWeight");
        bool hasScale = _rd.HasColumn("LHEScaleWeight");
        bool hasGen   = _rd.HasColumn("genWeight");

        if (!hasGen) {
            std::cerr << "[ERROR] applyWeights: Missing genWeight branch in MC sample!" << std::endl;
            throw std::runtime_error("[ERROR] applyWeights: Missing genWeight branch in MC sample!");
        }

        // Case A: Standard NLO MC with all 3 weight branches (TTbar, ST, DY, WJets, Signal)
        // Consolidated into a single event loop to reduce I/O overhead by 3x while preserving exact physics.
        if (hasPDF && hasPS && hasScale) {
            auto storeAllWeights = [this](const floats& pdfw, const floats& psw, const floats& scalew, float gen) {
                float sign = (gen != 0.0f) ? (gen / std::abs(gen)) : 0.0f;
                for (size_t i = 0; i < std::min(pdfw.size(), PDFWeights.size()); i++)
                    PDFWeights[i] += sign * pdfw[i];
                for (size_t i = 0; i < std::min(psw.size(), PSWeights.size()); i++) {
                    if (i > 3) continue; // JME Nano stores all PS
                    PSWeights[i] += sign * psw[i];
                }
                for (size_t i = 0; i < std::min(scalew.size(), ScaleWeights.size()); i++)
                    ScaleWeights[i] += sign * scalew[i];
            };
            try {
                _rlm.Foreach(storeAllWeights, {"LHEPdfWeight", "PSWeight", "LHEScaleWeight", "genWeight"});
            } catch (const std::exception& e) {
                std::cerr << "[ERROR] applyWeights: Failed executing consolidated weight loop: " << e.what() << std::endl;
                throw;
            }
        }
        // Case B: Selective processing for samples missing one or more theory weight branches (e.g. LO Pythia samples)
        else {
            auto storePDFWeights = [this](const floats& weights, float gen) {
                float sign = (gen != 0.0f) ? (gen / std::abs(gen)) : 0.0f;
                for (size_t i = 0; i < std::min(weights.size(), PDFWeights.size()); i++)
                    PDFWeights[i] += sign * weights[i];
            };
            auto storePSWeights = [this](const floats& weights, float gen) {
                float sign = (gen != 0.0f) ? (gen / std::abs(gen)) : 0.0f;
                for (size_t i = 0; i < std::min(weights.size(), PSWeights.size()); i++) {
                    if (i > 3) continue;
                    PSWeights[i] += sign * weights[i];
                }
            };
            auto storeScaleWeights = [this](const floats& weights, float gen) {
                float sign = (gen != 0.0f) ? (gen / std::abs(gen)) : 0.0f;
                for (size_t i = 0; i < std::min(weights.size(), ScaleWeights.size()); i++)
                    ScaleWeights[i] += sign * weights[i];
            };

            if (hasPDF) {
                try {
                    _rlm.Foreach(storePDFWeights, {"LHEPdfWeight", "genWeight"});
                } catch (const std::exception& e) {
                    std::cerr << "[ERROR] applyWeights: Failed reading LHEPdfWeight: " << e.what() << std::endl;
                }
            } else {
                std::cout << "[INFO] applyWeights: No LHEPdfWeight branch in this ROOT file." << std::endl;
            }

            if (hasPS) {
                try {
                    _rlm.Foreach(storePSWeights, {"PSWeight", "genWeight"});
                } catch (const std::exception& e) {
                    std::cerr << "[ERROR] applyWeights: Failed reading PSWeight: " << e.what() << std::endl;
                }
            } else {
                std::cout << "[INFO] applyWeights: No PSWeight branch in this ROOT file." << std::endl;
            }

            if (hasScale) {
                try {
                    _rlm.Foreach(storeScaleWeights, {"LHEScaleWeight", "genWeight"});
                } catch (const std::exception& e) {
                    std::cerr << "[ERROR] applyWeights: Failed reading LHEScaleWeight: " << e.what() << std::endl;
                }
            } else {
                std::cout << "[INFO] applyWeights: No LHEScaleWeight branch in this ROOT file." << std::endl;
            }
        }

        ////Check Normalisation issue for genWeight
        _rlm = _rlm.Redefine("unitGenWeight","genWeight != 0 ? genWeight/abs(genWeight) : 0");
        if (_outfilename.find("WtoLNu") != std::string::npos || _outfilename.find("WtoENu") != std::string::npos
            || _outfilename.find("WtoMuNu") != std::string::npos || _outfilename.find("WtoTauNu") != std::string::npos) {
            std::cout << "WtoLNu" << std::endl;
            _rlm = _rlm.Redefine("lhereweight","LHEWeight_originalXWGTUP/abs(LHEWeight_originalXWGTUP)");
            _rlm = _rlm.Redefine("unitGenWeight","LHEWeight_originalXWGTUP/abs(LHEWeight_originalXWGTUP)");
        };

        
        const std::string puJsonPath = (pileFile.find(".json.gz") != std::string::npos)
            ? ("data/LUM/" + pileFile)
            : ("data/LUM/" + pileFile + "/puWeights.json.gz");
        std::shared_ptr<correction::CorrectionSet> puWeightreader;
        correction::Correction::Ref _puweight;
        try {
            puWeightreader = loadCorrectionSet(puJsonPath);
            _puweight = puWeightreader->at(map);
        } catch (const std::exception& e) {
            std::string msg = "[ERROR] applyWeights: Failed to load pileup correction '" + map +
                              "' from " + puJsonPath + ": " + e.what();
            std::cerr << msg << std::endl;
            throw std::runtime_error(msg);
        }

        auto PuWeight = [this, _puweight, map, puJsonPath](float x) -> floats {
            floats out;
            try {
                out.emplace_back(_puweight->evaluate({x, "nominal"}));
                out.emplace_back(_puweight->evaluate({x, "up"}));
                out.emplace_back(_puweight->evaluate({x, "down"}));
            } catch (const std::exception& e) {
                std::string msg = "[ERROR] applyWeights: Pileup weight evaluation failed for '" + map +
                                  "' (" + puJsonPath + ") at Pileup_nTrueInt=" + std::to_string(x) + ": " + e.what();
                std::cerr << msg << std::endl;
                throw std::runtime_error(msg);
            }

            return out;
        };

        _rlm = _rlm.Define("puWeight", PuWeight, {"Pileup_nTrueInt"});
    }
}


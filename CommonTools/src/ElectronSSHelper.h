#ifndef ElectronSSHelper_h
#define ElectronSSHelper_h

#include <cmath>
#include <vector>
#include <string>
#include <memory>
#include <iostream>
#include <algorithm>
#include "TMath.h"
#include "correction.h"

// Helper class encapsulating CMS EGM POG Energy Scale & Resolution Smearing (S&S)
// Following official CMS EGM Run 3 recommendation (CMS_EgammSFandSSRun3.pdf).
class ElectronSSHelper {
public:
    ElectronSSHelper(std::shared_ptr<correction::CorrectionSet> cset)
        : _cset(cset)
    {
        if (_cset) {
            try {
                _smearCorr = _cset->at("SmearAndSyst");
            } catch (...) {
                _smearCorr = nullptr;
            }

            for (const auto& name : {"EGMScale_ElePTsplit_2022preEE", "EGMScale_ElePTsplit_2022postEE",
                                     "EGMScale_ElePTsplit_2023preBPIX", "EGMScale_ElePTsplit_2023postBPIX",
                                     "EGMScale_ElePTsplit_2024"}) {
                try {
                    _dataScaleCorr = _cset->at(name);
                    if (_dataScaleCorr) break;
                } catch (...) {
                    _dataScaleCorr = nullptr;
                }
            }
        }
    }

    // Deterministic Gaussian pseudo-random generator based on 64-bit event hash
    static inline double get_gaussian_rndm(unsigned long long evt, unsigned int lumi, double phi, int idx) {
        uint64_t z = evt ^ (static_cast<uint64_t>(lumi) * 0x9e3779b97f4a7c15ULL)
                         ^ (static_cast<uint64_t>(std::abs(phi) * 1.0e6) * 0xbf58476d1ce4e5b9ULL)
                         ^ (static_cast<uint64_t>(idx + 1) * 0x94d049bb133111ebULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        z = z ^ (z >> 31);
        double u = (z >> 11) * 0x1.0p-53;
        if (u <= 1e-7) u = 1e-7;
        if (u >= 1.0 - 1e-7) u = 1.0 - 1e-7;
        return std::sqrt(2.0) * TMath::ErfInverse(2.0 * u - 1.0);
    }

    double pt_corr(bool is_data, const std::string& syst, double pt, double eta, double phi, double deltaEtaSC, float r9,
                   unsigned long long evt, unsigned int lumi, int idx) const {
        // Per EGM POG prescription: scale and smearing corrections are not tuned below 15 GeV
        if (pt < 15.0) return pt;

        double scEta = eta + deltaEtaSC;
        double r9_val = static_cast<double>(r9);

        if (is_data) {
            if (_dataScaleCorr) {
                double scale = _dataScaleCorr->evaluate({std::string("scale"), pt, r9_val, scEta});
                return pt * scale;
            }
            return pt;
        }

        if (!_smearCorr) return pt;

        double g = get_gaussian_rndm(evt, lumi, phi, idx);

        double pt_res = pt;
        if (syst.find("elesmearup") != std::string::npos) {
            double sigma_up = _smearCorr->evaluate({std::string("smear_up"), pt, r9_val, scEta});
            pt_res = pt * (1.0 + g * sigma_up);
        } else if (syst.find("elesmeardown") != std::string::npos) {
            double sigma_down = _smearCorr->evaluate({std::string("smear_down"), pt, r9_val, scEta});
            pt_res = pt * (1.0 + g * sigma_down);
        } else if (syst.find("elescaleup") != std::string::npos) {
            // Apply scale systematic variation on top of nominal smeared MC (per EGM POG recommendation)
            double sigma = _smearCorr->evaluate({std::string("smear"), pt, r9_val, scEta});
            double pt_smeared = pt * (1.0 + g * sigma);
            double scale_up = _smearCorr->evaluate({std::string("scale_up"), pt, r9_val, scEta});
            pt_res = pt_smeared * scale_up;
        } else if (syst.find("elescaledown") != std::string::npos) {
            // Apply scale systematic variation on top of nominal smeared MC
            double sigma = _smearCorr->evaluate({std::string("smear"), pt, r9_val, scEta});
            double pt_smeared = pt * (1.0 + g * sigma);
            double scale_down = _smearCorr->evaluate({std::string("scale_down"), pt, r9_val, scEta});
            pt_res = pt_smeared * scale_down;
        } else {
            // Nominal MC stochastic resolution smearing
            double sigma = _smearCorr->evaluate({std::string("smear"), pt, r9_val, scEta});
            pt_res = pt * (1.0 + g * sigma);
        }

        if (pt_res < 0.0) pt_res = 0.0;
        return pt_res;
    }

private:
    std::shared_ptr<correction::CorrectionSet> _cset;
    std::shared_ptr<const correction::Correction> _smearCorr;
    std::shared_ptr<const correction::Correction> _dataScaleCorr;
};

#endif

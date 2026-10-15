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
    ElectronSSHelper(std::shared_ptr<correction::CorrectionSet> cset, const std::string& scaleCorrName)
        : _cset(cset)
    {
        if (_cset) {
            _smearCorr = _cset->at("SmearAndSyst");
            if (!scaleCorrName.empty()) {
                _dataScaleCorr = _cset->at(scaleCorrName);
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
        double g = std::sqrt(2.0) * TMath::ErfInverse(2.0 * u - 1.0);
        return g;
    }

    double pt_corr(bool is_data, const std::string& syst, double pt, double eta, double phi, double deltaEtaSC, float r9,
                   unsigned long long evt, unsigned int lumi, int idx) const {
        // Per EGM POG prescription: scale and smearing corrections are not calibrated below 15 GeV
        if (pt < 15.0) return pt;

        double scEta = eta + deltaEtaSC;
        double r9_val = static_cast<double>(r9);
        double pt_val = pt;

        if (is_data) {
            if (_dataScaleCorr) {
                double scale = _dataScaleCorr->evaluate({scEta, r9_val});
                return pt * scale;
            }
            return pt;
        }

        if (!_smearCorr) return pt;

        double g = get_gaussian_rndm(evt, lumi, phi, idx);

        if (syst.find("elesmearup") != std::string::npos) {
            double sigma_up = _smearCorr->evaluate({std::string("smear_up"), pt_val, r9_val, scEta});
            return pt * (1.0 + g * sigma_up);
        } else if (syst.find("elesmeardown") != std::string::npos) {
            double sigma_down = _smearCorr->evaluate({std::string("smear_down"), pt_val, r9_val, scEta});
            return pt * (1.0 + g * sigma_down);
        } else if (syst.find("elescaleup") != std::string::npos) {
            double sigma = _smearCorr->evaluate({std::string("smear"), pt_val, r9_val, scEta});
            double pt_smeared = pt * (1.0 + g * sigma);
            double scale_up = _smearCorr->evaluate({std::string("scale_up"), pt_val, r9_val, scEta});
            return pt_smeared * scale_up;
        } else if (syst.find("elescaledown") != std::string::npos) {
            double sigma = _smearCorr->evaluate({std::string("smear"), pt_val, r9_val, scEta});
            double pt_smeared = pt * (1.0 + g * sigma);
            double scale_down = _smearCorr->evaluate({std::string("scale_down"), pt_val, r9_val, scEta});
            return pt_smeared * scale_down;
        } else {
            double sigma = _smearCorr->evaluate({std::string("smear"), pt_val, r9_val, scEta});
            return pt * (1.0 + g * sigma);
        }
    }

private:
    std::shared_ptr<correction::CorrectionSet> _cset;
    std::shared_ptr<const correction::Correction> _smearCorr;
    std::shared_ptr<const correction::Correction> _dataScaleCorr;
};

#endif

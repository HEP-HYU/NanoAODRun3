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
// Following official CMS EGM Run 3 recommendation (CMS_egmScalesAndSmearingsRun3.pdf & CMS_egmScaleAndSmearingExample.py).
class ElectronSSHelper {
public:
    ElectronSSHelper(std::shared_ptr<correction::CorrectionSet> cset)
        : _cset(cset)
    {
        if (!_cset) {
            std::cerr << "[ERROR] ElectronSSHelper: CorrectionSet is null!" << std::endl;
        } else {
            _smearCorr = _cset->at("SmearAndSyst");
            _dataCompoundScaleCorr = _cset->compound().at("Scale");
        }
    }

    // Deterministic Gaussian pseudo-random generator based on 64-bit event and electron phi hash
    static inline double get_gaussian_rndm(unsigned long long evt, unsigned int lumi, double phi) {
        uint64_t z = evt ^ (static_cast<uint64_t>(lumi) * 0x9e3779b97f4a7c15ULL)
                         ^ (static_cast<uint64_t>(std::abs(phi) * 1.0e6) * 0xbf58476d1ce4e5b9ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        z = z ^ (z >> 31);
        double u = (z >> 11) * 0x1.0p-53;
        if (u <= 1e-7) u = 1e-7;
        if (u >= 1.0 - 1e-7) u = 1.0 - 1e-7;
        double g = std::sqrt(2.0) * TMath::ErfInverse(2.0 * u - 1.0);
        return g;
    }

    double pt_corr(bool is_data, const std::string& syst, double pt, double eta, double phi, double deltaEtaSC, float r9,
                   unsigned int run, int seedGain, unsigned long long evt, unsigned int lumi) const {
        // Per EGM POG prescription: scale and smearing corrections are not calibrated below 15 GeV
        if (pt < 15.0) return pt;

        double scEta = eta + deltaEtaSC;
        double r9_val = static_cast<double>(r9);
        double pt_val = pt;

        if (is_data) {
            if (!_dataCompoundScaleCorr) {
                std::cerr << "[ERROR] ElectronSSHelper: _dataCompoundScaleCorr ('Scale') is null for data!" << std::endl;
                return pt;
            }
            double scale = _dataCompoundScaleCorr->evaluate({
                std::string("scale"),
                static_cast<double>(run),
                scEta,
                r9_val,
                pt_val,
                static_cast<double>(seedGain)
            });
            return pt * scale;
        }

        if (!_smearCorr) {
            std::cerr << "[ERROR] ElectronSSHelper: _smearCorr ('SmearAndSyst') is null for MC!" << std::endl;
            return pt;
        }

        double g = get_gaussian_rndm(evt, lumi, phi);

        if (syst.find("elesmearup") != std::string::npos) {
            double sigma_up = _smearCorr->evaluate({std::string("smear_up"), pt_val, r9_val, scEta});
            return pt * (1.0 + g * sigma_up);
        } else if (syst.find("elesmeardown") != std::string::npos) {
            // Note from CMS_egmScaleAndSmearingExample.py: in 2022, "smear_down" can be negative in some cases; clamp at 0.0
            double sigma_down = std::max(0.0, _smearCorr->evaluate({std::string("smear_down"), pt_val, r9_val, scEta}));
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
    std::shared_ptr<const correction::CompoundCorrection> _dataCompoundScaleCorr;
};

#endif

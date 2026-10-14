#ifndef MuonScaReHelper_h
#define MuonScaReHelper_h

#include <cmath>
#include <vector>
#include <string>
#include <memory>
#include <iostream>
#include <algorithm>
#include "TMath.h"
#include "correction.h"

// Crystal Ball distribution helper implementing analytical inversion for random smearing.
// Following official CMS Muon POG ScaReKIT (MuonScaRe.cc).
struct CrystalBall {
    double pi = 3.14159265358979323846;
    double sqrtPiOver2 = std::sqrt(pi / 2.0);
    double sqrt2 = std::sqrt(2.0);
    double m;
    double s;
    double a;
    double n;
    double B;
    double C;
    double D;
    double N;
    double NA;
    double Ns;
    double NC;
    double F;
    double G;
    double k;
    double cdfMa;
    double cdfPa;

    CrystalBall() : m(0.0), s(1.0), a(10.0), n(10.0) {
        init();
    }

    CrystalBall(double mean, double sigma, double alpha, double n_val)
        : m(mean), s(sigma), a(alpha), n(n_val) {
        init();
    }

    void init() {
        double fa = std::fabs(a);
        if (fa < 1e-6) fa = 1e-6;
        if (n <= 1.0) n = 1.0001;
        if (s <= 1e-6) s = 1e-6;

        double ex = std::exp(-fa * fa / 2.0);
        double A  = std::pow(n / fa, n) * ex;
        double C1 = (n / fa / (n - 1.0)) * ex; 
        double D1 = 2.0 * sqrtPiOver2 * std::erf(fa / sqrt2);
        B = n / fa - fa;
        C = (D1 + 2.0 * C1) / C1;   
        D = (D1 + 2.0 * C1) / 2.0;   
        N = 1.0 / s / (D1 + 2.0 * C1); 
        k = 1.0 / (n - 1.0);  
        NA = N * A;       
        Ns = N * s;       
        NC = Ns * C1;     
        F = 1.0 - fa * fa / n; 
        G = s * n / fa;    
        cdfMa = cdf(m - a * s);
        cdfPa = cdf(m + a * s);
    }

    double pdf(double x) const { 
        double d = (x - m) / s;
        if (d < -a) return NA * std::pow(B - d, -n);
        if (d > a)  return NA * std::pow(B + d, -n);
        return N * std::exp(-d * d / 2.0);
    }

    double cdf(double x) const {
        double d = (x - m) / s;
        if (d < -a) return NC / std::pow(F - s * d / G, n - 1.0);
        if (d > a)  return NC * (C - std::pow(F + s * d / G, 1.0 - n));
        return Ns * (D - sqrtPiOver2 * std::erf(-d / sqrt2));
    }

    double invcdf(double u) const {
        if (u <= 0.0) u = 1e-7;
        if (u >= 1.0) u = 1.0 - 1e-7;
        if (u < cdfMa) return m + G * (F - std::pow(NC / u, k));
        if (u > cdfPa) return m - G * (F - std::pow(C - u / NC, -k));
        double arg = (D - u / Ns) / sqrtPiOver2;
        if (arg <= -1.0) arg = -0.99999999999;
        if (arg >= 1.0)  arg = 0.99999999999;
        return m - sqrt2 * s * TMath::ErfInverse(arg);
    }
};

// Evaluator class encapsulating CMS Muon POG ScaReKIT scale & smearing corrections
class MuonScaReHelper {
public:
    MuonScaReHelper(std::shared_ptr<correction::CorrectionSet> cset)
        : _cset(cset) {}

    double get_rndm(double eta, double phi, float nL, unsigned long long evtNumber, unsigned int lumiNumber) const {
        double mean  = _cset->at("cb_params")->evaluate({std::abs(eta), (double)nL, 0});
        double sigma = _cset->at("cb_params")->evaluate({std::abs(eta), (double)nL, 1});
        double n     = _cset->at("cb_params")->evaluate({std::abs(eta), (double)nL, 2});
        double alpha = _cset->at("cb_params")->evaluate({std::abs(eta), (double)nL, 3});

        CrystalBall cb(mean, sigma, alpha, n);
        double rndm = _cset->at("RandomSmearing")->evaluate({static_cast<int>(evtNumber), static_cast<int>(lumiNumber), (double)phi});
        return cb.invcdf(rndm);
    }

    double get_std(double pt, double eta, float nL) const {
        double p0 = _cset->at("poly_params")->evaluate({std::abs(eta), (double)nL, 0});
        double p1 = _cset->at("poly_params")->evaluate({std::abs(eta), (double)nL, 1});
        double p2 = _cset->at("poly_params")->evaluate({std::abs(eta), (double)nL, 2});

        double sigma = p0 + p1 * pt + p2 * pt * pt;
        if (sigma < 0.0) sigma = 0.0;
        return sigma;
    }

    double get_k(double eta, const std::string& var = "nom") const {
        double k_data = _cset->at("k_data")->evaluate({std::abs(eta), var});
        double k_mc   = _cset->at("k_mc")->evaluate({std::abs(eta), var});
        double k = 0.0;
        if (k_mc < k_data) k = std::sqrt(k_data * k_data - k_mc * k_mc);
        return k;
    }

    double pt_scale(bool is_data, double pt, double eta, double phi, int charge, double low_pt_threshold = 26.0) const {
        if (pt < low_pt_threshold) return pt;
        std::string dtmc = is_data ? "data" : "mc";
        double a = _cset->at("a_" + dtmc)->evaluate({(double)eta, (double)phi, "nom"});
        double m = _cset->at("m_" + dtmc)->evaluate({(double)eta, (double)phi, "nom"});
        double denom = m / pt + (double)charge * a;
        if (denom <= 0.0) return pt;
        return 1.0 / denom;
    }

    double pt_scale_var(double pt, double eta, double phi, int charge, const std::string& updn) const {
        double stat_a = _cset->at("a_mc")->evaluate({(double)eta, (double)phi, "stat"});
        double stat_m = _cset->at("m_mc")->evaluate({(double)eta, (double)phi, "stat"});
        double stat_rho = _cset->at("m_mc")->evaluate({(double)eta, (double)phi, "rho_stat"});

        double syst_a = _cset->at("a_mc")->evaluate({(double)eta, (double)phi, "syst"});
        double syst_m = _cset->at("m_mc")->evaluate({(double)eta, (double)phi, "syst"});

        double tot_a = std::sqrt(stat_a * stat_a + syst_a * syst_a);
        double tot_m = std::sqrt(stat_m * stat_m + syst_m * syst_m);

        double unc = pt * pt * std::sqrt(tot_m * tot_m / (pt * pt) + tot_a * tot_a + 2.0 * (double)charge * stat_rho * tot_m / pt * tot_a);

        double pt_var = pt;
        if (updn == "up") {
            pt_var = pt + unc;
        } else if (updn == "dn" || updn == "down") {
            pt_var = pt - unc;
        }
        if (pt_var <= 0.0 || pt_var / pt > 2.0 || pt_var / pt < 0.1) {
            pt_var = pt;
        }
        return pt_var;
    }

    double pt_resol(double pt, double eta, double phi, float nL, unsigned long long evtNumber, unsigned int lumiNumber, double low_pt_threshold = 26.0) const {
        if (pt < low_pt_threshold || pt > 200.0) return pt;
        double rndm = get_rndm(eta, phi, nL, evtNumber, lumiNumber);
        double std_val = get_std(pt, eta, nL);
        double k = get_k(eta, "nom");

        double ptc = pt * (1.0 + k * std_val * rndm);
        if (std::isnan(ptc) || ptc <= 0.0 || ptc / pt > 2.0 || ptc / pt < 0.1) {
            ptc = pt;
        }
        return ptc;
    }

    double pt_resol_var(double pt_woresol, double pt_wresol, double eta, const std::string& updn) const {
        double k = get_k(eta, "nom");
        if (k == 0.0) return pt_wresol;

        double k_stat = _cset->at("k_mc")->evaluate({std::abs(eta), "stat"});
        double k_syst = _cset->at("k_mc")->evaluate({std::abs(eta), "syst"});
        double k_unc = std::sqrt(k_stat * k_stat + k_syst * k_syst);

        double std_x_rndm = (pt_wresol / pt_woresol - 1.0) / k;
        double pt_var = pt_wresol;

        if (updn == "up") {
            pt_var = pt_woresol * (1.0 + (k + k_unc) * std_x_rndm);
        } else if (updn == "dn" || updn == "down") {
            pt_var = pt_woresol * (1.0 + (k - k_unc) * std_x_rndm);
        }
        if (std::isnan(pt_var) || pt_var <= 0.0 || pt_var / pt_woresol > 2.0 || pt_var / pt_woresol < 0.1) {
            pt_var = pt_woresol;
        }
        return pt_var;
    }

private:
    std::shared_ptr<correction::CorrectionSet> _cset;
};

#endif

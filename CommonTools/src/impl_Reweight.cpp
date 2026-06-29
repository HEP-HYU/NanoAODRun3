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



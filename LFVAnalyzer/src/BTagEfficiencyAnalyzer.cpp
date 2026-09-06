/*
 * BTagEfficiencyAnalyzer.cpp
 *
 *  BTV Method 1a: C++ RDataFrame-based b-tagging efficiency map extractor.
 *  Inherits from TopLFVAnalyzer to reuse full object selections (muon, elec, tau, jet overlap removal)
 *  while keeping the pre-tag baseline phase space (Njets >= 3, without requiring b-tags).
 */

#include "BTagEfficiencyAnalyzer.h"
#include "utility.h"
#include <iostream>
#include <vector>

BTagEfficiencyAnalyzer::BTagEfficiencyAnalyzer(TTree *t, std::string outfilename, std::string year, std::string ch, std::string syst, std::string jsonfname, bool applytauFF, std::string globaltag, int nthreads)
: TopLFVAnalyzer(t, outfilename, year, ch, syst, jsonfname, applytauFF, globaltag, nthreads)
{
    _chName = _isMuonCh ? "muon" : "electron";
    std::cout << "=========================================================" << std::endl;
    std::cout << "[BTagEfficiencyAnalyzer] Initialized" << std::endl;
    std::cout << "  Channel : " << _chName << std::endl;
    std::cout << "  Year    : " << _year << std::endl;
    std::cout << "  Systematic : " << _syst << std::endl;

    // Determine btag cut and column name based on era
    // 2022: 0.2450, 2022EE: 0.2605, 2023: 0.1919, 2023BPix: 0.1917, 2024: 0.1272
    if (_year == "2022" || _year == "2022postEE" || _year == "2022EE") {
        _btagCol = "Jet_btagPNetB";
        _btagcut = (_year == "2022") ? 0.2450f : 0.2605f;
    } else if (_year == "2023" || _year == "2023BPix") {
        _btagCol = "Jet_btagPNetB";
        _btagcut = (_year == "2023") ? 0.1919f : 0.1917f;
    } else if (_year.find("2024") != std::string::npos) {
        _btagCol = "Jet_btagUParTAK4B";
        _btagcut = 0.1272f;
    } else {
        _btagCol = "Jet_btagPNetB";
        _btagcut = 0.2450f;
    }
    std::cout << "  B-tag Col  : " << _btagCol << std::endl;
    std::cout << "  Medium WP  : " << _btagcut << std::endl;
    std::cout << "=========================================================" << std::endl;
}

void BTagEfficiencyAnalyzer::defineObjectSelection(std::vector<std::string> jes_var) {
    // Reuse full TopLFVAnalyzer object selection (muons, electrons, taus, jets, overlap removal)
    TopLFVAnalyzer::defineObjectSelection(jes_var);
}

void BTagEfficiencyAnalyzer::defineCuts() {
    // Pre-tag baseline cuts: identical to Signal Region but WITHOUT ncleanbjetspass cut!
    std::string lepCut = _isMuonCh ?
        "nmuonpass == 1 && nvetoelepass == 0 && nvetomuons == 0 && PV_npvsGood > 0" :
        "nelepass == 1 && nvetoelepass == 0 && nvetomuons == 0 && PV_npvsGood > 0";
    std::string chargeCut = "leptau_charge < 0";

    // S1: 1 Good Lepton + Lepton Veto + Good PV
    addCuts(lepCut, "0");
    // S2: Exactly 1 Clean Tau
    addCuts("ncleantaupass == 1", "00");
    // S3: Opposite-Sign Charge (q_lep * q_tau < 0)
    addCuts(chargeCut, "000");
    // S4: Clean Jet Multiplicity >= 3 (Pre-tag Baseline)
    addCuts("ncleanjetspass >= 3", "0000");
}

void BTagEfficiencyAnalyzer::defineMoreVars() {
    // 1. Charge product for OS selection
    if (_isMuonCh) {
        addVar({"leptau_charge", "Muon_charge[0] * Tau_charge[0]", ""});
    } else {
        addVar({"leptau_charge", "Electron_charge[0] * Tau_charge[0]", ""});
    }

    // 2. Absolute eta for clean jets (Jet_pt, Jet_eta, Jet_hadronFlavour are already cleaned in selectJets)
    addVar({"Jet_abseta", "abs(Jet_eta)", ""});

    // 3. Flavor masks (hadronFlavour: 5 = b, 4 = c, 0 = light/gluon)
    addVar({"b_mask", "Jet_hadronFlavour == 5", ""});
    addVar({"c_mask", "Jet_hadronFlavour == 4", ""});
    addVar({"l_mask", "Jet_hadronFlavour != 5 && Jet_hadronFlavour != 4", ""});

    // 4. Tagged mask (Medium WP)
    std::string tagExpr = _btagCol + " >= " + std::to_string(_btagcut) + "f";
    addVar({"tag_mask", tagExpr, ""});

    // 5. Kinematics per flavor (Denominator)
    addVar({"b_pt", "Jet_pt[b_mask]", ""});
    addVar({"b_abseta", "Jet_abseta[b_mask]", ""});

    addVar({"c_pt", "Jet_pt[c_mask]", ""});
    addVar({"c_abseta", "Jet_abseta[c_mask]", ""});

    addVar({"l_pt", "Jet_pt[l_mask]", ""});
    addVar({"l_abseta", "Jet_abseta[l_mask]", ""});

    // 6. Kinematics per flavor passing Medium WP (Numerator)
    addVar({"b_tagged_pt", "Jet_pt[b_mask && tag_mask]", ""});
    addVar({"b_tagged_abseta", "Jet_abseta[b_mask && tag_mask]", ""});

    addVar({"c_tagged_pt", "Jet_pt[c_mask && tag_mask]", ""});
    addVar({"c_tagged_abseta", "Jet_abseta[c_mask && tag_mask]", ""});

    addVar({"l_tagged_pt", "Jet_pt[l_mask && tag_mask]", ""});
    addVar({"l_tagged_abseta", "Jet_abseta[l_mask && tag_mask]", ""});

    // 7. Event weights (excluding b-tag SF, including gen, pu, top pt, lepton, tau ID/iso SFs)
    if (_isData) {
        addVar({"eventWeight", "1.0", ""});
    } else {
        std::string lepWeight = _isMuonCh ?
            "muonWeightId[0] * muonWeightIso[0] * muonWeightTrg[0]" :
            "elecWeightReco[0] * elecWeightId[0] * elecWeightTrg[0]";

        std::string weightExpr = "unitGenWeight * TopPtWeight[0] * puWeight[0] * " + lepWeight +
                                 " * tauWeightIdVsJet[0][0] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]";
        addVar({"eventWeight", weightExpr, ""});
    }

    // 8. Vector weights per jet matching collection sizes
    addVar({"b_weight", "ROOT::RVecF(b_pt.size(), eventWeight)", ""});
    addVar({"b_tagged_weight", "ROOT::RVecF(b_tagged_pt.size(), eventWeight)", ""});
    addVar({"c_weight", "ROOT::RVecF(c_pt.size(), eventWeight)", ""});
    addVar({"c_tagged_weight", "ROOT::RVecF(c_tagged_pt.size(), eventWeight)", ""});
    addVar({"l_weight", "ROOT::RVecF(l_pt.size(), eventWeight)", ""});
    addVar({"l_tagged_weight", "ROOT::RVecF(l_tagged_pt.size(), eventWeight)", ""});
}

void BTagEfficiencyAnalyzer::bookHists() {
    maxstep = "0000";

    // POG Recommended binning:
    // pT:  [40, 50, 70, 100, 140, 200, 300, 600, 1000]
    // eta: [0.0, 2.5]
    const std::vector<double> pt_bins = {40.0, 50.0, 70.0, 100.0, 140.0, 200.0, 300.0, 600.0, 1000.0};
    const std::vector<double> eta_bins = {0.0, 2.5};

    TH2DModel h2_denom_b_model("h2_denom_b", ";Jet p_{T} [GeV];Jet |#eta|", pt_bins.size() - 1, pt_bins.data(), eta_bins.size() - 1, eta_bins.data());
    TH2DModel h2_num_b_model("h2_num_b", ";Jet p_{T} [GeV];Jet |#eta|", pt_bins.size() - 1, pt_bins.data(), eta_bins.size() - 1, eta_bins.data());

    TH2DModel h2_denom_c_model("h2_denom_c", ";Jet p_{T} [GeV];Jet |#eta|", pt_bins.size() - 1, pt_bins.data(), eta_bins.size() - 1, eta_bins.data());
    TH2DModel h2_num_c_model("h2_num_c", ";Jet p_{T} [GeV];Jet |#eta|", pt_bins.size() - 1, pt_bins.data(), eta_bins.size() - 1, eta_bins.data());

    TH2DModel h2_denom_l_model("h2_denom_light", ";Jet p_{T} [GeV];Jet |#eta|", pt_bins.size() - 1, pt_bins.data(), eta_bins.size() - 1, eta_bins.data());
    TH2DModel h2_num_l_model("h2_num_light", ";Jet p_{T} [GeV];Jet |#eta|", pt_bins.size() - 1, pt_bins.data(), eta_bins.size() - 1, eta_bins.data());

    // Book 2D Histograms for Denominator & Numerator
    add2DHist(h2_denom_b_model, "b_pt", "b_abseta", "b_weight", "", "0000", maxstep);
    add2DHist(h2_num_b_model, "b_tagged_pt", "b_tagged_abseta", "b_tagged_weight", "", "0000", maxstep);

    add2DHist(h2_denom_c_model, "c_pt", "c_abseta", "c_weight", "", "0000", maxstep);
    add2DHist(h2_num_c_model, "c_tagged_pt", "c_tagged_abseta", "c_tagged_weight", "", "0000", maxstep);

    add2DHist(h2_denom_l_model, "l_pt", "l_abseta", "l_weight", "", "0000", maxstep);
    add2DHist(h2_num_l_model, "l_tagged_pt", "l_tagged_abseta", "l_tagged_weight", "", "0000", maxstep);

    // Book 1D control histograms for event counting & verification
    add1DHist({"h_nevents", ";Pre-tag Baseline Events;Events", 2, -0.5, 1.5}, "one", "eventWeight", "", "0000", maxstep);
    add1DHist({"h_ncleanjetspass", ";Clean Jet multiplicity;Events", 10, 0.0, 10.0}, "ncleanjetspass", "eventWeight", "", "0000", maxstep);
}

/*
 * TopLFVAnalyzer.cpp
 *
 *  Created on: April 9, 2020
 *      Author: Tae Jeong Kim
 *      Refactored: Jiwon Park
 */

#include "TopLFVAnalyzer.h"
#include "utility.h"
#include <algorithm>
#include <fstream>

TopLFVAnalyzer::TopLFVAnalyzer(TTree *t, std::string outfilename, std::string year, std::string ch, std::string syst, std::string jsonfname, bool applytauFF, string globaltag, int nthreads)
:NanoAODAnalyzerrdframe(t, outfilename, year, ch, syst, jsonfname, globaltag, nthreads), _applytauFF(applytauFF)
{
    cout << "<< Start Process NanoAOD >>" << endl;

    if(syst.find("jes") != std::string::npos or syst.find("jer") != std::string::npos or syst.find("metUnclust") != std::string::npos or
            syst.find("tes") != std::string::npos or syst.find("hdamp") != std::string::npos or syst.find("tune") != std::string::npos or
            syst.find("muonhighscale") != std::string::npos or syst.find("btag") != std::string::npos) {
        ext_syst = true;
    }
    cout << "syst: " << _syst << endl;

    tauYear = _year;

    if (_outfilename.find("Tau-LFV-") != std::string::npos) {
        _isSignal = true;
        cout << "Input file is LFV signal" <<endl;
    } else {
        _isSignal = false;
    }
    if (_ch.find("muon") != std::string::npos){
        cout << "muon channel" << endl;
        _isMuonCh = true;
    } else{
        cout << "electron channel" << endl;
        _isMuonCh = false;
    }
}

void TopLFVAnalyzer::defineBTagNormalization(){
    if (_isData) return;
    std::string lepWeight = "";
    std::string preBTagCut = "";
    if (_isMuonCh) {
        lepWeight = "muonWeightId[0] * muonWeightIso[0] * muonWeightTrg[0]";
        preBTagCut = "nmuonpass == 1 && nvetoelepass == 0 && nvetomuons == 0 && "
                     "ncleantaupass == 1 && Muon_charge[0] * Tau_charge[0] < 0 && "
                     "ncleanjetspass >= 3 && PV_npvsGood > 0";
    } else {
        lepWeight = "elecWeightReco[0] * elecWeightId[0] * elecWeightTrg[0]";
        preBTagCut = "nelepass == 1 && nvetoelepass == 0 && nvetomuons == 0 && "
                     "ncleantaupass == 1 && Electron_charge[0] * Tau_charge[0] < 0 && "
                     "ncleanjetspass >= 3 && PV_npvsGood > 0";
    }

    const std::string baseWeight = "double(unitGenWeight * TopPtWeight[0] * puWeight[0] * " +
                                   lepWeight +
                                   " * tauWeightIdVsJet[0][0] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0])";

    auto df_tmp = _rlm.Define("_btagNormBaseWeight", baseWeight).Filter(preBTagCut);
    auto njet_col = df_tmp.Take<int>("ncleanjetspass");
 
    auto w_nom_col = df_tmp.Take<double>("_btagNormBaseWeight");
    auto btag_col = df_tmp.Take<floats>("btagWeight");


    const size_t minBTagVars = 17;
    std::map<int, float> sum_nominal;
    std::map<int, floats> sum_btag;
    for (size_t i = 0; i < njet_col->size(); ++i) {
        const int njet = (*njet_col)[i];
        const float w_nom = (*w_nom_col)[i];
        const floats &btag_weights = (*btag_col)[i];
        const size_t nBTagVars = std::max(minBTagVars, size_t(btag_weights.size()));

        sum_nominal[njet] += w_nom;
        if (sum_btag[njet].size() < nBTagVars){
            sum_btag[njet].resize(nBTagVars, 0.0);
        }
        for (size_t ivar = 0; ivar < nBTagVars; ++ivar){
            const float btag_weight = ivar < btag_weights.size() ? btag_weights[ivar] : 1.0;
            sum_btag[njet][ivar] += w_nom * btag_weight;
        }
    }

    std::map<int, floats> norm_factors;
    for (auto &kv : sum_nominal) {
        const int njet = kv.first;
        const float numerator = sum_nominal[njet];
        floats norms;
        norms.reserve(sum_btag[njet].size());

        for(size_t ivar=0; ivar<sum_btag[njet].size(); ++ivar){
            const float denom = sum_btag[njet][ivar];
            const float norm = denom != 0.0 ? numerator / denom : 1.0;
            norms.emplace_back(norm);
        }

        norm_factors[njet] = norms;
        std::cout << "BTag norm factor njet = " << njet;
        for (size_t ivar=0; ivar<norms.size(); ++ivar){
            std::cout << " [" << ivar << "]=" << norms[ivar];
        }
        std::cout << std::endl;
    }

    auto getBTagNorm = [norm_factors, minBTagVars](int njet)->floats {
        auto it = norm_factors.find(njet);
        if (it != norm_factors.end()) return it->second;

        floats out(minBTagVars, 1.0f);
        return out;
    };

    auto applyBTagNorm = [minBTagVars](floats &weights, floats &norms)->floats {
        floats out;

        const size_t nBTagVars = std::max({size_t(weights.size()), size_t(norms.size()), minBTagVars});
        out.reserve(nBTagVars);
        
        for (size_t ivar=0; ivar<nBTagVars; ++ivar){
            const float weight = ivar < weights.size() ? weights[ivar] : 1.0;
            const float norm = ivar < norms.size() ? norms[ivar] : 1.0;
            out.emplace_back(weight * norm);
        }

        return out;
    };

    defineVar("btagNorm", getBTagNorm, {"ncleanjetspass"});
    defineVar("btagWeightNorm", applyBTagNorm, {"btagWeight", "btagNorm"});
}

void TopLFVAnalyzer::defineObjectSelection(std::vector<std::string> jes_var){
    // ── Load analysis configuration (static cache, read once per process) ────
    // All physics thresholds are defined in data/config/analysis_config.json.
    // Values were verified against previous hardcoded literals (Python PASS).
    static bool cfg_loaded = false;
    static std::string s_vetomuon, s_vetoelec, s_jetcut, s_taucut;
    static std::string s_mu_vsjet, s_mu_vsmu, s_mu_vse;
    static std::string s_el_vsjet, s_el_vsmu, s_el_vse;
    static std::map<std::string, float> s_btagcut_map;
    static std::string s_btagMap, s_btagMapLight;

    if (!cfg_loaded) {
        const std::string cfgPath = "data/config/analysis_config.json";
        std::ifstream fin(cfgPath);
        if (!fin.good()) {
            std::cerr << "ERROR defineObjectSelection: cannot open " << cfgPath
                      << " — aborting" << std::endl;
            std::abort();
        }
        Json::Value root;
        fin >> root;

        const auto &sel  = root["object_selection"];
        const auto &dtau = root["deeptau_wp"];
        const auto &bt   = root["btag"];

        s_vetomuon     = sel["muon_veto_cut"].asString();
        s_vetoelec     = sel["electron_veto_cut"].asString();
        s_jetcut       = sel["jet_cut"].asString();
        s_taucut       = sel["tau_cut"].asString();
        s_mu_vsjet     = dtau["muon_ch"]["vsjet"].asString();
        s_mu_vsmu      = dtau["muon_ch"]["vsmu"].asString();
        s_mu_vse       = dtau["muon_ch"]["vse"].asString();
        s_el_vsjet     = dtau["electron_ch"]["vsjet"].asString();
        s_el_vsmu      = dtau["electron_ch"]["vsmu"].asString();
        s_el_vse       = dtau["electron_ch"]["vse"].asString();
        if (bt["wp_medium"].isObject()) {
            for (auto const& id : bt["wp_medium"].getMemberNames()) {
                s_btagcut_map[id] = bt["wp_medium"][id].asFloat();
            }
        }
        s_btagMap      = bt["tagger_run3"].asString();
        s_btagMapLight = bt["tagger_light_run3"].asString();
        cfg_loaded = true;
    }

    std::string cut  = "onlyveto";

    std::string btagYear = "";
    std::string btagEra  = "";
    std::string btagMap      = s_btagMap;
    std::string btagMapLight = s_btagMapLight;
    if (_isRun22) {
        tauYear  = "2022_preEE";
        btagYear = "2022_Summer22";
        btagEra  = "2022_Summer22";
    } else if (_isRun22EE) {
        tauYear  = "2022_postEE";
        btagYear = "2022_Summer22EE";
        btagEra  = "2022_Summer22EE";
    } else if (_isRun23) {
        tauYear  = "2023_preBPix";
        btagYear = "2023_Summer23";
        btagEra  = "2023_Summer23";
    } else if (_isRun23BPix) {
        tauYear  = "2023_postBPix";
        btagYear = "2023_Summer23BPix";
        btagEra  = "2023_Summer23BPix";
    } else if (_isRun24){
        tauYear  = "2024";
        // TODO: Run3Summer24 UParTAK4 shape SF not yet available.
        // Temporarily use 2023BPix particleNet_shape correction.
        btagYear = "2023_Summer23BPix";
        btagEra  = "2024_Summer24";
        // btagMap stays "particleNet_shape" (default from config)
    }
    // NanoAODv12 btag JSONs do not have "UParTAK4_light".
    // Use "particleNet_shape" for all eras (exists in all v12 JSONs and in the
    // 2023BPix JSON used as the 2024 temporary fallback).
    btagMapLight = "particleNet_shape";
    float btagcut = s_btagcut_map.count(btagEra) ? s_btagcut_map[btagEra] : 0.1272f;

    // MET unclustered energy variation:
    // NanoAOD pre-computes PuppiMET_ptUnclusteredUp/Down and PuppiMET_phiUnclusteredUp/Down.
    // Redefine the nominal branches BEFORE any downstream variable (lepMET_mt, st_met, etc.)
    // so that all derived kinematics automatically pick up the shifted MET.
    if (_syst.find("metUnclustup") != std::string::npos) {
        _rlm = _rlm.Redefine("PuppiMET_pt",  "PuppiMET_ptUnclusteredUp")
                   .Redefine("PuppiMET_phi", "PuppiMET_phiUnclusteredUp");
        cout << "[metUnclust] Swapped PuppiMET_pt/phi -> Unclustered UP" << endl;
    } else if (_syst.find("metUnclustdown") != std::string::npos) {
        _rlm = _rlm.Redefine("PuppiMET_pt",  "PuppiMET_ptUnclusteredDown")
                   .Redefine("PuppiMET_phi", "PuppiMET_phiUnclusteredDown");
        cout << "[metUnclust] Swapped PuppiMET_pt/phi -> Unclustered DOWN" << endl;
    }

    if (_isMuonCh){
        selectElectrons(cut, s_vetoelec);
        selectTaus(s_taucut, tauYear, s_mu_vsjet, s_mu_vsmu, s_mu_vse); //VSjet VVTight, VSmu Tight, VSe VVLoose
    } else {
        selectMuons(cut, s_vetomuon);
        selectTaus(s_taucut, tauYear, s_el_vsjet, s_el_vsmu, s_el_vse); //VSjet VVTight, VSmu VLoose, VSe VTight
    }
    selectJets(jes_var, jes_var_flav, s_jetcut);
    if (!_isData){
        topPtReweight();
        applyBSFs(jes_var, btagYear, btagMap, btagMapLight, btagcut);
    }
}


// Define your cuts here
void TopLFVAnalyzer::defineCuts() {
    if (_isMuonCh) {
        cout << "muon channel" << endl;
        addCuts("nmuonpass == 1 && nvetoelepass == 0 && nvetomuons == 0 && PV_npvsGood > 0", "0");
    } else {
        cout << "electron channel" << endl;
        addCuts("nelepass == 1 && nvetoelepass == 0 && nvetomuons == 0 && PV_npvsGood > 0", "0");
    }

    addCuts("ncleantaupass == 1", "00");

    // genuine taus
    //if (_isSignal or _syst == "data") {
    //    addCuts("ncleantaupass == 1", "00");
    //} else {
    //    addCuts("ncleantaupass == 1 && Tau_pt_gen.size()>0", "00");
    //}
    // fake taus
    //addCuts("ncleantaupass == 1 && Tau_pt_gen.size()==0", "00");


    addCuts("leptau_charge < 0", "000");
    addCuts("ncleanjetspass >= 3", "0000");
    addCuts("ncleanbjetspass == 1", "00000");
}

void TopLFVAnalyzer::defineMoreVars() {
    defineKinematicVars();     // 4-vectors, ΔR, MT, ST^MET, top reco, object branches
    defineWeightVars();        // nominal & systematic eventWeight definitions
    storeOutputBranches();     // addVartoStore calls
}

void TopLFVAnalyzer::defineKinematicVars() {
    defineVar("lepvec", ::select_leadingvec, {"lep4vecs"});
    defineVar("lepMET_mt", ::calculate_MT, {"lep4vecs","PuppiMET_pt","PuppiMET_phi"});
    
    defineVar("tauvec", ::select_leadingvec, {"cleantau4vecs"});
    defineVar("leptau_mass", ::calculate_invMass, {"lepvec","tauvec"});
    defineVar("leptau_dEta", ::calculate_deltaEta, {"lepvec","tauvec"});
    defineVar("leptau_dPhi", ::calculate_deltaPhi, {"lepvec","tauvec"});
    defineVar("leptau_dR", ::calculate_deltaR, {"lepvec","tauvec"});

    // Temporary solution to blind low mutau mass
    if (_syst == "data") {
        addVar({"leptau_mass_blind", "(leptau_mass < 100) ? -1.0 : leptau_mass"});
    } else {
        addVar({"leptau_mass_blind", "leptau_mass"});
    }
    
    //Adding unc for muon SF to top phase space: done in plotit and combine tool
    //defineVar("muonWeightId", ::addMuonUnc, {"muonWeightId"});
    //defineVar("muonWeightIso", ::addMuonUnc, {"muonWeightIso"});
    //defineVar("muonWeightTrg", ::addMuonUnc, {"muonWeightTrg"});

    // isFakeTau: true when the selected tau is not matched to a generator-level hadronic tau.
    // In data _isData guard below keeps it safe; in MC Tau_pt_gen is filled by selectTaus.
    if (!_isData) {
        addVar({"isFakeTau", "!(Tau_pt_gen.size()>0)"});
    } else {
        addVar({"isFakeTau", "false"});
    }

    // Tau fake factor branches — always define (stub = 1.0) so addVartoStore("tauFF.*")
    // doesn't fail on the no-FF code path.  The real functors are wired below only when
    // _applytauFF is true.
    if (!_applytauFF || _isData || _isSignal) {
        addVar({"tauFF",        "1.0", ""});
        addVar({"tauFFstatup",  "1.0", ""});
        addVar({"tauFFstatdown","1.0", ""});
        addVar({"tauFFsystup",  "1.0", ""});
        addVar({"tauFFsystdown","1.0", ""});
    } else {
        auto tauFF_nom      = tauFFfunctor(tauYear, "nom",  0);
        auto tauFF_statup   = tauFFfunctor(tauYear, "stat", 1);
        auto tauFF_statdown = tauFFfunctor(tauYear, "stat",-1);
        auto tauFF_systup   = tauFFfunctor(tauYear, "syst", 1);
        auto tauFF_systdown = tauFFfunctor(tauYear, "syst",-1);
        defineVar("tauFF",        tauFF_nom,      {"Tau_pt", "Tau_pt_gen", "Tau_decayMode"});
        defineVar("tauFFstatup",  tauFF_statup,   {"Tau_pt", "Tau_pt_gen", "Tau_decayMode"});
        defineVar("tauFFstatdown",tauFF_statdown, {"Tau_pt", "Tau_pt_gen", "Tau_decayMode"});
        defineVar("tauFFsystup",  tauFF_systup,   {"Tau_pt", "Tau_pt_gen", "Tau_decayMode"});
        defineVar("tauFFsystdown",tauFF_systdown, {"Tau_pt", "Tau_pt_gen", "Tau_decayMode"});
    }

    // unitGenWeightFF: generator weight including the tau fake factor.
    // Note: UFO_reweight was removed — no Run3 UFO reweighting prescription.
    addVar({"unitGenWeightFF", "unitGenWeight * tauFF", ""});

    // Safe-index helper: "(col.size()>n) ? float(col[n]) : fallback"
    // Prevents UB / heap corruption when the vector is empty at early cut levels.
    auto si  = [](const std::string& c, int n, const std::string& fb="-1.f") -> std::string {
        return "(" + c + ".size()>" + std::to_string(n) + ") ? float(" + c + "[" + std::to_string(n) + "]) : " + fb;
    };
    auto sii = [](const std::string& c, int n, const std::string& fb="-1") -> std::string {
        return "(" + c + ".size()>" + std::to_string(n) + ") ? int(" + c + "[" + std::to_string(n) + "]) : " + fb;
    };

    if (_isMuonCh) {
        addVar({"Muon1_pt",     si("Muon_pt",     0),        ""});
        addVar({"Muon1_eta",    si("Muon_eta",    0, "-99.f"), ""});
        addVar({"Muon1_mass",   si("Muon_mass",   0),        ""});
        addVar({"Muon1_charge", sii("Muon_charge", 0, "0"),   ""});
    } else {
        addVar({"Electron1_pt",     si("Electron_pt",     0),        ""});
        addVar({"Electron1_eta",    si("Electron_eta",    0, "-99.f"), ""});
        addVar({"Electron1_mass",   si("Electron_mass",   0),        ""});
        addVar({"Electron1_charge", sii("Electron_charge", 0, "0"),   ""});
    }

    addVar({"Tau1_pt",        si("Tau_pt",       0),        ""});
    addVar({"Tau1_eta",       si("Tau_eta",      0, "-99.f"), ""});
    addVar({"Tau1_mass",      si("Tau_mass",     0),        ""});
    addVar({"Tau1_charge",    sii("Tau_charge",  0, "0"),   ""});
    addVar({"Tau1_decayMode", sii("Tau_decayMode",0, "-1"), ""});

    if (_isMuonCh) {
        addVar({"leptau_charge", "Muon1_charge * Tau1_charge", ""});
    } else {
        addVar({"leptau_charge", "Electron1_charge * Tau1_charge", ""});
    }

    addVar({"Jet1_pt",   si("Jet_pt",   0),        ""});
    addVar({"Jet1_eta",  si("Jet_eta",  0, "-99.f"), ""});
    addVar({"Jet1_mass", si("Jet_mass", 0),        ""});
    addVar({"Jet1_btagPNetB", si("Jet_btagPNetB", 0), ""});

    // NanoAODv12 (22/23 eras) does not have UParTAK4 tagger branches in the skim output.
    // Use PNet equivalents for pre-2024; keep UParTAK4 for v15_2024.
    const std::string bTagCol    = _isRun24 ? "Jet_btagUParTAK4B"   : "Jet_btagPNetB";
    const std::string bTagColCvB = _isRun24 ? "Jet_btagUParTAK4CvB" : "Jet_btagPNetCvB";
    const std::string bTagColCvL = _isRun24 ? "Jet_btagUParTAK4CvL" : "Jet_btagPNetCvL";

    // Read btag medium WP from config per era
    std::string btagEra = "";
    if (_isRun22) btagEra = "2022_Summer22";
    else if (_isRun22EE) btagEra = "2022_Summer22EE";
    else if (_isRun23) btagEra = "2023_Summer23";
    else if (_isRun23BPix) btagEra = "2023_Summer23BPix";
    else if (_isRun24) btagEra = "2024_Summer24";

    static std::map<std::string, float> s_btagcut_kin_map;
    static bool btagkin_loaded = false;
    if (!btagkin_loaded) {
        const std::string cfgPath = "data/config/analysis_config.json";
        std::ifstream fin(cfgPath);
        if (fin.good()) {
            Json::Value root; fin >> root;
            if (root["btag"]["wp_medium"].isObject()) {
                for (auto const& id : root["btag"]["wp_medium"].getMemberNames()) {
                    s_btagcut_kin_map[id] = root["btag"]["wp_medium"][id].asFloat();
                }
            }
        }
        btagkin_loaded = true;
    }
    float _btagcutKin = s_btagcut_kin_map.count(btagEra) ? s_btagcut_kin_map[btagEra] : 0.1272f;
    const std::string bWP = std::to_string(_btagcutKin);

    // Safe boolean index helper: "(col.size()>n) ? (col[n] op thr) : false"
    auto sib = [](const std::string& c, int n, const std::string& op, const std::string& thr) -> std::string {
        return "(" + c + ".size()>" + std::to_string(n) + ") ? (" + c + "[" + std::to_string(n) + "] " + op + " " + thr + ") : false";
    };

    addVar({"Jet1_btagUParTAK4B",   si(bTagCol,    0), ""});
    addVar({"Jet1_btagUParTAK4CvB", si(bTagColCvB, 0), ""});
    addVar({"Jet1_btagUParTAK4CvL", si(bTagColCvL, 0), ""});
    addVar({"Jet1_btag", sib(bTagCol, 0, ">", bWP), ""});

    addVar({"Jet2_pt",   si("Jet_pt",   1),        ""});
    addVar({"Jet2_eta",  si("Jet_eta",  1, "-99.f"), ""});
    addVar({"Jet2_mass", si("Jet_mass", 1),        ""});
    addVar({"Jet2_btagPNetB",        si("Jet_btagPNetB", 1), ""});
    addVar({"Jet2_btagUParTAK4B",   si(bTagCol,    1), ""});
    addVar({"Jet2_btagUParTAK4CvB", si(bTagColCvB, 1), ""});
    addVar({"Jet2_btagUParTAK4CvL", si(bTagColCvL, 1), ""});
    addVar({"Jet2_btag", sib(bTagCol, 1, ">", bWP), ""});

    addVar({"Jet3_pt",   si("Jet_pt",   2),        ""});
    addVar({"Jet3_eta",  si("Jet_eta",  2, "-99.f"), ""});
    addVar({"Jet3_mass", si("Jet_mass", 2),        ""});
    addVar({"Jet3_btagPNetB",        si("Jet_btagPNetB", 2), ""});
    addVar({"Jet3_btagUParTAK4B",   si(bTagCol,    2), ""});
    addVar({"Jet3_btagUParTAK4CvB", si(bTagColCvB, 2), ""});
    addVar({"Jet3_btagUParTAK4CvL", si(bTagColCvL, 2), ""});
    addVar({"Jet3_btag", sib(bTagCol, 2, ">", bWP), ""});

    addVar({"bJet1_pt",   si("bJet_pt",   0), ""});
    addVar({"bJet1_eta",  si("bJet_eta",  0, "-99.f"), ""});
    addVar({"bJet1_mass", si("bJet_mass", 0), ""});


    //if (_syst == "data") {
    //    defineVar("tauWeightIdVsJetFIX", ::fixtausf, {"tauWeightIdVsJet", "Tau_pt"});
    //}

    // Reconstruction
    defineVar("top_reco_whad", ::top_reconstruction_STLFV, {"cleanjet4vecs","cleanbjet4vecs","lep4vecs","cleantau4vecs"});
    addVar({"chi2", "top_reco_whad[0]",""});
    addVar({"chi2_SMW_mass", "top_reco_whad[1]",""});
    addVar({"chi2_SMTop_mass", "top_reco_whad[2]",""});
    addVar({"chi2_wjet1_idx", "top_reco_whad[3]",""});
    addVar({"chi2_wjet2_idx", "top_reco_whad[4]",""});
    addVar({"chi2_SMW", "top_reco_whad[5]",""});
    addVar({"chi2_SMTop", "top_reco_whad[6]",""});
    addVar({"chi2_SMTop_pt", "top_reco_whad[7]",""});

    defineVar("top_reco_prod", ::top_reco_products_STLFV, {"cleanjet4vecs","lep4vecs","cleantau4vecs","top_reco_whad"});
    addVar({"chi2_wqq_dEta","top_reco_prod[0]",""});
    addVar({"chi2_wqq_dPhi","top_reco_prod[1]",""});
    addVar({"chi2_wqq_dR","top_reco_prod[2]",""});

    // S_T^MET: in EXO-19-016, defined by pt1+pt2+ptj1+met where pt1,2 come from taus
    if (_isMuonCh){
        defineVar("st_met", ::st_met, {"Muon_pt", "Tau_pt", "Jet_pt", "PuppiMET_pt"});
    } else {
        defineVar("st_met", ::st_met, {"Electron_pt", "Tau_pt", "Jet_pt", "PuppiMET_pt"});
    }

}

void TopLFVAnalyzer::defineWeightVars() {
    // EventWeights
    // Calculate product of weights and store for systematic study
    // External systs: JES (+btag) (14), JER(2), TES(2), hdamp(2), TuneCP5(2)
    // Weights systs: genWeight(1), PU(2), btag(16), muon Id(2)/Iso(2)/Trg(2), tauId(18+2*2), ME scale(6), PS scale(ISR 2, FSR 2), PDF(102, or 2)
    // Not implemented: EEprefire, top pt reweighting,
    // eventWeight_xx : xxweight
    // eventWeight__xx: xx unc.



    addVar({"eventWeight_no", "1.0"});
    if (_syst == "data") {
        addVar({"eventWeight", "1.0"});
        addVar({"eventWeight_notau", "1.0"});
        addVar({"eventWeight_notoppt", "1.0"});
        addVar({"eventWeight_notau_nobtag", "1.0"}); //didn't want to duplicate entry...
        addVar({"eventWeight_nobtag", "1.0"});
    } else {
        defineBTagNormalization();
        addVar({"eventWeight_genpu", "unitGenWeight * TopPtWeight[0] * puWeight[0]"});
        if (_isMuonCh){
            addVar({"eventWeight_lep", "muonWeightId[0] * muonWeightIso[0] * muonWeightTrg[0]"});
        } else {
            addVar({"eventWeight_lep", "elecWeightReco[0] * elecWeightId[0] * elecWeightTrg[0]"});
        }
        addVar({"eventWeight_tau", "tauWeightIdVsJet[0][0] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
        addVar({"eventWeight_genpumu", "eventWeight_genpu * eventWeight_lep"});
        addVar({"eventWeight_notau_nobtag", "eventWeight_genpumu"}); //didn't want to duplicate entry...
        addVar({"eventWeight_genputau", "eventWeight_genpu * eventWeight_tau"});
        addVar({"eventWeight_nobtag", "eventWeight_genpu * eventWeight_lep * eventWeight_tau"});
        addVar({"eventWeight_nopu", "unitGenWeight * TopPtWeight[0] * eventWeight_lep * eventWeight_tau * btagWeightNorm[0]"});
        addVar({"eventWeight_noprefire", "unitGenWeight * TopPtWeight[0] * puWeight[0] * eventWeight_lep * eventWeight_tau * btagWeightNorm[0]"});
        addVar({"eventWeight_notoppt", "unitGenWeight * puWeight[0] * eventWeight_lep * eventWeight_tau * btagWeightNorm[0]"});

        //TODO
        // eventWeight, eventWeight_notau
        // defined here for ALL non-data modes (nosyst, all, theory, ext_syst).
        // The if-else branches below only add ADDITIONAL syst-specific variations.
        addVar({"eventWeight",         "eventWeight_nobtag * btagWeightNorm[0]"});
        addVar({"eventWeight_notau",   "eventWeight_genpumu * btagWeightNorm[0]"});

        if (_syst == "" or _syst == "nosyst" or ext_syst) {
            // nosyst / ext_syst: only nominal weights needed, no additional weight variations.
            // No additional addVar calls needed here.
        } else {

            addVar({"eventWeight__puup", "eventWeight_nopu * puWeight[1]"});
            addVar({"eventWeight__pudown", "eventWeight_nopu * puWeight[2]"});
            addVar({"eventWeight__topptup", "eventWeight_nobtag * btagWeightNorm[0] * TopPtWeight[1] / TopPtWeight[0]"});
            addVar({"eventWeight__topptdown", "eventWeight_nobtag * btagWeightNorm[0] * TopPtWeight[2] / TopPtWeight[0]"});
            // Lepton SF variations — muon channel
            if (_isMuonCh) {
                addVar({"eventWeight__muidup",     "eventWeight_genputau * muonWeightId[1] * muonWeightIso[0] * muonWeightTrg[0] * btagWeightNorm[0]"});
                addVar({"eventWeight__muiddown",   "eventWeight_genputau * muonWeightId[2] * muonWeightIso[0] * muonWeightTrg[0] * btagWeightNorm[0]"});
                addVar({"eventWeight__muisoup",    "eventWeight_genputau * muonWeightId[0] * muonWeightIso[1] * muonWeightTrg[0] * btagWeightNorm[0]"});
                addVar({"eventWeight__muisodown",  "eventWeight_genputau * muonWeightId[0] * muonWeightIso[2] * muonWeightTrg[0] * btagWeightNorm[0]"});
                addVar({"eventWeight__mutrgup",    "eventWeight_genputau * muonWeightId[0] * muonWeightIso[0] * muonWeightTrg[1] * btagWeightNorm[0]"});
                addVar({"eventWeight__mutrgdown",  "eventWeight_genputau * muonWeightId[0] * muonWeightIso[0] * muonWeightTrg[2] * btagWeightNorm[0]"});
            } else {
                // Electron SF variations — electron channel
                addVar({"eventWeight__elecrecoup",    "eventWeight_genputau * elecWeightReco[1] * elecWeightId[0] * elecWeightTrg[0] * btagWeightNorm[0]"});
                addVar({"eventWeight__elecrecodown",  "eventWeight_genputau * elecWeightReco[2] * elecWeightId[0] * elecWeightTrg[0] * btagWeightNorm[0]"});
                addVar({"eventWeight__elecidup",      "eventWeight_genputau * elecWeightReco[0] * elecWeightId[1] * elecWeightTrg[0] * btagWeightNorm[0]"});
                addVar({"eventWeight__eleciddown",    "eventWeight_genputau * elecWeightReco[0] * elecWeightId[2] * elecWeightTrg[0] * btagWeightNorm[0]"});
                addVar({"eventWeight__electrgup",     "eventWeight_genputau * elecWeightReco[0] * elecWeightId[0] * elecWeightTrg[1] * btagWeightNorm[0]"});
                addVar({"eventWeight__electrgdown",   "eventWeight_genputau * elecWeightReco[0] * elecWeightId[0] * elecWeightTrg[2] * btagWeightNorm[0]"});
            }
            // muonHighPtAddUncUp/Dn removed: no official Run3 high-pT muon additional uncertainty prescription.

            addVar({"eventWeight__tauidjetup", "eventWeight_notau * tauWeightIdVsJet[0][1] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauidjetdown", "eventWeight_notau * tauWeightIdVsJet[0][2] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauidjetStat1dm0"+tauYear+"up", "eventWeight_notau * tauWeightIdVsJet[0][3] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauidjetStat1dm0"+tauYear+"down", "eventWeight_notau * tauWeightIdVsJet[0][4] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauidjetStat1dm1"+tauYear+"up", "eventWeight_notau * tauWeightIdVsJet[0][5] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauidjetStat1dm1"+tauYear+"down", "eventWeight_notau * tauWeightIdVsJet[0][6] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauidjetStat1dm10"+tauYear+"up", "eventWeight_notau * tauWeightIdVsJet[0][7] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauidjetStat1dm10"+tauYear+"down", "eventWeight_notau * tauWeightIdVsJet[0][8] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauidjetStat1dm11"+tauYear+"up", "eventWeight_notau * tauWeightIdVsJet[0][9] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauidjetStat1dm11"+tauYear+"down", "eventWeight_notau * tauWeightIdVsJet[0][10] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauidjetStat2dm0"+tauYear+"up", "eventWeight_notau * tauWeightIdVsJet[0][11] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauidjetStat2dm0"+tauYear+"down", "eventWeight_notau * tauWeightIdVsJet[0][12] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauidjetStat2dm1"+tauYear+"up", "eventWeight_notau * tauWeightIdVsJet[0][13] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauidjetStat2dm1"+tauYear+"down", "eventWeight_notau * tauWeightIdVsJet[0][14] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauidjetStat2dm10"+tauYear+"up", "eventWeight_notau * tauWeightIdVsJet[0][15] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauidjetStat2dm10"+tauYear+"down", "eventWeight_notau * tauWeightIdVsJet[0][16] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauidjetStat2dm11"+tauYear+"up", "eventWeight_notau * tauWeightIdVsJet[0][17] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauidjetStat2dm11"+tauYear+"down", "eventWeight_notau * tauWeightIdVsJet[0][18] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauidjetSystallerasup", "eventWeight_notau * tauWeightIdVsJet[0][19] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauidjetSystallerasdown", "eventWeight_notau * tauWeightIdVsJet[0][20] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauidjetSyst"+tauYear+"up", "eventWeight_notau * tauWeightIdVsJet[0][21] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauidjetSyst"+tauYear+"down", "eventWeight_notau * tauWeightIdVsJet[0][22] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauidjetSystTESdm0"+tauYear+"up", "eventWeight_notau * tauWeightIdVsJet[0][23] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauidjetSystTESdm0"+tauYear+"down", "eventWeight_notau * tauWeightIdVsJet[0][24] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauidjetSystTESdm1"+tauYear+"up", "eventWeight_notau * tauWeightIdVsJet[0][25] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauidjetSystTESdm1"+tauYear+"down", "eventWeight_notau * tauWeightIdVsJet[0][26] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauidjetSystTESdm10"+tauYear+"up", "eventWeight_notau * tauWeightIdVsJet[0][27] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauidjetSystTESdm10"+tauYear+"down", "eventWeight_notau * tauWeightIdVsJet[0][28] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauidjetSystTESdm11"+tauYear+"up", "eventWeight_notau * tauWeightIdVsJet[0][29] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauidjetSystTESdm11"+tauYear+"down", "eventWeight_notau * tauWeightIdVsJet[0][30] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauidjetHighptstat_bin1up", "eventWeight_notau * tauWeightIdVsJet[0][31] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauidjetHighptstat_bin1down", "eventWeight_notau * tauWeightIdVsJet[0][32] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauidjetHighptstat_bin2up", "eventWeight_notau * tauWeightIdVsJet[0][33] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauidjetHighptstat_bin2down", "eventWeight_notau * tauWeightIdVsJet[0][34] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauidjetHighptsystup", "eventWeight_notau * tauWeightIdVsJet[0][35] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauidjetHighptsystdown", "eventWeight_notau * tauWeightIdVsJet[0][36] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauidjetHighptextrapup", "eventWeight_notau * tauWeightIdVsJet[0][37] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauidjetHighptextrapdown", "eventWeight_notau * tauWeightIdVsJet[0][38] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][0]"});
            // Legacy aliases for backward compatibility
            addVar({"eventWeight__tauidjetUncert0up", "eventWeight__tauidjetStat1dm0"+tauYear+"up"});
            addVar({"eventWeight__tauidjetUncert0down", "eventWeight__tauidjetStat1dm0"+tauYear+"down"});
            addVar({"eventWeight__tauidjetUncert1up", "eventWeight__tauidjetStat1dm1"+tauYear+"up"});
            addVar({"eventWeight__tauidjetUncert1down", "eventWeight__tauidjetStat1dm1"+tauYear+"down"});
            addVar({"eventWeight__tauidjetSystdm0"+tauYear+"up", "eventWeight__tauidjetStat1dm0"+tauYear+"up"});
            addVar({"eventWeight__tauidjetSystdm0"+tauYear+"down", "eventWeight__tauidjetStat1dm0"+tauYear+"down"});
            addVar({"eventWeight__tauidjetSystdm1"+tauYear+"up", "eventWeight__tauidjetStat1dm1"+tauYear+"up"});
            addVar({"eventWeight__tauidjetSystdm1"+tauYear+"down", "eventWeight__tauidjetStat1dm1"+tauYear+"down"});
            addVar({"eventWeight__tauidjetSystdm10"+tauYear+"up", "eventWeight__tauidjetStat1dm10"+tauYear+"up"});
            addVar({"eventWeight__tauidjetSystdm10"+tauYear+"down", "eventWeight__tauidjetStat1dm10"+tauYear+"down"});
            addVar({"eventWeight__tauidjetSystdm11"+tauYear+"up", "eventWeight__tauidjetStat1dm11"+tauYear+"up"});
            addVar({"eventWeight__tauidjetSystdm11"+tauYear+"down", "eventWeight__tauidjetStat1dm11"+tauYear+"down"});
            addVar({"eventWeight__tauidelup", "eventWeight_notau * tauWeightIdVsJet[0][0] * tauWeightIdVsEl[0][1] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauideldown", "eventWeight_notau * tauWeightIdVsJet[0][0] * tauWeightIdVsEl[0][2] * tauWeightIdVsMu[0][0]"});
            addVar({"eventWeight__tauidmuup", "eventWeight_notau * tauWeightIdVsJet[0][0] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][1]"});
            addVar({"eventWeight__tauidmudown", "eventWeight_notau * tauWeightIdVsJet[0][0] * tauWeightIdVsEl[0][0] * tauWeightIdVsMu[0][2]"});
            addVar({"eventWeight__btaghfup", "eventWeight_nobtag * btagWeightNorm[1]"});
            addVar({"eventWeight__btaghfdown", "eventWeight_nobtag * btagWeightNorm[2]"});
            addVar({"eventWeight__btaglfup", "eventWeight_nobtag * btagWeightNorm[3]"});
            addVar({"eventWeight__btaglfdown", "eventWeight_nobtag * btagWeightNorm[4]"});
            addVar({"eventWeight__btaghfstats1up", "eventWeight_nobtag * btagWeightNorm[5]"});
            addVar({"eventWeight__btaghfstats1down", "eventWeight_nobtag * btagWeightNorm[6]"});
            addVar({"eventWeight__btaghfstats2up", "eventWeight_nobtag * btagWeightNorm[7]"});
            addVar({"eventWeight__btaghfstats2down", "eventWeight_nobtag * btagWeightNorm[8]"});
            addVar({"eventWeight__btaglfstats1up", "eventWeight_nobtag * btagWeightNorm[9]"});
            addVar({"eventWeight__btaglfstats1down", "eventWeight_nobtag * btagWeightNorm[10]"});
            addVar({"eventWeight__btaglfstats2up", "eventWeight_nobtag * btagWeightNorm[11]"});
            addVar({"eventWeight__btaglfstats2down", "eventWeight_nobtag * btagWeightNorm[12]"});
            addVar({"eventWeight__btagcferr1up", "eventWeight_nobtag * btagWeightNorm[13]"});
            addVar({"eventWeight__btagcferr1down", "eventWeight_nobtag * btagWeightNorm[14]"});
            addVar({"eventWeight__btagcferr2up", "eventWeight_nobtag * btagWeightNorm[15]"});
            addVar({"eventWeight__btagcferr2down", "eventWeight_nobtag * btagWeightNorm[16]"});

            if (_applytauFF) {
                addVar({"eventWeight__tauFFstatup", "eventWeight * tauFFstatup"});
                addVar({"eventWeight__tauFFstatdown", "eventWeight * tauFFstatdown"});
                addVar({"eventWeight__tauFFsystup", "eventWeight * tauFFsystup"});
                addVar({"eventWeight__tauFFsystdown", "eventWeight * tauFFsystdown"});
                //addVar({"eventWeight__tauFFptdepup", "eventWeight * tauSFAddUncUp"});
                //addVar({"eventWeight__tauFFptdepdown", "eventWeight * tauSFAddUncDn"});
            }

            // no tau - nominal is eventWeight_notau
            // Use unitGenWeight (NOT unitGenWeightFF) because these weights explicitly
            // exclude tau SFs; tauFF is a tau-related correction.
            addVar({"eventWeight_notau_nopu", "unitGenWeight * TopPtWeight[0] * eventWeight_lep * btagWeightNorm[0]"});
            addVar({"eventWeight_notau__puup", "unitGenWeight * TopPtWeight[0] * puWeight[1] * eventWeight_lep * btagWeightNorm[0]"});
            addVar({"eventWeight_notau__pudown", "unitGenWeight * TopPtWeight[0] * puWeight[2] * eventWeight_lep * btagWeightNorm[0]"});
            addVar({"eventWeight_notau__topptup", "unitGenWeight * TopPtWeight[1] * puWeight[0] * eventWeight_lep * btagWeightNorm[0]"});
            addVar({"eventWeight_notau__topptdown", "unitGenWeight * TopPtWeight[2] * puWeight[0] * eventWeight_lep * btagWeightNorm[0]"});

            // Lepton SF variations (no-tau denominator) — channel-dependent
            if (_isMuonCh) {
                addVar({"eventWeight_notau__muidup",    "eventWeight_genpu * muonWeightId[1] * muonWeightIso[0] * muonWeightTrg[0] * btagWeightNorm[0]"});
                addVar({"eventWeight_notau__muiddown",  "eventWeight_genpu * muonWeightId[2] * muonWeightIso[0] * muonWeightTrg[0] * btagWeightNorm[0]"});
                addVar({"eventWeight_notau__muisoup",   "eventWeight_genpu * muonWeightId[0] * muonWeightIso[1] * muonWeightTrg[0] * btagWeightNorm[0]"});
                addVar({"eventWeight_notau__muisodown", "eventWeight_genpu * muonWeightId[0] * muonWeightIso[2] * muonWeightTrg[0] * btagWeightNorm[0]"});
                addVar({"eventWeight_notau__mutrgup",   "eventWeight_genpu * muonWeightId[0] * muonWeightIso[0] * muonWeightTrg[1] * btagWeightNorm[0]"});
                addVar({"eventWeight_notau__mutrgdown", "eventWeight_genpu * muonWeightId[0] * muonWeightIso[0] * muonWeightTrg[2] * btagWeightNorm[0]"});
            } else {
                addVar({"eventWeight_notau__elecrecoup",   "eventWeight_genpu * elecWeightReco[1] * elecWeightId[0] * elecWeightTrg[0] * btagWeightNorm[0]"});
                addVar({"eventWeight_notau__elecrecodown", "eventWeight_genpu * elecWeightReco[2] * elecWeightId[0] * elecWeightTrg[0] * btagWeightNorm[0]"});
                addVar({"eventWeight_notau__elecidup",     "eventWeight_genpu * elecWeightReco[0] * elecWeightId[1] * elecWeightTrg[0] * btagWeightNorm[0]"});
                addVar({"eventWeight_notau__eleciddown",   "eventWeight_genpu * elecWeightReco[0] * elecWeightId[2] * elecWeightTrg[0] * btagWeightNorm[0]"});
                addVar({"eventWeight_notau__electrgup",    "eventWeight_genpu * elecWeightReco[0] * elecWeightId[0] * elecWeightTrg[1] * btagWeightNorm[0]"});
                addVar({"eventWeight_notau__electrgdown",  "eventWeight_genpu * elecWeightReco[0] * elecWeightId[0] * elecWeightTrg[2] * btagWeightNorm[0]"});
            }
            // muonHighPtAddUncUp/Dn removed: no official Run3 high-pT muon additional uncertainty prescription.

            addVar({"eventWeight_notau__btaghfup", "eventWeight_genpumu * btagWeightNorm[1]"});
            addVar({"eventWeight_notau__btaghfdown", "eventWeight_genpumu * btagWeightNorm[2]"});
            addVar({"eventWeight_notau__btaglfup", "eventWeight_genpumu * btagWeightNorm[3]"});
            addVar({"eventWeight_notau__btaglfdown", "eventWeight_genpumu * btagWeightNorm[4]"});
            addVar({"eventWeight_notau__btaghfstats1up", "eventWeight_genpumu * btagWeightNorm[5]"});
            addVar({"eventWeight_notau__btaghfstats1down", "eventWeight_genpumu * btagWeightNorm[6]"});
            addVar({"eventWeight_notau__btaghfstats2up", "eventWeight_genpumu * btagWeightNorm[7]"});
            addVar({"eventWeight_notau__btaghfstats2down", "eventWeight_genpumu * btagWeightNorm[8]"});
            addVar({"eventWeight_notau__btaglfstats1up", "eventWeight_genpumu * btagWeightNorm[9]"});
            addVar({"eventWeight_notau__btaglfstats1down", "eventWeight_genpumu * btagWeightNorm[10]"});
            addVar({"eventWeight_notau__btaglfstats2up", "eventWeight_genpumu * btagWeightNorm[11]"});
            addVar({"eventWeight_notau__btaglfstats2down", "eventWeight_genpumu * btagWeightNorm[12]"});
            addVar({"eventWeight_notau__btagcferr1up", "eventWeight_genpumu * btagWeightNorm[13]"});
            addVar({"eventWeight_notau__btagcferr1down", "eventWeight_genpumu * btagWeightNorm[14]"});
            addVar({"eventWeight_notau__btagcferr2up", "eventWeight_genpumu * btagWeightNorm[15]"});
            addVar({"eventWeight_notau__btagcferr2down", "eventWeight_genpumu * btagWeightNorm[16]"});

            if (_syst == "theory") {
                // ME: [0] is renscfact=0.5d0 facscfact=0.5d0 ; [1] is renscfact=0.5d0 facscfact=1d0 ; [3] is renscfact=1d0 facscfact=0.5d0 ;
                //     [5] is renscfact=1d0 facscfact=2d0 ; [7] is renscfact=2d0 facscfact=1d0 ; [8] is renscfact=2d0 facscfact=2d0
                addVar({"eventWeight__mescaledown", "eventWeight * " + si("LHEScaleWeight", 0, "1.f")});
                addVar({"eventWeight__renscaledown", "eventWeight * " + si("LHEScaleWeight", 1, "1.f")});
                addVar({"eventWeight__facscaledown", "eventWeight * " + si("LHEScaleWeight", 3, "1.f")});
                addVar({"eventWeight__facscaleup",   "eventWeight * " + si("LHEScaleWeight", 5, "1.f")});
                addVar({"eventWeight__renscaleup",   "eventWeight * " + si("LHEScaleWeight", 7, "1.f")});
                addVar({"eventWeight__mescaleup",     "eventWeight * " + si("LHEScaleWeight", 8, "1.f")});
                // notau
                addVar({"eventWeight_notau__mescaledown", "eventWeight_notau * " + si("LHEScaleWeight", 0, "1.f")});
                addVar({"eventWeight_notau__renscaledown", "eventWeight_notau * " + si("LHEScaleWeight", 1, "1.f")});
                addVar({"eventWeight_notau__facscaledown", "eventWeight_notau * " + si("LHEScaleWeight", 3, "1.f")});
                addVar({"eventWeight_notau__facscaleup",   "eventWeight_notau * " + si("LHEScaleWeight", 5, "1.f")});
                addVar({"eventWeight_notau__renscaleup",   "eventWeight_notau * " + si("LHEScaleWeight", 7, "1.f")});
                addVar({"eventWeight_notau__mescaleup",     "eventWeight_notau * " + si("LHEScaleWeight", 8, "1.f")});
                // PS: [0] is ISR=2 FSR=1; [1] is ISR=1 FSR=2; [2] is ISR=0.5 FSR=1; [3] is ISR=1 FSR=0.5;
                addVar({"eventWeight__isrup",   "eventWeight * " + si("PSWeight", 0, "1.f")});
                addVar({"eventWeight__fsrup",   "eventWeight * " + si("PSWeight", 1, "1.f")});
                addVar({"eventWeight__isrdown", "eventWeight * " + si("PSWeight", 2, "1.f")});
                addVar({"eventWeight__fsrdown", "eventWeight * " + si("PSWeight", 3, "1.f")});
                // notau
                addVar({"eventWeight_notau__isrup",   "eventWeight_notau * " + si("PSWeight", 0, "1.f")});
                addVar({"eventWeight_notau__fsrup",   "eventWeight_notau * " + si("PSWeight", 1, "1.f")});
                addVar({"eventWeight_notau__isrdown", "eventWeight_notau * " + si("PSWeight", 2, "1.f")});
                addVar({"eventWeight_notau__fsrdown", "eventWeight_notau * " + si("PSWeight", 3, "1.f")});
                // PDF: LHA IDs 306000 - 306102. 306000 = nominal.
                for (int i=1; i<=102; i++) {
                    addVar({"eventWeight__pdf" + std::to_string(i), "eventWeight * " + si("LHEPdfWeight", i, "1.f")});
                    addVar({"eventWeight_notau__pdf" + std::to_string(i), "eventWeight_notau * " + si("LHEPdfWeight", i, "1.f")});
                }
                addVar({"eventWeight__pdfalphasup",     "eventWeight * " + si("LHEPdfWeight", 102, "1.f")});
                addVar({"eventWeight__pdfalphasdown",   "eventWeight * " + si("LHEPdfWeight", 101, "1.f")});
                addVar({"eventWeight_notau__pdfalphasup",   "eventWeight_notau * " + si("LHEPdfWeight", 102, "1.f")});
                addVar({"eventWeight_notau__pdfalphasdown", "eventWeight_notau * " + si("LHEPdfWeight", 101, "1.f")});
            }
        }
    }
    
}

void TopLFVAnalyzer::storeOutputBranches() {
    // define variables that you want to store
    addVartoStore("run");
    addVartoStore("event");
    addVartoStore("luminosityBlock");
    addVartoStore("eventWeight.*");
    addVartoStore("nmuonpass");
    addVartoStore("nelepass");
    addVartoStore("ncleanjetspass");
    addVartoStore("ncleanbjetspass");
    addVartoStore("ncleantaupass");
    addVartoStore("Jet_pt");
    addVartoStore("Jet1_.*");
    addVartoStore("Jet2_.*");
    addVartoStore("Jet3_.*");
    addVartoStore("bJet1_.*");
    addVartoStore("Tau1.*");
    if (_isMuonCh) {
        addVartoStore("Muon1.*");
    } else {
        addVartoStore("Electron1.*");
    }
    addVartoStore("leptau.*");
    //addVartoStore("MET_pt");
    //addVartoStore("MET_phi");
    addVartoStore("PuppiMET_pt");
    addVartoStore("PuppiMET_phi");
    addVartoStore("chi2.*");
    addVartoStore("btag.*");
    addVartoStore("GenPart_top_pt");
    addVartoStore("TopPtWeight");
    addVartoStore("LHEPart_pt");
    addVartoStore("LHEPart_pdgId");
    addVartoStore("tauWeight.*");
    addVartoStore("tauFF.*");
    addVartoStore("isFakeTau");
}

void TopLFVAnalyzer::bookHists() {
    std::string lep = "e", title_tmp = "";
    if (_isMuonCh) lep = "#mu"; 
    std::vector<std::string> init_weight = {""};
    // Base weight suffixes common to both channels
    std::vector<std::string> sf_weight = {
        "",
        "_nobtag", "_nopu", "_notau", "_notoppt",
        "__puup",   "__pudown",
        "__topptup", "__topptdown",
        // btag shape variations (all 16 components)
        "__btaghfup",       "__btaghfdown",      // = btaghf (always defined)
        "__btaglfup",       "__btaglfdown",
        "__btaghfstats1up",  "__btaghfstats1down",
        "__btaghfstats2up",  "__btaghfstats2down",
        "__btaglfstats1up",  "__btaglfstats1down",
        "__btaglfstats2up",  "__btaglfstats2down",
        "__btagcferr1up", "__btagcferr1down",
        "__btagcferr2up", "__btagcferr2down",
    };
    // Channel-dependent lepton SF suffixes
    if (_isMuonCh) {
        sf_weight.insert(sf_weight.end(), {
            "__muidup",   "__muiddown",
            "__muisoup",  "__muisodown",
            "__mutrgup",  "__mutrgdown",
            // Note: __muhighptup/down removed — no official Run3 high-pT muon prescription
        });
    } else {
        sf_weight.insert(sf_weight.end(), {
            "__elecrecoup",  "__elecrecodown",
            "__elecidup",    "__eleciddown",
            "__electrgup",   "__electrgdown",
        });
    }



    std::vector<std::string> sf_weight_tau = {"__tauidjetup", "__tauidjetdown",
                                              "__tauidjetStat1dm0"+tauYear+"up", "__tauidjetStat1dm0"+tauYear+"down",
                                              "__tauidjetStat1dm1"+tauYear+"up", "__tauidjetStat1dm1"+tauYear+"down",
                                              "__tauidjetStat1dm10"+tauYear+"up", "__tauidjetStat1dm10"+tauYear+"down",
                                              "__tauidjetStat1dm11"+tauYear+"up", "__tauidjetStat1dm11"+tauYear+"down",
                                              "__tauidjetStat2dm0"+tauYear+"up", "__tauidjetStat2dm0"+tauYear+"down",
                                              "__tauidjetStat2dm1"+tauYear+"up", "__tauidjetStat2dm1"+tauYear+"down",
                                              "__tauidjetStat2dm10"+tauYear+"up", "__tauidjetStat2dm10"+tauYear+"down",
                                              "__tauidjetStat2dm11"+tauYear+"up", "__tauidjetStat2dm11"+tauYear+"down",
                                              "__tauidjetUncert0up", "__tauidjetUncert0down",
                                              "__tauidjetUncert1up", "__tauidjetUncert1down",
                                              "__tauidjetSystallerasup", "__tauidjetSystallerasdown",
                                              "__tauidjetSyst"+tauYear+"up", "__tauidjetSyst"+tauYear+"down",
                                              "__tauidjetSystdm0"+tauYear+"up", "__tauidjetSystdm0"+tauYear+"down",
                                              "__tauidjetSystdm1"+tauYear+"up", "__tauidjetSystdm1"+tauYear+"down",
                                              "__tauidjetSystdm10"+tauYear+"up", "__tauidjetSystdm10"+tauYear+"down",
                                              "__tauidjetSystdm11"+tauYear+"up", "__tauidjetSystdm11"+tauYear+"down",
                                              "__tauidjetSystTESdm0"+tauYear+"up", "__tauidjetSystTESdm0"+tauYear+"down",
                                              "__tauidjetSystTESdm1"+tauYear+"up", "__tauidjetSystTESdm1"+tauYear+"down",
                                              "__tauidjetSystTESdm10"+tauYear+"up", "__tauidjetSystTESdm10"+tauYear+"down",
                                              "__tauidjetSystTESdm11"+tauYear+"up", "__tauidjetSystTESdm11"+tauYear+"down",
                                              "__tauidjetHighptstat_bin1up", "__tauidjetHighptstat_bin1down",
                                              "__tauidjetHighptstat_bin2up", "__tauidjetHighptstat_bin2down",
                                              "__tauidjetHighptsystup", "__tauidjetHighptsystdown",
                                              "__tauidjetHighptextrapup", "__tauidjetHighptextrapdown",
                                              "__tauidelup", "__tauideldown", "__tauidmuup", "__tauidmudown"};

    std::vector<std::string> sf_weight_FF ={"__tauFFstatup", "__tauFFstatdown", "__tauFFsystup", "__tauFFsystdown"};

    std::vector<std::string> theory_weight = {"__isrup", "__fsrup", "__isrdown", "__fsrdown",
                   "__mescaleup", "__mescaledown", "__renscaleup", "__renscaledown", "__facscaleup", "__facscaledown",
                   "__pdfalphasup", "__pdfalphasdown"};
                   //Not including pdf eigenvectors for plots

    //Step definition. To replace for tauFF
    //This implies histos w/o FF must be drawn, for bSF and other weights renormalization
    std::string minstep_S1 = "0";
    std::string minstep_S2 = "00";
    std::string minstep_S3 = "000";
    std::string minstep_S4 = "0000";
    std::string minstep_S5 = "00000";

    if (_applytauFF) {
        minstep_S1 = "00000";
        minstep_S2 = "00000";
        //minstep_S3 = "00000";
        minstep_S4 = "00000";
        minstep_S5 = "00000";
    }

    // TODO refine this. Too heavy? we will see.
    std::vector<std::string> syst_weight;
    if (_syst == "" or _syst == "nosyst" or _syst == "data" or ext_syst) {
        syst_weight = init_weight;
        if (_syst != "data") {
            //We anyway need this for bSF rescaling
            add1DHist({"h_nevents", ";Number of events w/o b SF;Events", 2, -0.5, 1.5}, "one", "eventWeight", "_nobtag", minstep_S1, "");
            add1DHist({"h_nevents_notausf", ";Number of events w/o b and tau SF;Events", 2, -0.5, 1.5}, "one", "eventWeight_notau", "_nobtag", minstep_S1, "00");
            add1DHist({"h_jet_ht", ";Jet HT w/o b SF (GeV);Events", 48, 40, 1000}, "Jet_HT", "eventWeight", "_nobtag", minstep_S1, "");
        }
    }
    else {
        syst_weight = sf_weight;
        if (_syst == "theory") syst_weight.insert(syst_weight.end(), theory_weight.begin(), theory_weight.end());
    }

    // S1 w/o tau SF
    maxstep = "00"; //Must be +1 step than its cut
    for (std::string weightstr : syst_weight) {
        if (weightstr.find("notau") != std::string::npos) continue;
        if (weightstr.find("notoppt") != std::string::npos) continue;

        add1DHist({"h_nevents_notausf", ";Number of events;Events", 2, -0.5, 1.5}, "one", "eventWeight_notau", weightstr, minstep_S1, maxstep);
        add1DHist({"h_nvtx_notausf", ";Number of primary vertex;Events", 70, 0.0, 70.0}, "PV_npvsGood", "eventWeight_notau", weightstr, minstep_S1, maxstep);

        add1DHist({"h_puppimet_pt_notausf", ";PuppiMET (GeV);Events", 20, 0, 400}, "PuppiMET_pt", "eventWeight_notau", weightstr, minstep_S1, maxstep);
        add1DHist({"h_puppisum_et_notausf", ";Sum ET;Events", 50, 0.0, 5000.0}, "PuppiMET_sumEt", "eventWeight_notau", weightstr, minstep_S1, maxstep);
        add1DHist({"h_puppimet_phi_notausf", ";PuppiMET #phi;Events", 20, -3.2, 3.2}, "PuppiMET_phi", "eventWeight_notau", weightstr, minstep_S1, maxstep);

        add1DHist({"h_ncleantaupass_notausf", ";#tau_{h} multiplicity;Events", 5, 0.0, 5.0}, "ncleantaupass", "eventWeight_notau", weightstr, minstep_S1, maxstep);
        add1DHist({"h_ncleanjetspass_notausf", ";Jet multiplicity;Events", 10, 0.0, 10.0}, "ncleanjetspass", "eventWeight_notau", weightstr, minstep_S1, maxstep);
        add1DHist({"h_ncleanbjetspass_notausf", ";b-tagged jet multiplicity;Events", 5, 0.0, 5.0}, "ncleanbjetspass", "eventWeight_notau", weightstr, minstep_S1, maxstep);

        if (_isMuonCh) { 
            add1DHist({"h_nmuonpass_notausf", ";Muon multiplicity;Events", 5, 0.0, 5.0}, "nmuonpass", "eventWeight_notau", weightstr, minstep_S1, maxstep);
            add1DHist({"h_muon1_pt_notausf", ";Muon p_{T} (GeV);Events", 10, 0, 200}, "Muon1_pt", "eventWeight_notau", weightstr, minstep_S1, maxstep);
            add1DHist({"h_muon1_eta_notausf", ";Muon #eta;Events", 20, -2.4, 2.4}, "Muon1_eta", "eventWeight_notau", weightstr, minstep_S1, maxstep);
        } else {
            add1DHist({"h_nelepass_notausf", ";electron multiplicity;Events", 5, 0.0, 5.0}, "nelepass", "eventWeight_notau", weightstr, minstep_S1, maxstep);
            add1DHist({"h_electron1_pt_notausf", ";Electron p_{T} (GeV);Events", 10, 0, 200}, "Electron1_pt", "eventWeight_notau", weightstr, minstep_S1, maxstep);
            add1DHist({"h_electron1_eta_notausf", ";Electron #eta;Events", 20, -2.4, 2.4}, "Electron1_eta", "eventWeight_notau", weightstr, minstep_S1, maxstep);
        }
        title_tmp = std::string(";m_{T}(")+lep+", MET) (GeV);Events";
        add1DHist({"h_lepMET_mt_notausf", title_tmp.c_str(), 20, 0, 400}, "lepMET_mt", "eventWeight_notau", weightstr, minstep_S1, maxstep);

        add1DHist({"h_jet1_pt_notausf", ";Leading jet p_{T} (GeV);Events", 20, 0, 400}, "Jet1_pt", "eventWeight_notau", weightstr, minstep_S1, maxstep);
        add1DHist({"h_jet1_eta_notausf", ";Leading jet #eta;Events", 20, -2.4, 2.4}, "Jet1_eta", "eventWeight_notau", weightstr, minstep_S1, maxstep);
        add1DHist({"h_jet1_mass_notausf", ";Leading jet mass (GeV);Events", 20, 0, 100}, "Jet1_mass", "eventWeight_notau", weightstr, minstep_S1, maxstep);
        add1DHist({"h_jet1_btag_notausf",";PNet score of leading jet;Events", 20, 0, 1.0}, "Jet1_btagPNetB", "eventWeight_notau", weightstr, minstep_S1, maxstep);
    }

    //for all the other nominal histograms with tauSF
    if (_syst == "all" or _syst == "theory") {
        syst_weight.insert(syst_weight.end(), sf_weight_tau.begin(), sf_weight_tau.end());
        if (_applytauFF) syst_weight.insert(syst_weight.end(), sf_weight_FF.begin(), sf_weight_FF.end());
    }

    cout << "Variations to take care :";
    for (auto i : syst_weight) cout << i << " ";
    cout << endl;

    maxstep = "";

    for (std::string weightstr : syst_weight) {
        if (weightstr.find("notoppt") != std::string::npos) continue;

        add1DHist({"h_nevents", ";Number of events;Events", 2, -0.5, 1.5}, "one", "eventWeight", weightstr, minstep_S1, maxstep);
        add1DHist({"h_nvtx", ";Number of primary vertex;Events", 70, 0.0, 70.0}, "PV_npvsGood", "eventWeight", weightstr, minstep_S1, maxstep);
        add1DHist({"h_puppimet_pt", ";PuppiMET (GeV);Events", 20, 0, 400}, "PuppiMET_pt", "eventWeight", weightstr, minstep_S1, maxstep);
        add1DHist({"h_puppisum_et", ";Sum ET;Events", 50, 0.0, 5000.0}, "PuppiMET_sumEt", "eventWeight", weightstr, minstep_S1, maxstep);
        add1DHist({"h_puppimet_phi", ";PuppiMET #phi;Events", 20, -3.2, 3.2}, "PuppiMET_phi", "eventWeight", weightstr, minstep_S1, maxstep);

        add1DHist({"h_ncleantaupass", ";#tau_{h} multiplicity;Events", 5, 0.0, 5.0}, "ncleantaupass", "eventWeight", weightstr, minstep_S1, maxstep);
        add1DHist({"h_ncleanjetspass", ";Jet multiplicity;Events", 10, 0.0, 10.0}, "ncleanjetspass", "eventWeight", weightstr, minstep_S1, maxstep);
        add1DHist({"h_ncleanbjetspass", ";b-tagged jet multiplicity;Events", 5, 0.0, 5.0}, "ncleanbjetspass", "eventWeight", weightstr, minstep_S1, maxstep);

        if (_isMuonCh) {
            add1DHist({"h_nmuonpass", ";Muon multiplicity;Events", 5, 0.0, 5.0}, "nmuonpass", "eventWeight", weightstr, minstep_S1, maxstep);
            add1DHist({"h_muon1_pt", ";Muon p_{T} (GeV);Events", 10, 0, 200}, "Muon1_pt", "eventWeight", weightstr, minstep_S1, maxstep);
            add1DHist({"h_muon1_eta", ";Muon #eta;Events", 20, -2.4, 2.4}, "Muon1_eta", "eventWeight", weightstr, minstep_S1, maxstep);
        } else {
            add1DHist({"h_nelepass", ";electron multiplicity;Events", 5, 0.0, 5.0}, "nelepass", "eventWeight", weightstr, minstep_S1, maxstep);
            add1DHist({"h_electron1_pt", ";Electron p_{T} (GeV);Events", 10, 0, 200}, "Electron1_pt", "eventWeight", weightstr, minstep_S1, maxstep);
            add1DHist({"h_electron1_eta", ";Electron #eta;Events", 20, -2.4, 2.4}, "Electron1_eta", "eventWeight", weightstr, minstep_S1, maxstep);
        }

        title_tmp =  std::string(";m_{T}(")+lep+", MET) (GeV);Events";
        add1DHist({"h_lepMET_mt", title_tmp.c_str(), 20, 0, 400}, "lepMET_mt", "eventWeight", weightstr, minstep_S1, maxstep);

        add1DHist({"h_tau1_pt", ";#tau_{h} p_{T} (GeV);Events", 10, 0, 200}, "Tau1_pt", "eventWeight", weightstr, minstep_S2, maxstep);
        add1DHist({"h_tau1_eta", ";#tau_{h} #eta;Events", 20, -2.3, 2.3}, "Tau1_eta", "eventWeight", weightstr, minstep_S2, maxstep);
        add1DHist({"h_tau1_mass", ";m_{#tau_{h}} (GeV);Events", 20, 0, 100}, "Tau1_mass", "eventWeight", weightstr, minstep_S2, maxstep);
        add1DHist({"h_tau1_decayMode", ";Decaymode of #tau_{h};Events", 15, 0, 15}, "Tau1_decayMode", "eventWeight", weightstr, minstep_S2, maxstep);

        title_tmp = std::string(";#Delta#eta(")+lep+", #tau_{h});Events";
        add1DHist({"h_leptau_dEta", title_tmp.c_str(), 25, -5, 5}, "leptau_dEta", "eventWeight", weightstr, minstep_S2, maxstep);
        title_tmp = std::string(";#Delta#phi(")+lep+", #tau_{h});Events";
        add1DHist({"h_leptau_dPhi", title_tmp.c_str(), 20, -3.2, 3.2}, "leptau_dPhi", "eventWeight", weightstr, minstep_S2, maxstep);
        title_tmp = std::string(";#Delta R(")+lep+", #tau_{h});Events";
        add1DHist({"h_leptau_dR", title_tmp.c_str(), 20, 0, 4.0}, "leptau_dR", "eventWeight", weightstr, minstep_S2, maxstep);
        title_tmp = std::string(";Mass of ")+lep+" + #tau_{h} (GeV);Events";
        add1DHist({"h_leptau_mass", title_tmp.c_str(), 20, 0, 500}, "leptau_mass", "eventWeight", weightstr, minstep_S2, maxstep);
        add1DHist({"h_leptau_mass_blind", title_tmp.c_str(), 20, 0, 500}, "leptau_mass_blind", "eventWeight", weightstr, minstep_S2, maxstep);

        add1DHist({"h_jet1_pt", ";Leading jet p_{T} (GeV);Events", 20, 0, 400}, "Jet1_pt", "eventWeight", weightstr, minstep_S1, maxstep);
        add1DHist({"h_jet1_eta", ";Leading jet #eta;Events", 20, -2.4, 2.4}, "Jet1_eta", "eventWeight", weightstr, minstep_S1, maxstep);
        add1DHist({"h_jet1_mass", ";Leading jet mass (GeV);Events", 20, 0, 100}, "Jet1_mass", "eventWeight", weightstr, minstep_S1, maxstep);
        add1DHist({"h_jet1_btag",";PNet score of leading jet;Events", 20, 0, 1.0}, "Jet1_btagPNetB", "eventWeight", weightstr, minstep_S1, maxstep);

        add1DHist({"h_jet2_pt", ";Subleading jet p_{T} (GeV);Events", 20, 0, 400}, "Jet2_pt", "eventWeight", weightstr, minstep_S4, maxstep);
        add1DHist({"h_jet2_eta", ";Subleading jet #eta;Events", 20, -2.4, 2.4}, "Jet2_eta", "eventWeight", weightstr, minstep_S4, maxstep);
        add1DHist({"h_jet2_mass", ";Subleading jet mass (GeV);Events", 20, 0, 100}, "Jet2_mass", "eventWeight", weightstr, minstep_S4, maxstep);
        add1DHist({"h_jet2_btag",";PNet score of subleading jet;Events", 20, 0, 1.0}, "Jet2_btagPNetB", "eventWeight", weightstr, minstep_S4, maxstep);

        add1DHist({"h_jet3_pt", ";Subsubleading jet p_{T} (GeV);Events", 20, 0, 400}, "Jet3_pt", "eventWeight", weightstr, minstep_S4, maxstep);
        add1DHist({"h_jet3_eta", ";Subsubleading jet #eta;Events", 20, -2.4, 2.4}, "Jet3_eta", "eventWeight", weightstr, minstep_S4, maxstep);
        add1DHist({"h_jet3_mass", ";Subsubleading jet mass (GeV);Events", 20, 0, 100}, "Jet3_mass", "eventWeight", weightstr, minstep_S4, maxstep);
        add1DHist({"h_jet3_btag",";PNet score of subsubleading jet;Events", 20, 0, 1.0}, "Jet3_btagPNetB", "eventWeight", weightstr, minstep_S4, maxstep);

        add1DHist({"h_bjet1_pt", ";b-tagged jet p_{T} (GeV);Events", 20, 0, 400}, "bJet1_pt", "eventWeight", weightstr, minstep_S5, maxstep);
        add1DHist({"h_bjet1_eta", ";b-tagged jet #eta;Events", 20, -2.4, 2.4}, "bJet1_eta", "eventWeight", weightstr, minstep_S5, maxstep);
        add1DHist({"h_bjet1_mass", ";b-tagged jet mass (GeV);Events", 20, 0, 100}, "bJet1_mass", "eventWeight", weightstr, minstep_S5, maxstep);

        add1DHist({"h_jet_ht", ";Jet HT (GeV);Events", 48, 40, 1000}, "Jet_HT", "eventWeight", weightstr, minstep_S1, maxstep);
        add1DHist({"h_st_met", ";S_{T}^{MET} (GeV);Events", 28, 100, 1500}, "st_met", "eventWeight", weightstr, minstep_S4, maxstep);
        // Fiilld with the same st_met, to be drawn with "WIDTH" option.
        add1DHist({"h_st_met_binwidth", ";S_{T}^{MET};Events / GeV", 28, 100, 1500}, "st_met", "eventWeight", weightstr, minstep_S4, maxstep);

        // Histogram of Top mass reconstruction
        add1DHist({"h_chi2", ";Minimum #chi^{2};Events", 20, 0, 1000}, "chi2", "eventWeight", weightstr, minstep_S5, maxstep);
        add1DHist({"h_chi2_SMTop_mass", ";SM top quark mass (GeV);Events", 20, 80, 880}, "chi2_SMTop_mass", "eventWeight", weightstr, minstep_S5, maxstep);
        add1DHist({"h_chi2_SMTop_pt", ";SM top quark pt (GeV);Events", 20, 0, 400}, "chi2_SMTop_pt", "eventWeight", weightstr, minstep_S5, maxstep);
        add1DHist({"h_chi2_SMTop_pt_notoppt", ";SM top quark pt (GeV);Events", 20, 0, 400}, "chi2_SMTop_pt", "eventWeight_notoppt", "", minstep_S5, maxstep);
        add1DHist({"h_chi2_SMW_mass", ";SM W_{had} mass (GeV);Events", 20, 0, 400}, "chi2_SMW_mass", "eventWeight", weightstr, minstep_S5, maxstep);
        add1DHist({"h_chi2_wqq_dEta", ";#Delta#eta of jets from W;Events", 25, -5, 5}, "chi2_wqq_dEta", "eventWeight", weightstr, minstep_S5, maxstep);
        add1DHist({"h_chi2_wqq_dPhi", ";#Delta#phi of jets from W;Events;", 20, -4, 4}, "chi2_wqq_dPhi", "eventWeight", weightstr, minstep_S5, maxstep);
        add1DHist({"h_chi2_wqq_dR", ";#Delta R of jets from W;Events", 20, 0, 4.0}, "chi2_wqq_dR", "eventWeight", weightstr, minstep_S5, maxstep);
    }
}

double TopLFVAnalyzer::tauFF(std::string year_, std::string unc_, int direction_,
                              floats &tau_pt_, floats &tau_gen_pt_, ints &tau_dm_) {
    // For genuine tau: FF and uncertainty are always 1.0.
    if (tau_gen_pt_.size() > 0) return 1.0;

    if (abs(direction_) != 1 && direction_ != 0) return 1.0;

    // ── Static JSON-backed cache (loaded once per process) ────────
    // Values are identical to the previous hardcoded table; the JSON
    // at data/TauFF/tau_fake_factors.json is the authoritative source.
    // Key layout: map_ff[dm_key][year][unc] where unc ∈ {nom, stat, syst}.
    using YearMap = std::map<std::string, std::map<std::string, double>>;
    static std::map<std::string, YearMap> map_ff;
    static bool loaded = false;

    if (!loaded) {
        const std::string ffPath = "data/TauFF/tau_fake_factors.json";
        std::ifstream fin(ffPath);
        if (!fin.good()) {
            std::cerr << "ERROR tauFF: cannot open " << ffPath
                      << " — falling back to 1.0" << std::endl;
            return 1.0;
        }
        Json::Value root;
        fin >> root;
        for (const auto &dmKey : root.getMemberNames()) {
            if (dmKey.rfind("_comment", 0) == 0) continue; // skip comments
            const Json::Value &dmNode = root[dmKey];
            for (const auto &yr : dmNode.getMemberNames()) {
                const Json::Value &yrNode = dmNode[yr];
                if (yrNode.isNull()) continue; // Run 3 placeholder
                for (const auto &unc : yrNode.getMemberNames()) {
                    map_ff[dmKey][yr][unc] = yrNode[unc].asDouble();
                }
            }
        }
        loaded = true;
    }

    // ── Decay-mode key ────────────────────────────────────────────
    std::string dmKey;
    const int dm = tau_dm_[0];
    if      (dm == 0)  dmKey = "dm0";
    else if (dm == 1)  dmKey = "dm1";
    else if (dm == 10) dmKey = "dm10";
    else if (dm == 11) dmKey = "dm11";
    else {
        std::cout << "tauFF: unknown decay mode " << dm << std::endl;
        return 1.0;
    }

    // ── Era key normalization ──────────────────────────────────
    std::string yrKey = year_;
    if      (year_ == "2022_Summer22"    || year_ == "2022")      yrKey = "2022_preEE";
    else if (year_ == "2022_Summer22EE"  || year_ == "2022EE")    yrKey = "2022_postEE";
    else if (year_ == "2023_Summer23"    || year_ == "2023")      yrKey = "2023_preBPix";
    else if (year_ == "2023_Summer23BPix"|| year_ == "2023BPix")  yrKey = "2023_postBPix";
    else if (year_ == "2024_Summer24")                          yrKey = "2024";

    // Guard for years with null (Run 3) entries
    if (map_ff.find(dmKey) == map_ff.end() ||
        map_ff.at(dmKey).find(yrKey) == map_ff.at(dmKey).end()) {
        std::cout << "tauFF: no entry for dm=" << dmKey
                  << " year=" << yrKey << " — returning 1.0" << std::endl;
        return 1.0;
    }

    // ── Return value ──────────────────────────────────────────────
    const auto &entry = map_ff.at(dmKey).at(yrKey);
    if (unc_ == "nom") {
        return entry.at("nom");
    } else {
        const double frac = entry.at(unc_);          // fractional uncertainty
        if (direction_ ==  1) return 1.0 + frac;
        if (direction_ == -1) return 1.0 - frac;
    }
    return 1.0;
}


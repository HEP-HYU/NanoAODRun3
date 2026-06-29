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

void NanoAODAnalyzerrdframe::calculateTauES(string tauYear, string tauid_vsjet, string tauid_vsmu, string tauid_vse) {
    // Tau SF
    cout << "Loading Tau SF" << endl;

    cout << "Tau ID WP vsJet : " << tauid_vsjet << endl;
    cout << "Tau ID WP vsMuon : " << tauid_vsmu << endl;
    cout << "Tau ID WP vsElectron : " << tauid_vse << endl;


    auto tauSFreader = correction::CorrectionSet::from_file("data/TauIDSFs/tau_DeepTau2018v2p5_" + tauYear + ".json.gz");
    auto _testool  = tauSFreader->at("tau_energy_scale");
    auto _tauVsJet = tauSFreader->at("DeepTau2018v2p5VSjet");
    auto _tauVsEle = tauSFreader->at("DeepTau2018v2p5VSe");
    auto _tauVsMu  = tauSFreader->at("DeepTau2018v2p5VSmu");

    // Tau ID SF
    cout << "Applying Tau ID SF" << endl;
    auto tauIdSF = [this, _tauVsJet, _tauVsEle, _tauVsMu, tauid_vsjet, tauid_vse, tauid_vsmu](floats &pts, floats &etas, uchars &dms, uchars &genids)->floats{
        floats xout;
        xout.reserve(pts.size());
        for (size_t i=0; i<pts.size(); i++){
            float sf = 1.0;
            if (int(dms[i]) != 5 && int(dms[i]) != 6) {
                if (int(genids[i]) == 5) { //genuine tau
                    string flag = "dm";
                    if (pts[i] > 140.) flag = "pt";
                    sf = _tauVsJet->evaluate({pts[i], int(dms[i]), int(genids[i]), tauid_vsjet, tauid_vse, "nom", flag});
                } else if (int(genids[i]) == 1 || int(genids[i]) == 3){ //genuine electron
                    sf = _tauVsEle->evaluate({etas[i], int(dms[i]), int(genids[i]), tauid_vse, "nom"});
                } else if (int(genids[i]) == 2 || int(genids[i]) == 4) { //genuine muon
                    sf = _tauVsMu->evaluate({etas[i], int(genids[i]), tauid_vsmu, tauid_vse, tauid_vsjet, "nom"});
                }
            }
            xout.emplace_back(sf);
        }
        return xout;
    };

    // Tau ES
    cout<<"Applying TauES on Genuine taus"<<endl;
    auto tauES = [this, _testool, tauid_vsjet, tauid_vse](std::string var) {
        return [this, _testool, tauid_vsjet, tauid_vse, var](floats &pts, floats &etas, uchars &dms, uchars &genids)->floats {
            floats xout;
            xout.reserve(pts.size());
            for (size_t i=0; i<pts.size(); i++) {
                float es = 1.0;
                int dm_tmp = int (dms[i]);
                try{
                    if (dm_tmp!=5 && dm_tmp!=6 && int(genids[i])!=2 && int(genids[i])!=4 && int(genids[i])!=6){
                        if (dm_tmp == 2) dm_tmp = 1;
                        es = _testool->evaluate({pts[i], etas[i], dm_tmp, int(genids[i]), "DeepTau2018v2p5", tauid_vsjet, tauid_vse, var});
                    }
                } catch (const exception &e) {
                    std::cout << "e: " << e.what() <<std::endl;
                    std::cout << "pt : " << pts[i] << " eta " << etas[i] << " dm " << int(dms[i]) << " genid " << int(genids[i]) << std::endl;
                }
                xout.emplace_back(es);
            }
            return xout;
        };
    };

    _rlm = _rlm.Define("Tau_pt_uncor", "Tau_pt")
           .Define("Tau_ES_nom", tauES("nom"), {"Tau_pt_uncor", "Tau_eta", "Tau_decayMode", "Tau_genPartFlav"})
           .Define("Tau_ES_up", tauES("up"), {"Tau_pt_uncor", "Tau_eta", "Tau_decayMode", "Tau_genPartFlav"})
           .Define("Tau_ES_down", tauES("down"), {"Tau_pt_uncor", "Tau_eta", "Tau_decayMode", "Tau_genPartFlav"})
           .Redefine("Tau_pt", "Tau_pt_uncor*Tau_ES_nom")
           .Redefine("Tau_mass", "Tau_mass*Tau_ES_nom");
}

void NanoAODAnalyzerrdframe::skimJets(string cut) {
    // input vector: vec[pt][vars]
    // Note: do not skim with exact value of pt!
    //auto skimCol = [this](floatsVec toSkim, ints cut)->floatsVec {
    //    floatsVec out;
    //    for (size_t i=0; i<toSkim.size(); i++) {
    //        if (cut[i] > 0) out.emplace_back(toSkim[i]);
    //    }
    //    return out;
    //};

    // skim jet collection
    // https://twiki.cern.ch/twiki/bin/view/CMS/JetID13p6TeV
    // nanoAOD Flags
    _rlm = _rlm.Define("jetcuts", cut)
               .Redefine("Jet_pt", "Jet_pt[jetcuts]")
               .Redefine("Jet_eta", "Jet_eta[jetcuts]")
               .Redefine("Jet_phi", "Jet_phi[jetcuts]")
               .Redefine("Jet_mass", "Jet_mass[jetcuts]")
               .Redefine("Jet_passJetIdTightLepVeto", "Jet_passJetIdTightLepVeto[jetcuts]")
               .Redefine("Jet_area", "Jet_area[jetcuts]")
               .Redefine("Jet_rawFactor", "Jet_rawFactor[jetcuts]")
               .Redefine("Jet_pt_uncorr", "Jet_pt_uncorr[jetcuts]")
               .Redefine("Jet_btagPNetB", "Jet_btagPNetB[jetcuts]")
               .Redefine("Jet_btagPNetCvB", "Jet_btagPNetCvB[jetcuts]")
               .Redefine("Jet_btagPNetCvL", "Jet_btagPNetCvL[jetcuts]")
               .Redefine("nJet", "int(Jet_pt.size())");
    if (_isRun24){
        _rlm = _rlm.Redefine("Jet_btagUParTAK4B", "Jet_btagUParTAK4B[jetcuts]");
        _rlm = _rlm.Redefine("Jet_btagUParTAK4CvB", "Jet_btagUParTAK4CvB[jetcuts]");
        _rlm = _rlm.Redefine("Jet_btagUParTAK4CvL", "Jet_btagUParTAK4CvL[jetcuts]");
    }
    if (!_isData) {
        //_rlm = _rlm.Redefine("Jet_pt_unc", skimCol, {"Jet_pt_unc", "jetcuts"})
        //           .Redefine("Jet_jer", skimCol, {"Jet_jer", "jetcuts"})
        _rlm = _rlm.Redefine("Jet_hadronFlavour","Jet_hadronFlavour[jetcuts]")
                   .Redefine("Jet_genJetIdx","Jet_genJetIdx[jetcuts]");
    }
}

void NanoAODAnalyzerrdframe::matchGenReco() {
    if (_isMuonCh){
        _rlm = _rlm.Define("FinalGenPart_idx", ::FinalGenPart_idx, {"GenPart_pdgId", "GenPart_genPartIdxMother"});
    } else {
        _rlm = _rlm.Define("FinalGenPart_idx", ::FinalGenPart_idx_elec, {"GenPart_pdgId", "GenPart_genPartIdxMother"});
    }
    _rlm = _rlm.Define("GenPart_SMb_idx", "FinalGenPart_idx[3]")
               .Define("GenPart_SMW1_idx", "FinalGenPart_idx[4]")
               .Define("GenPart_SMW2_idx", "FinalGenPart_idx[5]")
               .Define("GenPart_SMtop_idx", "FinalGenPart_idx[7]")
               .Define("GenPart_d_SMb_idx", "FinalGenPart_idx[8]")
               .Define("GenPart_d_SMW_idx", "FinalGenPart_idx[9]")
               .Define("GenPart_d_SMW1_idx", "FinalGenPart_idx[10]")
               .Define("GenPart_d_SMW2_idx", "FinalGenPart_idx[11]")
               .Define("GenPart_d_SMtop_idx", "FinalGenPart_idx[12]");

    _rlm = _rlm.Define("drmax1", "float(0.15)")
               .Define("drmax2", "float(0.4)")
               .Define("GenJet_SMb_idx",::dRmatching,{"GenPart_SMb_idx","drmax2","GenPart_pt","GenPart_eta","GenPart_phi","GenPart_mass","GenJet_pt","GenJet_eta","GenJet_phi","GenJet_mass"})
               .Define("GenJet_SMW1_idx",::dRmatching,{"GenPart_SMW1_idx","drmax2","GenPart_pt","GenPart_eta","GenPart_phi","GenPart_mass","GenJet_pt","GenJet_eta","GenJet_phi","GenJet_mass"})
               .Define("GenJet_SMW2_idx",::dRmatching,{"GenPart_SMW2_idx","drmax2", "GenPart_pt","GenPart_eta","GenPart_phi","GenPart_mass","GenJet_pt","GenJet_eta","GenJet_phi","GenJet_mass"})
               .Define("Jet_SMb_idx",::dRmatching,{"GenJet_SMb_idx","drmax2","GenJet_pt","GenJet_eta","GenJet_phi","GenJet_mass","Jet_pt","Jet_eta","Jet_phi","Jet_mass"})
               .Define("Jet_SMW1_idx",::dRmatching,{"GenJet_SMW1_idx","drmax2","GenJet_pt","GenJet_eta","GenJet_phi","GenJet_mass","Jet_pt","Jet_eta","Jet_phi","Jet_mass"})
               .Define("Jet_SMW2_idx",::dRmatching,{"GenJet_SMW2_idx","drmax2", "GenJet_pt","GenJet_eta","GenJet_phi","GenJet_mass","Jet_pt","Jet_eta","Jet_phi","Jet_mass"});

    _rlm = _rlm.Define("GenPart_SMbmatched_4vecs", ::gen4vec_withidx, {"GenPart_pt","GenPart_eta","GenPart_phi","GenPart_mass", "GenPart_SMb_idx"})
               .Define("GenPart_SMW1matched_4vecs", ::gen4vec_withidx, {"GenPart_pt","GenPart_eta","GenPart_phi","GenPart_mass", "GenPart_SMW1_idx"})
               .Define("GenPart_SMW2matched_4vecs", ::gen4vec_withidx, {"GenPart_pt","GenPart_eta","GenPart_phi","GenPart_mass", "GenPart_SMW2_idx"});

    _rlm = _rlm.Define("GenPart_d_SMbmatched_4vecs", ::gen4vec_withidx, {"GenPart_pt","GenPart_eta","GenPart_phi","GenPart_mass", "GenPart_d_SMb_idx"})
               .Define("GenPart_d_SMWmatched_4vecs", ::gen4vec_withidx, {"GenPart_pt","GenPart_eta","GenPart_phi","GenPart_mass", "GenPart_d_SMW_idx"})
               .Define("GenPart_d_SMW1matched_4vecs", ::gen4vec_withidx, {"GenPart_pt","GenPart_eta","GenPart_phi","GenPart_mass", "GenPart_d_SMW1_idx"})
               .Define("GenPart_d_SMW2matched_4vecs", ::gen4vec_withidx, {"GenPart_pt","GenPart_eta","GenPart_phi","GenPart_mass", "GenPart_d_SMW2_idx"})
               .Define("GenPart_d_SMtopmatched_4vecs", ::gen4vec_withidx, {"GenPart_pt","GenPart_eta","GenPart_phi","GenPart_mass", "GenPart_d_SMtop_idx"})

               .Define("GenJet_SMbmatched_4vecs", ::gen4vec_withidx, {"GenJet_pt","GenJet_eta","GenJet_phi","GenJet_mass", "GenJet_SMb_idx"})
               .Define("GenJet_SMW1matched_4vecs", ::gen4vec_withidx, {"GenJet_pt","GenJet_eta","GenJet_phi","GenJet_mass", "GenJet_SMW1_idx"})
               .Define("GenJet_SMW2matched_4vecs", ::gen4vec_withidx, {"GenJet_pt","GenJet_eta","GenJet_phi","GenJet_mass", "GenJet_SMW2_idx"})

               .Define("Jet_SMbmatched_4vecs", ::gen4vec_withidx, {"Jet_pt","Jet_eta","Jet_phi","Jet_mass", "Jet_SMb_idx"})
               .Define("Jet_SMW1matched_4vecs", ::gen4vec_withidx, {"Jet_pt","Jet_eta","Jet_phi","Jet_mass", "Jet_SMW1_idx"})
               .Define("Jet_SMW2matched_4vecs", ::gen4vec_withidx, {"Jet_pt","Jet_eta","Jet_phi","Jet_mass", "Jet_SMW2_idx"});

    cout << "match Gen Reco ended" << endl;
}


/*
 * TopLFVAnalyzer.h
 *
 *  Created on: April 9, 2020
 *      Author: Tae Jeong Kim 
 */

#ifndef TopLFVAnalyzer_H_
#define TopLFVAnalyzer_H_

#include "NanoAODAnalyzerrdframe.h"

class TopLFVAnalyzer: public NanoAODAnalyzerrdframe {

public:
    TopLFVAnalyzer(TTree *t, std::string outfilename, std::string year="", std::string ch="", std::string syst="", std::string jsonfname="", bool applytauFF=false, string globaltag="", int nthreads=1);
    void defineCuts() override;
    void defineMoreVars() override; // define higher-level variables from
    void bookHists() override;
    void defineObjectSelection(std::vector<std::string> jes_var) override;

protected:
    bool ext_syst = false;  ///< true when _syst encodes an external variation (JES/JER/TES/btag)
    std::string maxstep;
    std::string tauYear = "";
    bool _applytauFF;
    void defineBTagNormalization();
    /// @name defineMoreVars helpers (called in order; RDataFrame column ordering preserved)
    ///@{
    void defineKinematicVars();     ///< 4-vectors, ΔR, M_T, S_T^MET, top reco, object branches
    void defineWeightVars();        ///< nominal and systematic eventWeight definitions
    void storeOutputBranches();     ///< addVartoStore calls (flat ntuple branch list)
    ///@}
    static double tauFF(std::string year_, std::string ch_, std::string unc_, int direction_, floats &tau_pt_, floats &tau_gen_pt_, ints &tau_dm_);

    struct tauFFfunctor {
        std::string year_;
        std::string ch_;  // "muon" or "electron"
        std::string unc_; // "nom", "stat", or "syst"
        int direction_ = 0;

        tauFFfunctor(std::string year_tmp, std::string ch_tmp, std::string unc_tmp, int direction_tmp):
            year_(year_tmp), ch_(ch_tmp), unc_(unc_tmp), direction_(direction_tmp) {}

        double operator()(floats &tau_pt_, floats &tau_gen_pt_, ints &tau_dm_) {
           return TopLFVAnalyzer::tauFF(year_, ch_, unc_, direction_, tau_pt_, tau_gen_pt_, tau_dm_);
        }
    };
};

#endif /* TopLFVAnalyzer_H_ */

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
    void defineCuts();
    void defineMoreVars(); // define higher-level variables from
    void bookHists();
    bool ext_syst = false;

private:
    std::string _year;
    std::string _ch;
    std::string _syst;
    std::string _outfilename;
    std::string maxstep;
    std::string tauYear = "";
    bool _applytauFF;
    bool _isSignal;
    bool _isMuonCh;
    static double tauFF(std::string year_, std::string unc_, int direction_, floats &tau_pt_, floats &tau_gen_pt_, ints &tau_dm_);

    struct tauFFfunctor {
        std::string year_;
        std::string unc_; // "nom", "stat", or "syst"
        int direction_ = 0;

        tauFFfunctor(std::string year_tmp, std::string unc_tmp, int direction_tmp):
            year_(year_tmp), unc_(unc_tmp), direction_(direction_tmp) {}

        double operator()(floats &tau_pt_, floats &tau_gen_pt_, ints &tau_dm_) {
           return TopLFVAnalyzer::tauFF(year_, unc_, direction_, tau_pt_, tau_gen_pt_, tau_dm_);
        }
    };
};

#endif /* TopLFVAnalyzer_H_ */

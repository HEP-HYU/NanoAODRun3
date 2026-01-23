/*
 * TauFakeFactorAnalyzer.h
 *
 *  Created on: June 14, 2023
 *      Author: Jiwon Park
 */

#ifndef TauFakeFactorAnalyzer_H_
#define TauFakeFactorAnalyzer_H_

#include "NanoAODAnalyzerrdframe.h"

class TauFakeFactorAnalyzer: public NanoAODAnalyzerrdframe {

public:
    TauFakeFactorAnalyzer(TTree *t, std::string outfilename, std::string year="", std::string ch="", std::string syst="", std::string jsonfname="", string globaltag="", int nthreads=1, std::string mode="");
    void defineCuts();
    void defineMoreVars(); // define higher-level variables from
    void bookHists();

private:
    std::string _year;
    std::string _ch;
    std::string _syst;
    std::string _mode;
    std::string maxstep;
};

#endif /* TauFakeFactorAnalyzer_H_ */

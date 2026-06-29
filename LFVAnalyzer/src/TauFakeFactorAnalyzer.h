/*
 * TauFakeFactorAnalyzer.h
 *
 *  Created on: June 14, 2023
 *      Author: Jiwon Park
 */

#ifndef TauFakeFactorAnalyzer_H_
#define TauFakeFactorAnalyzer_H_

#include "NanoAODAnalyzerrdframe.h"
#include "TopLFVAnalyzer.h"

class TauFakeFactorAnalyzer: public TopLFVAnalyzer{

public:
    TauFakeFactorAnalyzer(TTree *t, std::string outfilename, std::string year="", std::string ch="", std::string syst="", std::string jsonfname="", bool applytauFF=false, string globaltag="", int nthreads=1, std::string mode="");
    void defineCuts() override;
    void defineMoreVars() override; // define higher-level variables from
    void bookHists() override;

private:
    std::string _mode;
    std::string maxstep;
};

#endif /* TauFakeFactorAnalyzer_H_ */

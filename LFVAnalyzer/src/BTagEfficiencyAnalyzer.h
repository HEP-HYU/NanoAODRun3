/*
 * BTagEfficiencyAnalyzer.h
 *
 *  BTV Method 1a: C++ RDataFrame-based b-tagging efficiency map extractor.
 *  Inherits from TopLFVAnalyzer to reuse full object selections (muon, elec, tau, jet overlap removal)
 *  while keeping the pre-tag baseline phase space (Njets >= 3, without requiring b-tags).
 */

#ifndef BTagEfficiencyAnalyzer_H_
#define BTagEfficiencyAnalyzer_H_

#include "NanoAODAnalyzerrdframe.h"
#include "TopLFVAnalyzer.h"
#include <string>
#include <vector>

class BTagEfficiencyAnalyzer : public TopLFVAnalyzer {

public:
    BTagEfficiencyAnalyzer(TTree *t, std::string outfilename, std::string year="", std::string ch="", std::string syst="", std::string jsonfname="", bool applytauFF=false, std::string globaltag="", int nthreads=1);
    ~BTagEfficiencyAnalyzer() override = default;

    void defineObjectSelection(std::vector<std::string> jes_var) override;
    void defineCuts() override;
    void defineMoreVars() override;
    void bookHists() override;

private:
    float _btagcut = 0.2450f;
    std::string _btagCol = "Jet_btagPNetB";
    std::string _chName = "muon";
    std::string maxstep = "0000";
};

#endif /* BTagEfficiencyAnalyzer_H_ */

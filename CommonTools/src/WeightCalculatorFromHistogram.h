#ifndef PhysicsTools_NanoAODTools_WeightCalculatorFromHistogram_h
#define PhysicsTools_NanoAODTools_WeightCalculatorFromHistogram_h

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include "TH1.h"

class WeightCalculatorFromHistogram {
 public:
  WeightCalculatorFromHistogram() {}
  // get the weight from the bin content of the passed histogram
  WeightCalculatorFromHistogram(TH1 *histogram, bool verbose=false) : histogram_(histogram), verbose_(verbose) {}
  // get the weight from the bin content of the ratio hist/targethist
  WeightCalculatorFromHistogram(TH1 *hist, TH1* targethist, bool norm=true, bool fixLargeWeights=true, bool verbose=false);
  ~WeightCalculatorFromHistogram() {}
  
  float getWeight(float x, float y=0) const;
  float getWeightErr(float x, float y=0) const;
  
 private:
  std::vector<double> loadVals(TH1 *hist, bool norm=true);
  TH1* ratio(TH1 *hist, TH1* targethist, bool fixLargeWgts);
  void fixLargeWeights(std::vector<double> &weights, double maxshift=0.0025,double hardmax=3);
  double checkIntegral(std::vector<double> wgt1, std::vector<double> wgt2);

  TH1* histogram_;
  std::vector<double> refvals_,targetvals_;
  bool verbose_;
  bool norm_;
};

#endif

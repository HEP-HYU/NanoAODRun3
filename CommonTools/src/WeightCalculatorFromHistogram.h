#ifndef PhysicsTools_NanoAODTools_WeightCalculatorFromHistogram_h
#define PhysicsTools_NanoAODTools_WeightCalculatorFromHistogram_h

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include "TH1.h"

/**
 * @brief Utility class for computing per-event reweighting factors from TH1/TH2 histograms.
 *
 * Supports two construction modes:
 *  - **Direct histogram**: pass a pre-computed ratio histogram. `getWeight()` returns
 *    the bin content at (x) or (x, y).
 *  - **Ratio mode**: pass a numerator and denominator histogram. The class computes
 *    the ratio internally, with optional normalisation and large-weight capping.
 *
 * ### Usage in LFV analysis
 * Used for pileup reweighting (`applyWeights()` in NanoAODAnalyzerrdframe.cpp).
 *
 * ### Known issue (fixed)
 * `getWeightErr()` in 2D mode previously used `GetXaxis()` for the Y bin lookup.
 * This is corrected as of 2026-06. Note: `getWeightErr()` is not called anywhere
 * in the current analysis \u2014 the fix is a correctness improvement for future use.
 */
class WeightCalculatorFromHistogram {
 public:
  WeightCalculatorFromHistogram() {}
  /// @brief Construct from a pre-computed ratio histogram.
  WeightCalculatorFromHistogram(TH1 *histogram, bool verbose=false) : histogram_(histogram), verbose_(verbose) {}
  /// @brief Construct from numerator/denominator pair; computes ratio internally.
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

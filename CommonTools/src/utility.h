/*
 * utility.h
 *
 *  Created on: Dec 4, 2018
 *      Author: suyong
 */

#ifndef UTILITY_H_
#define UTILITY_H_

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "Math/Vector4D.h"
#include <string>
#include <map>

/// @file utility.h
/// @brief Physics utility functions and type aliases for the LFV NanoAOD analysis.
///
/// This header defines:
/// - RVec type aliases used throughout the analysis (floats, ints, bools, …)
/// - Gen-level particle tracing functions (FinalGenPart_idx)
/// - 4-vector construction helpers
/// - Top-quark reconstruction functions (chi-squared minimisation)
/// - Kinematic variable calculators (deltaR, M_T, S_T)
///
/// @warning Many functions here are used inside RDataFrame lambda captures.
///          Their signatures and return types are part of the flat-ntuple contract
///          and must NOT be changed without updating all downstream Define() calls.

extern std::map<int,float> btagNormFactors;
float getBTagNormFactor(int njet);

/// @name RVec type aliases
/// Shorthand names for ROOT::VecOps::RVec specialisations used as column types
/// in RDataFrame. These aliases appear in Define() lambdas and Snapshot() calls.
/// Changing them would break the flat-ntuple column schema.
///@{
using floats    = ROOT::VecOps::RVec<float>;
using floatsVec = ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>>;
using doubles   = ROOT::VecOps::RVec<double>;
using doublesVec= ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>>;
using ints      = ROOT::VecOps::RVec<int>;
using shorts    = ROOT::VecOps::RVec<short>;
using bools     = ROOT::VecOps::RVec<bool>;
using uchars    = ROOT::VecOps::RVec<unsigned char>;
using strings   = ROOT::VecOps::RVec<std::string>;

using FourVector    = ROOT::Math::PtEtaPhiMVector;
using FourVectorVec = std::vector<FourVector>;
using FourVectorRVec= ROOT::VecOps::RVec<FourVector>;
///@}


struct hist1dinfo
{
  ROOT::RDF::TH1DModel hmodel;
  std::string varname;
  std::string weightname;
  std::string systname;
  std::string mincutstep;
  std::string maxcutstep;
} ;

struct hist2dinfo
{
  ROOT::RDF::TH2DModel hmodel;
  std::string varnameX;
  std::string varnameY;
  std::string weightname;
  std::string systname;
  std::string mincutstep;
  std::string maxcutstep;
} ;


struct varinfo
{
  std::string varname;
  std::string vardefinition;
  std::string mincutstep;
};

struct cutinfo
{
  std::string cutdefinition;
  std::string idx;
};


// generates vectors of 4 vectors given vectors of pt, eta, phi, mass
FourVectorVec gen4vec(floats &pt, floats &eta, floats &phi, floats &mass);

FourVectorVec gen4vec_withidx(floats &pt, floats &eta, floats &phi, floats &mass, int idx);

FourVectorVec genmet4vec(float met_pt, float met_phi);

// return a vector size equal to length of x all filled with evWeight value
floats weightv(floats &x, float evWeight);

floats sphericity(FourVectorVec &p);

double foxwolframmoment(int l, FourVectorVec &p, int minj=0, int maxj=-1);

ints good_idx(ints good);

floats lqtop_reconstruction( FourVectorVec &cjet, FourVectorVec &mu, FourVectorVec &tau);

/// @name Top-quark reconstruction
/// Chi-squared-based kinematic top reconstruction for ST-LFV and TT-LFV topologies.
/// Mass constants (MT, MW) are hardcoded in utility.cpp and calibrated to the
/// analysis samples. Do NOT change them without a physics justification and
/// corresponding AN update.
///@{
floats top_reconstruction_STLFV(FourVectorVec &jets, FourVectorVec &bjets, FourVectorVec &muons, FourVectorVec &taus);
floats top_reconstruction_TTLFV(FourVectorVec &jets, FourVectorVec &bjets, FourVectorVec &muons, FourVectorVec &taus);
floats top_reco_products_STLFV(FourVectorVec &jets, FourVectorVec &muons, FourVectorVec &taus, floats topreco);
floats top_reco_products_TTLFV(FourVectorVec &jets, FourVectorVec &muons, FourVectorVec &taus, floats topreco);
///@}

float st_met(floats &pt1, floats &pt2, floats &pt3, float pt4);

float calculate_deltaEta( FourVector &p1, FourVector &p2);

float calculate_deltaPhi( FourVector &p1, FourVector &p2);

float calculate_deltaR( FourVector &p1, FourVector &p2);

float calculate_invMass( FourVector &p1, FourVector &p2);

float calculate_MT( FourVectorVec &muons, float met, float metphi);

FourVector sum_4vec( FourVector &p1, FourVector &p2);

/// |Δφ| between a leading-lepton vector and MET
float dPhi_MET(FourVector &lep, float metphi);

/// Invariant mass of three four-vectors (lep + tau + bjet)
float invMass3(FourVector &p1, FourVector &p2, FourVector &p3);

///@}

floats sort_discriminant( floats discr, floats obj );

ints find_element(ints vec, int a);

ints find_element_binary( ints vec, int a);

/// @name Generator-level particle tracing
/// Traverse the NanoAOD GenPart tree to identify final-state particles
/// from the LFV decay chain.
///@{

/// @brief Find the last copy of a particle with given PDG ID in the gen record.
ints LastGenPart_idx(int target_id, ints GenPart_pdgId, shorts GenPart_genPartIdxMother);
int  lastgenpart_idx(int target_i,  ints GenPart_pdgId, shorts GenPart_genPartIdxMother);
int  firstgenpart_idx(int target_i, ints GenPart_pdgId, shorts GenPart_genPartIdxMother);

/**
 * @brief Trace the muon-channel LFV decay chain in the generator record.
 *
 * Returns a 13-element vector. Each index has a fixed physics meaning:
 *
 * | Index | Particle         | Notes                        |
 * |-------|------------------|------------------------------|
 * |  0    | up-type quark    | from LFV top vertex          |
 * |  1    | muon             | from LFV top vertex          |
 * |  2    | tau              | from LFV top vertex          |
 * |  3    | b quark          | from LFV top vertex          |
 * |  4    | W daughter q1    | from SM top W decay          |
 * |  5    | W daughter q2    | from SM top W decay          |
 * |  6    | LFV top          | t -> q mu tau                |
 * |  7    | SM top           | t -> W b                     |
 * |  8    | W from SM top    |                              |
 * |  9    | b from SM top    |                              |
 * | 10    | lepton from W    |                              |
 * | 11    | neutrino from W  |                              |
 * | 12    | quark1 from W    |                              |
 *
 * @warning Index semantics are FIXED. The flat ntuple and DNN feature list
 *          depend on these positions. Do NOT reorder or resize the return vector.
 */
ints FinalGenPart_idx(ints GenPart_pdgId, shorts GenPart_genPartIdxMother);

/**
 * @brief Electron-channel variant of FinalGenPart_idx.
 * Same index layout as FinalGenPart_idx but traces the e+tau LFV decay chain.
 * @warning Same fixity constraint as FinalGenPart_idx applies.
 */
ints FinalGenPart_idx_elec(ints GenPart_pdgId, shorts GenPart_genPartIdxMother);
///@}

int dRmatching( int origin_i,float maxdR,  floats origin_pt, floats origin_eta, floats origin_phi, floats origin_mass, floats target_pt, floats target_eta, floats target_phi, floats target_mass);

ints dRmatching_binary( int origin_i,float maxdR,  floats origin_pt, floats origin_eta, floats origin_phi, floats origin_mass, floats target_pt, floats target_eta, floats target_phi, floats target_mass);

float SumMass( FourVectorVec object1, FourVectorVec object2, FourVectorVec object3);

FourVector select_leadingvec( FourVectorVec &v );

floats addMuonUnc( floats &input );

floatsVec fixtausf ( floatsVec &input, floats &pts);
#endif /* UTILITY_H_ */

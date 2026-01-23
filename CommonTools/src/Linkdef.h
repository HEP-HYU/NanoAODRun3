/*
 * Linkdef.h
 *
 *  Created on: Oct 16, 2018
 *      Author: suyong
 */

#ifndef LINKDEF_H_
#define LINKDEF_H_

#pragma link C++ class std::vector<ROOT::Math::PtEtaPhiMVector>+;
#pragma link C++ class std::vector<std::vector<float>>+;
#pragma link C++ class ROOT::VecOps::RVec<float>+;
#pragma link C++ class ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>>+;
#pragma link C++ class ROOT::VecOps::RVec<double>+;
#pragma link C++ class ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>>+;
#pragma link C++ class std::vector<ROOT::Math::PtEtaPhiMVector>+;
#pragma link C++ class NanoAODAnalyzerrdframe +;
#pragma link C++ class TopLFVAnalyzer +;
#pragma link C++ class TauFakeFactorAnalyzer +;
#pragma link C++ class SkimEvents +;

#endif /* LINKDEF_H_ */

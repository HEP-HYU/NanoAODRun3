#include <iostream>
#include "GEScaleSyst.h"

using namespace std;

int GEScaleSyst::GetCopyIndex(int mode, int icopy) {
  if (era == "2016_UL_HIPM") {
    if (mode == 0)
      return 81600+icopy;
  } else if (era == "2016_UL") {
    if (mode == 0)
      return 91600+icopy;
  } else if (era == "2017_UL") {
    if (mode == 0)
      return 91700+icopy;
  } else if (era == "2018_UL") {
    if (mode == 0)
      return 91800+icopy;
  }

  return -1;
}

float GEScaleSyst::GetGEScaleKappa(float eta, float phi, int mode, int icopy) {
  int eta_idx = -999;
  int phi_idx = -999;

  for(int ieta=0; ieta<neta; ++ieta) {
    if( eta >= etaBinEdge[ieta] && eta < etaBinEdge[ieta+1] ) {
      eta_idx = ieta;
      break;
    }
  }

  for(int iphi=0; iphi<nphi; ++iphi) {
    if( phi >= phiBinEdge[iphi] && phi < phiBinEdge[iphi+1] ) {
      phi_idx = iphi;
      break;
    }
  }

  if(eta_idx < 0 || phi_idx < 0) {
    cout << "ERROR   GEScaleSyst::GetGEScaleKappa   correction value not found... icopy= " << icopy << ", eta= " << eta << ", phi= " << phi << endl;
    return -1e9;
  }

  return _matrix[GetCopyIndex(mode, icopy)][eta_idx][phi_idx];
}

float GEScaleSyst::GEScaleCorrPt(float pt, float eta, float phi, int charge, int mode, int icopy, bool doOpp) {
  if( _matrix.find(GetCopyIndex(mode, icopy)) == _matrix.end() ) {
    //cout << "ERROR GEScaleSyst::GEScaleCorrPt invalid icopy= " << icopy << endl;
    cout << "ERROR GEScaleSyst::GEScaleCorrPt invalid icopy= " << icopy << endl;
    return -1e9;
  }
  if( pt < 0 ) {
    cout << "ERROR GEScaleSyst::GEScaleCorrPt invalid muon pt= " << pt << endl;
    return -1e9;
  }
  if( fabs(eta) > 2.4 ) {
    cout << "ERROR GEScaleSyst::GEScaleCorrPt invalid muon eta= " << eta << endl;
    return -1e9;
  }
  if( fabs(phi) > M_PI ) {
    //This creates too many lines in log files
    //cout << "ERROR GEScaleSyst::GEScaleCorrPt invalid muon phi= " << phi << endl;
    //return -1e9;
    return pt;
  }
  if( fabs(charge) != 1 ) {
    cout << "ERROR GEScaleSyst::GEScaleCorrPt invalid muon charge= " << charge << endl;
    return -1e9;
  }

  float kappa = GetGEScaleKappa(eta, phi, mode, icopy);

  if(verbose>=1) {
    cout << "\n Start GE pT correction:" << endl;
    cout << "    initial pT= " << pt << ", eta= " << eta << ", phi= " << phi << ", q= " << charge << endl;
    cout << "    matrix for mode= " << mode << ", icopy= " << icopy << ", copy idx= " << GetCopyIndex(mode, icopy) << endl;
    for(int ieta=0; ieta<neta; ++ieta) {
      for(int iphi=0; iphi<nphi; ++iphi) {
        printf("%3.5f    ", _matrix[GetCopyIndex(mode, icopy)][ieta][iphi]);
      }
      cout << endl;
    }
    cout << "\t --> kappa= " << kappa << endl;
  }

  float kappaGeV = kappa * 1.0e-3; // 1/TeV -> 1/GeV

  float corr_pt = doOpp ? (pt / (1 - charge * kappaGeV * pt)) : (pt / (1 + charge * kappaGeV * pt));
  float corr_p = corr_pt * cosh(eta);

  if(verbose>=1)  cout << "\t --> corr pT= " << corr_pt << endl;
  if(verbose>=1)  cout << "\t -->  corr p= " << corr_p << endl;
  if(verbose>=1)  cout << "\t -->    diff= " << (corr_pt-pt) << endl;

  // if corrected pt is greater then 6.5 TeV or even becomes negative, return 6.5 TeV
  if( (corr_pt > 6500.) || (corr_pt <= 0) ) {
    corr_pt = 6500.;
    if(verbose>=1)  cout << "\t\t --> re-corr pT= " << corr_pt << endl;
    if(verbose>=1)  cout << "\t\t -->       diff= " << (corr_pt-pt) << endl;
  }

  return corr_pt;
}

TLorentzVector GEScaleSyst::GEScaleCorrLvec(float pt, float eta, float phi, int charge, int mode, int icopy, bool doOpp) {
  TLorentzVector Lvec;
  float corr_pt = GEScaleCorrPt(pt, eta, phi, charge, mode, icopy, doOpp);

  if(corr_pt < 0) {
    cout << "ERROR   GEScaleSyst::GEScaleCorrLvec   invalid corrected pt= " << corr_pt << endl;
    Lvec.SetPtEtaPhiM(1e99, 0, 0, 0);
    return Lvec;
  }

  Lvec.SetPtEtaPhiM(corr_pt, eta, phi, mu_mass);

  if(verbose>=1)  cout << "Final muon pT, eta, phi" << endl;
  if(verbose>=1)  cout << "pT=  " << Lvec.Pt() << endl;
  if(verbose>=1)  cout << "eta= " << Lvec.Eta() << endl;
  if(verbose>=1)  cout << "phi= " << Lvec.Phi() << endl;
  if(verbose>=1)  cout << "E=   " << Lvec.E() << endl;
  if(verbose>=1)  cout << "m=   " << Lvec.M() << endl;

  return Lvec;
}

void GEScaleSyst::Init() {
  if(verbose>=1)  cout << "Initiating era: " << era << endl;

  if(era == "2016_UL_HIPM") { //https://twiki.cern.ch/twiki/bin/view/CMS/MuonUL2016
    _matrix[81600][0][0] = -0.075;       _matrix[81600][0][1] = -0.049;       _matrix[81600][0][2] = -0.007;
    _matrix[81600][1][0] = -0.053;       _matrix[81600][1][1] = -0.100;       _matrix[81600][1][2] =  0.003;
    _matrix[81600][2][0] = -0.020;       _matrix[81600][2][1] =  0.072;       _matrix[81600][2][2] = -0.015;
    _matrix[81600][3][0] = -0.009;       _matrix[81600][3][1] =  0.057;       _matrix[81600][3][2] = -0.003;
    _matrix[81600][4][0] = -0.010;       _matrix[81600][4][1] = -0.036;       _matrix[81600][4][2] =  0.060;
    _matrix[81600][5][0] =  0.078;       _matrix[81600][5][1] = -0.073;       _matrix[81600][5][2] = -0.196;

    _matrix[81601][0][0] = -0.075+0.088; _matrix[81601][0][1] = -0.049+0.101; _matrix[81601][0][2] = -0.007+0.060;
    _matrix[81601][1][0] = -0.053+0.055; _matrix[81601][1][1] = -0.100+0.040; _matrix[81601][1][2] =  0.003+0.038;
    _matrix[81601][2][0] = -0.020+0.029; _matrix[81601][2][1] =  0.072+0.030; _matrix[81601][2][2] = -0.015+0.030;
    _matrix[81601][3][0] = -0.009+0.029; _matrix[81601][3][1] =  0.057+0.034; _matrix[81601][3][2] = -0.003+0.036;
    _matrix[81601][4][0] = -0.010+0.040; _matrix[81601][4][1] = -0.036+0.043; _matrix[81601][4][2] =  0.060+0.040;
    _matrix[81601][5][0] =  0.078+0.067; _matrix[81601][5][1] = -0.073+0.102; _matrix[81601][5][2] = -0.196+0.084;

    _matrix[81602][0][0] = -0.075-0.088; _matrix[81602][0][1] = -0.049-0.101; _matrix[81602][0][2] = -0.007-0.060;
    _matrix[81602][1][0] = -0.053-0.055; _matrix[81602][1][1] = -0.100-0.040; _matrix[81602][1][2] =  0.003-0.038;
    _matrix[81602][2][0] = -0.020-0.029; _matrix[81602][2][1] =  0.072-0.030; _matrix[81602][2][2] = -0.015-0.030;
    _matrix[81602][3][0] = -0.009-0.029; _matrix[81602][3][1] =  0.057-0.034; _matrix[81602][3][2] = -0.003-0.036;
    _matrix[81602][4][0] = -0.010-0.040; _matrix[81602][4][1] = -0.036-0.043; _matrix[81602][4][2] =  0.060-0.040;
    _matrix[81602][5][0] =  0.078-0.067; _matrix[81602][5][1] = -0.073-0.102; _matrix[81602][5][2] = -0.196-0.084;
  } else if(era == "2016_UL") {
    _matrix[91600][0][0] =  0.002;       _matrix[91600][0][1] = -0.005;       _matrix[91600][0][2] =  0.026;
    _matrix[91600][1][0] =  0.049;       _matrix[91600][1][1] = -0.064;       _matrix[91600][1][2] = -0.142;
    _matrix[91600][2][0] = -0.045;       _matrix[91600][2][1] =  0.023;       _matrix[91600][2][2] =  0.039;
    _matrix[91600][3][0] = -0.061;       _matrix[91600][3][1] = -0.004;       _matrix[91600][3][2] =  0.000;
    _matrix[91600][4][0] =  0.016;       _matrix[91600][4][1] =  0.041;       _matrix[91600][4][2] =  0.003;
    _matrix[91600][5][0] = -0.108;       _matrix[91600][5][1] = -0.023;       _matrix[91600][5][2] = -0.091;

    _matrix[91601][0][0] =  0.002+0.078; _matrix[91601][0][1] = -0.005+0.089; _matrix[91601][0][2] =  0.026+0.086;
    _matrix[91601][1][0] =  0.049+0.046; _matrix[91601][1][1] = -0.064+0.039; _matrix[91601][1][2] = -0.142+0.046;
    _matrix[91601][2][0] = -0.045+0.033; _matrix[91601][2][1] =  0.023+0.032; _matrix[91601][2][2] =  0.039+0.031;
    _matrix[91601][3][0] = -0.061+0.033; _matrix[91601][3][1] = -0.004+0.031; _matrix[91601][3][2] =  0.000+0.032;
    _matrix[91601][4][0] =  0.016+0.043; _matrix[91601][4][1] =  0.041+0.045; _matrix[91601][4][2] =  0.003+0.042;
    _matrix[91601][5][0] = -0.108+0.109; _matrix[91601][5][1] = -0.023+0.114; _matrix[91601][5][2] = -0.091+0.073;

    _matrix[91602][0][0] =  0.002+0.078; _matrix[91602][0][1] = -0.005+0.089; _matrix[91602][0][2] =  0.026+0.086;
    _matrix[91602][1][0] =  0.049+0.046; _matrix[91602][1][1] = -0.064+0.039; _matrix[91602][1][2] = -0.142+0.046;
    _matrix[91602][2][0] = -0.045+0.033; _matrix[91602][2][1] =  0.023+0.032; _matrix[91602][2][2] =  0.039+0.031;
    _matrix[91602][3][0] = -0.061+0.033; _matrix[91602][3][1] = -0.004+0.031; _matrix[91602][3][2] =  0.000+0.032;
    _matrix[91602][4][0] =  0.016+0.043; _matrix[91602][4][1] =  0.041+0.045; _matrix[91602][4][2] =  0.003+0.042;
    _matrix[91602][5][0] = -0.108+0.109; _matrix[91602][5][1] = -0.023+0.114; _matrix[91602][5][2] = -0.091+0.073;
  } else if(era == "2017_UL") { // https://twiki.cern.ch/twiki/bin/view/CMS/MuonUL2017#Momentum_Scale
    _matrix[91700][0][0] = -0.049;       _matrix[91700][0][1] = -0.048;       _matrix[91700][0][2] = -0.025;  
    _matrix[91700][1][0] = -0.032;       _matrix[91700][1][1] =  0.001;       _matrix[91700][1][2] = -0.019;  
    _matrix[91700][2][0] =  0.027;       _matrix[91700][2][1] =  0.039;       _matrix[91700][2][2] = -0.049;  
    _matrix[91700][3][0] = -0.026;       _matrix[91700][3][1] =  0.027;       _matrix[91700][3][2] =  0.018;  
    _matrix[91700][4][0] =  0.012;       _matrix[91700][4][1] = -0.007;       _matrix[91700][4][2] =  0.000;  
    _matrix[91700][5][0] = -0.151;       _matrix[91700][5][1] =  0.058;       _matrix[91700][5][2] =  0.004;  

    _matrix[91701][0][0] = -0.049+0.070; _matrix[91701][0][1] = -0.048+0.044; _matrix[91701][0][2] = -0.025+0.044;       
    _matrix[91701][1][0] = -0.032+0.029; _matrix[91701][1][1] =  0.001+0.031; _matrix[91701][1][2] = -0.019+0.030;       
    _matrix[91701][2][0] =  0.027+0.024; _matrix[91701][2][1] =  0.039+0.027; _matrix[91701][2][2] = -0.049+0.021;       
    _matrix[91701][3][0] = -0.026+0.021; _matrix[91701][3][1] =  0.027+0.022; _matrix[91701][3][2] =  0.018+0.025;       
    _matrix[91701][4][0] =  0.012+0.025; _matrix[91701][4][1] = -0.007+0.023; _matrix[91701][4][2] =  0.000+0.026;       
    _matrix[91701][5][0] = -0.151+0.087; _matrix[91701][5][1] =  0.058+0.057; _matrix[91701][5][2] =  0.004+0.063;

    _matrix[91702][0][0] = -0.049-0.070; _matrix[91702][0][1] = -0.048-0.044; _matrix[91702][0][2] = -0.025-0.044;
    _matrix[91702][1][0] = -0.032-0.029; _matrix[91702][1][1] =  0.001-0.031; _matrix[91702][1][2] = -0.019-0.030;
    _matrix[91702][2][0] =  0.027-0.024; _matrix[91702][2][1] =  0.039-0.027; _matrix[91702][2][2] = -0.049-0.021;
    _matrix[91702][3][0] = -0.026-0.021; _matrix[91702][3][1] =  0.027-0.022; _matrix[91702][3][2] =  0.018-0.025;
    _matrix[91702][4][0] =  0.012-0.025; _matrix[91702][4][1] = -0.007-0.023; _matrix[91702][4][2] =  0.000-0.026;
    _matrix[91702][5][0] = -0.151-0.087; _matrix[91702][5][1] =  0.058-0.057; _matrix[91702][5][2] =  0.004-0.063;
  } else if(era == "2018_UL") {
    _matrix[91800][0][0] = -0.089;       _matrix[91800][0][1] =  0.057;       _matrix[91800][0][2] =  0.149;
    _matrix[91800][1][0] = -0.005;       _matrix[91800][1][1] = -0.028;       _matrix[91800][1][2] =  0.021;
    _matrix[91800][2][0] = -0.014;       _matrix[91800][2][1] = -0.009;       _matrix[91800][2][2] = -0.028;
    _matrix[91800][3][0] =  0.020;       _matrix[91800][3][1] =  0.007;       _matrix[91800][3][2] =  0.000;
    _matrix[91800][4][0] =  0.042;       _matrix[91800][4][1] = -0.068;       _matrix[91800][4][2] = -0.059;
    _matrix[91800][5][0] = -0.012;       _matrix[91800][5][1] = -0.235;       _matrix[91800][5][2] =  0.025;

    _matrix[91801][0][0] = -0.089+0.050; _matrix[91801][0][1] =  0.057+0.052; _matrix[91801][0][2] =  0.149+0.057;
    _matrix[91801][1][0] = -0.005+0.025; _matrix[91801][1][1] = -0.028+0.025; _matrix[91801][1][2] =  0.021+0.021;
    _matrix[91801][2][0] = -0.014+0.018; _matrix[91801][2][1] = -0.009+0.017; _matrix[91801][2][2] = -0.028+0.017;
    _matrix[91801][3][0] =  0.020+0.017; _matrix[91801][3][1] =  0.007+0.019; _matrix[91801][3][2] =  0.000+0.020;
    _matrix[91801][4][0] =  0.042+0.031; _matrix[91801][4][1] = -0.068+0.027; _matrix[91801][4][2] = -0.059+0.027;
    _matrix[91801][5][0] = -0.012+0.044; _matrix[91801][5][1] = -0.235+0.054; _matrix[91801][5][2] =  0.025+0.039;

    _matrix[91802][0][0] = -0.089-0.050; _matrix[91802][0][1] =  0.057-0.052; _matrix[91802][0][2] =  0.149-0.057;
    _matrix[91802][1][0] = -0.005-0.025; _matrix[91802][1][1] = -0.028-0.025; _matrix[91802][1][2] =  0.021-0.021;
    _matrix[91802][2][0] = -0.014-0.018; _matrix[91802][2][1] = -0.009-0.017; _matrix[91802][2][2] = -0.028-0.017;
    _matrix[91802][3][0] =  0.020-0.017; _matrix[91802][3][1] =  0.007-0.019; _matrix[91802][3][2] =  0.000-0.020;
    _matrix[91802][4][0] =  0.042-0.031; _matrix[91802][4][1] = -0.068-0.027; _matrix[91802][4][2] = -0.059-0.027;
    _matrix[91802][5][0] = -0.012-0.044; _matrix[91802][5][1] = -0.235-0.054; _matrix[91802][5][2] =  0.025-0.039;
  }
}

GEScaleSyst::GEScaleSyst() {}

GEScaleSyst::GEScaleSyst(std::string era_) {
  era = era_;
  Init();
}

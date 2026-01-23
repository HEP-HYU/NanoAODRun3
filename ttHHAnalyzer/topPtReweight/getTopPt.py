from __future__ import print_function
from ROOT import *
import ROOT
import sys, os
from collections import OrderedDict
gROOT.SetBatch(True)
hostname = os.environ["HOSTNAME"]

years = {"2016pre" : "RunIISummer20UL16NanoAODAPVv9",
         "2016post": "RunIISummer20UL16NanoAODv9",
         "2017"    : "RunIISummer20UL17NanoAODv9",
         "2018"    : "RunIISummer20UL18NanoAODv9"}
dss = {"TTTo2L2Nu": "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8",
       "TT_LFV_TCMuTau_Scalar": "TT_LFV_TToCMuTau_Scalar_TuneCP5_13TeV-madgraph-pythia8",
       "TT_LFV_TCMuTau_Vector": "TT_LFV_TToCMuTau_Vector_TuneCP5_13TeV-madgraph-pythia8",
       "TT_LFV_TCMuTau_Tensor": "TT_LFV_TToCMuTau_Tensor_TuneCP5_13TeV-madgraph-pythia8",
       "TT_LFV_TUMuTau_Scalar": "TT_LFV_TToUMuTau_Scalar_TuneCP5_13TeV-madgraph-pythia8",
       "TT_LFV_TUMuTau_Vector": "TT_LFV_TToUMuTau_Vector_TuneCP5_13TeV-madgraph-pythia8",
       "TT_LFV_TUMuTau_Tensor": "TT_LFV_TToUMuTau_Tensor_TuneCP5_13TeV-madgraph-pythia8"}
basedir = "/data1/common/NanoAOD/v9_UL/mc"
datasets = OrderedDict()
filelists = OrderedDict()

for year, era in years.items():
    for ds, das in dss.items():

        datasets[year + '_' + ds] = []
        if not os.path.isdir(os.path.join(basedir, era, das, 'NANOAODSIM')): continue
        gt = os.listdir(os.path.join(basedir, era, das, 'NANOAODSIM')) #must be one GT
        dirs = os.listdir(os.path.join(basedir, era, das, 'NANOAODSIM', gt[0]))
        for dd in dirs:
            datasets[year + '_' + ds].append(os.path.join(basedir, era, das, 'NANOAODSIM', gt[0], dd))

    for name, dataset in datasets.items():
        filelists[name] = []
        for path in dataset:
            if not os.path.isdir(path): continue
            for tmps in os.listdir(path):
                filelists[name].append(os.path.join(path, tmps))

    for name, files in filelists.items():
      entries = 0
      h = TH1F('h', 'h', 40, 0, 2000)
      for file in files:
        try: f = TFile.Open(file,'READ')
        except: continue
        print("processing "+file.split("/")[-5] + " " + file.split("/")[-1], end="\n")
        tree = f.Get("Events")
        htmp1 = TH1F('htmp1', 'htmp1', 40, 0, 2000)
        tree.Draw("GenPart_pt>>htmp1", "abs(GenPart_pdgId) == 6 && GenPart_statusFlags & (1 << 13)")
        h.Add(htmp1, 1.0)
        entries += tree.GetEntries()
        f.Close()
        del htmp1

      tmp_string = '['
      for i in range(1,h.GetNbinsX()+1):
        tmp_string += str(h.GetBinContent(i)) + ', '
      tmp_string += ']\n'
      tmp_string += str(entries)
      with open(name+'.txt', 'w') as f: f.write(tmp_string)
      del h

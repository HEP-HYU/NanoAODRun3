#!/usr/bin/env python
import os, sys
from ROOT import *
from array import array
gROOT.SetBatch(True)
import numpy as np

basedir = ['/data1/common/NanoAOD/mc/RunIISummer20UL16NanoAODAPVv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_preVFP_v9-v1/',
           '/data1/common/NanoAOD/mc/RunIISummer20UL16NanoAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/',
           '/data1/common/NanoAOD/mc/RunIISummer20UL17NanoAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_mc2017_realistic_v8-v1/',
           '/data1/common/NanoAOD/mc/RunIISummer20UL18NanoAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v15_L1v1-v1/']

path = {}
for bd in basedir:
    tmp = []
    for (root, dirs, files) in os.walk(bd):
        if len(files) > 0:
            for ff in files:
                if not ff.endswith(".root"): continue
                tmp.append(os.path.join(root, ff))

    path[bd.split('/')[5]] = tmp


for year, files in path.items():

    h1 = TH1F("h1","h1", 300, 0, 300)
    h2 = TH1F("h2","h2", 300, 0, 300)
    h1.SetDirectory(0)
    h2.SetDirectory(0)

    for filepath in files:
        print(filepath.split('/')[-1]+'\r')
        data = TFile.Open(filepath)
        tree = data.Get('Events')

        h3 = TH1F("h3","h3", 300, 0, 300)
        h4 = TH1F("h4","h4", 300, 0, 300)

        hfull = tree.Draw("LHE_HT>>h3", "(genWeight/abs(genWeight))", "")
        hgood = tree.Draw("LHE_HT>>h4", "(genWeight/abs(genWeight))*(LHE_HT < 100)", "")

        h1.Add(h3)
        h2.Add(h4)
        data.Close()

    if h1.Integral() == 0 : ratio = 1.0
    else: ratio = float(h2.Integral())/h1.Integral()

    print(year, ratio)

    h1.Reset()
    h2.Reset()

#!/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 11:01:46 2018

@author: Suyong Choi (Department of Physics, Korea University suyong@korea.ac.kr)
@refactored for LFV analysis: Jiwon Park (jiwon.park@cern.ch)

This script applies nanoaod processing to one file
"""
import sys, os, re, argparse
import cppyy
import ROOT

if __name__=='__main__':

    parser = argparse.ArgumentParser(usage="%prog [options]")
    parser.add_argument("-I", "--infile",  dest="infile", type=str, default="", help="Input file name")
    parser.add_argument("-O", "--outputroot", dest="outputroot", type=str, default="", help="Output file in your working directory")
    parser.add_argument("-Y", "--year", dest="year", type=str, default="", help="Select 2016pre, 2016post, 2017, or 2018 runs")
    parser.add_argument("-C", "--ch",  dest="ch", type=str, default="", help="Select electron or muon")
    parser.add_argument("-S", "--syst", dest="syst", type=str, default="theory", help="Systematic: 'data' for Data, 'nosyst' for mc without uncertainties. Default is 'theory'. To run without theory unc for TT samples, put 'all'.")
    parser.add_argument("-D", "--dataset", dest="dataset", action="store", nargs="+", default=[], help="Put dataset folder name (eg. TTTo2L2Nu) to process specific one.")
    parser.add_argument("-F", "--dataOrMC", dest="dataOrMC", type=str, default="", help="data or mc flag, if you want to process data-only or mc-only")
    parser.add_argument("-M", "--mode", dest="mode", type=str, default="", help="Only for fake rate: lss, los, tss, tos")
    parser.add_argument("--ff", dest="ff", action="store_true", default=False, help="Apply tau fake factor for final selection")
    options = parser.parse_args()

    print("Input: " + options.infile + ", " + "Output: " + options.outputroot)
    outputroot = options.outputroot
    infile = options.infile
    year = options.year
    ch   = options.ch
    syst = options.syst
    mode = options.mode
    applytauFF = options.ff

    json = ""
    if 'data' in infile[5:]:
        if syst == "data":
            if "2016" in infile:
                json = os.path.join(workdir, "data/GoldenJSON/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt")
            elif "2017" in infile:
                json = os.path.join(workdir, "data/GoldenJSON/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt")
            elif "2018" in infile:
                json = os.path.join(workdir, "data/GoldenJSON/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt")
            elif "2022" in infile:
                json = os.path.join(workdir, "data/GoldenJSON/Cert_Collisions2022_355100_362760_Golden.json")
            elif "2023" in infile:
                json = os.path.join(workdir, "data/GoldenJSON/Cert_Collisions2023_366442_370790_Golden.json")
            elif "2024" in infile:
                json = os.path.join(workdir, "data/GoldenJSON/Cert_Collisions2024_378981_386951_Golden.json")

    print("Input: " + infile + ", " + "Output: " + outputroot + ", Syst: " + syst + ", Json: " + json + "\n")

    if not re.match('.*\.root', outputroot):
        print("Output file should be a root file! Quitting")
        exit(-1)

    # load compiled C++ library into ROOT/python
    cppyy.load_reflection_info("libnanoadrdframe.so")
    t = ROOT.TChain("Events")
    tmpf = ROOT.TFile.Open(infile)
    tmpt = tmpf.Get("Events")
    if tmpt.GetEntries() > 0:
        t.Add(infile)
    aproc = None
    if len(mode) > 1:
        aproc = ROOT.TauFakeFactorAnalyzer(t, outputroot, year, ch, syst, json, "", 1, mode)
    else:
        aproc = ROOT.TopLFVAnalyzer(t, outputroot, year, ch, syst, json, applytauFF, "", 1)
    aproc.setupAnalysis()
    aproc.run(False, "Events")

    # process input rootfiles to sum up all the counterhistograms
    counterhistogramsum = None
    LHEPdfWeightSumAll = None
    PSWeightSumAll = None
    ScaleWeightSumAll = None

    intf = ROOT.TFile(infile)
    counterhistogram = intf.Get("hcounter")
    counterhistogramsum = counterhistogram.Clone()
    counterhistogramsum.SetDirectory(0)

    if syst == "theory":
        LHEPdfWeightSum = intf.Get("LHEPdfWeightSum")
        PSWeightSum = intf.Get("PSWeightSum")
        ScaleWeightSum = intf.Get("ScaleWeightSum")

        LHEPdfWeightSumAll = LHEPdfWeightSum.Clone()
        LHEPdfWeightSumAll.SetDirectory(0)
        PSWeightSumAll = PSWeightSum.Clone()
        PSWeightSumAll.SetDirectory(0)
        ScaleWeightSumAll = ScaleWeightSum.Clone()
        ScaleWeightSumAll.SetDirectory(0)

    intf.Close()

    if counterhistogramsum != None:
        print("Updating with counter histogram")
        outf = ROOT.TFile(outputroot, "UPDATE")
        counterhistogramsum.Write()
        outf.Write("", ROOT.TObject.kOverwrite)
        outf.Close()
    else:
        print("counter histogram not found")

    if LHEPdfWeightSumAll != None:
        print("Updating with theory weight sum histograms")
        outf = ROOT.TFile(outputroot, "UPDATE")
        LHEPdfWeightSumAll.Write()
        PSWeightSumAll.Write()
        ScaleWeightSumAll.Write()
        outf.Write("", ROOT.TObject.kOverwrite)
        outf.Close()
    else:
        print("theory weight sum histograms not found")

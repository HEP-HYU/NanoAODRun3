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

    print("processonefile.py")
    print("Input: " + options.infile + ", " + "Output: " + options.outputroot)
    outputroot = options.outputroot
    infile = options.infile
    year = options.year
    ch   = options.ch
    syst = options.syst
    mode = options.mode
    applytauFF = options.ff

    workdir = os.getcwd()
    script_dir = os.path.dirname(os.path.abspath(__file__))
    default_json_dir = os.path.join(os.path.dirname(script_dir), "data/GoldenJSON")
    json_base = os.path.join(workdir, "data/GoldenJSON") if os.path.isdir(os.path.join(workdir, "data/GoldenJSON")) else default_json_dir

    golden_json_map = (
        ("2022", "Cert_Collisions2022_355100_362760_Golden.json"),
        ("2023", "Cert_Collisions2023_366442_370790_Golden.json"),
        ("2024", "Cert_Collisions2024_378981_386951_Golden.json"),
        ("2016", "Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"),
        ("2017", "Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt"),
        ("2018", "Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt"),
    )

    json = ""
    if syst == "data":
        # Check explicit --year before falling back to infile substring match
        for tag, fname in golden_json_map:
            if tag in year:
                json = os.path.join(json_base, fname)
                break
        if not json:
            for tag, fname in golden_json_map:
                if tag in infile:
                    json = os.path.join(json_base, fname)
                    break
        if not json:
            print(f"[ERROR] syst is data but cannot determine Golden JSON for year '{year}' and infile '{infile}'!")
            sys.exit(1)

        if not os.path.isfile(json):
            print(f"[ERROR] Golden JSON file does not exist: '{json}'")
            sys.exit(1)

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
    if mode == "btageff":
        aproc = ROOT.BTagEfficiencyAnalyzer(t, outputroot, year, ch, syst, json, False, "", 1)
    elif len(mode) > 1:
        aproc = ROOT.TauFakeFactorAnalyzer(t, outputroot, year, ch, syst, json, False, "", 1, mode)
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
        LHEPdfWeightSum  = intf.Get("LHEPdfWeightSum")
        PSWeightSum      = intf.Get("PSWeightSum")
        ScaleWeightSum   = intf.Get("ScaleWeightSum")

        if LHEPdfWeightSum is not None:
            LHEPdfWeightSumAll = LHEPdfWeightSum.Clone()
            LHEPdfWeightSumAll.SetDirectory(0)
        if PSWeightSum is not None:
            PSWeightSumAll = PSWeightSum.Clone()
            PSWeightSumAll.SetDirectory(0)
        if ScaleWeightSum is not None:
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

    if LHEPdfWeightSumAll is not None or PSWeightSumAll is not None or ScaleWeightSumAll is not None:
        print("Updating with theory weight sum histograms")
        outf = ROOT.TFile(outputroot, "UPDATE")
        if LHEPdfWeightSumAll is not None: LHEPdfWeightSumAll.Write()
        if PSWeightSumAll is not None: PSWeightSumAll.Write()
        if ScaleWeightSumAll is not None: ScaleWeightSumAll.Write()
        outf.Write("", ROOT.TObject.kOverwrite)
        outf.Close()
    else:
        print("theory weight sum histograms not found")

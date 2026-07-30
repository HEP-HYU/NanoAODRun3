import os
import sys
import ROOT
from ROOT import *
import numpy as np
import optparse
from optparse import OptionParser

# Set Runs
years = ['v12_2022', 'v12_2022EE', 'v12_2023', 'v12_2023BPix', 'v15_2024']

parser = OptionParser(usage="%prog [options]")
parser.add_option("-I", "--input",  dest="input", type="string", default="", help="Input file name")
(options, args) = parser.parse_args()
base_path = options.input

def collect_systhists(inpath, outfile, fname, year):
    fileList_wt_syst = [i for i in os.listdir(inpath) if fname in i]
    for sysfile in fileList_wt_syst:
        if ('SingleMuon' in sysfile or 'Muon' in sysfile) and sysfile not in ['hist_SingleMuon.root', 'hist_Muon.root']: continue
        tmpf = TFile.Open(os.path.join(inpath, sysfile), 'READ')
        hlists = [ h.GetName() for h in tmpf.GetListOfKeys() \
                   if any(i in h.GetName() for i in ['dnn_pred_S5', 'counter', 'muon1_pt', 'tau1_pt', 'jet1_pt', 'jet2_pt', 'mutau_mass', 'mutau_dr']) ]
        for histname in hlists:
            tmpf.cd()
            tmphist = tmpf.Get(histname)
            tmphist.SetDirectory(0)
            outfile.cd()
            if not "__" in sysfile:
                tmphist.Write()
            if "__" in sysfile:
                syst = sysfile.split("__")[1].split(".root")[0] 
                tmphist.Write(histname+"__"+syst)
        tmpf.Close()
    return

for year in years:

    inpath  = os.path.join(base_path, year + '_postprocess')
    outpath = os.path.join(base_path, year + '_postprocess_2')
    print("Directory created: ", outpath)
    os.makedirs(outpath, exist_ok=True)

    file_list = [i.replace('.root','') for i in os.listdir(inpath) if '.root' in i and not '__' in i]
    for fname in file_list:
        outfname = fname + ".root"
        print("Saving histograms at {}/{}".format(outpath, outfname))
        outfile = TFile.Open(os.path.join(outpath, outfname), 'RECREATE')
        collect_systhists(inpath, outfile, fname, year)
        outfile.Close()

    inpath  = os.path.join(base_path, year + '_postprocess', 'fake')
    if os.path.exists(inpath):
        file_list = [i.replace('.root','') for i in os.listdir(inpath) if '.root' in i and not '__' in i]
        for fname in file_list:
            fname_tmp = fname.replace('hist_', '')
            outfname = 'hist_fake_' + fname_tmp + ".root"
            print("Saving histograms at {}/{}".format(outpath, outfname))
            outfile = TFile.Open(os.path.join(outpath, outfname), 'RECREATE')
            collect_systhists(inpath, outfile, fname, year)
            outfile.Close()

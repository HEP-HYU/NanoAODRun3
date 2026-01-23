import os, sys, glob, re
import ROOT
from ROOT import *
import numpy as np
import subprocess
import optparse
import array
import re

from optparse import OptionParser
parser = OptionParser(usage="%prog [options]")
parser.add_option("-I", "--infile", dest="infile", type="string", default="", help="Input file name")
parser.add_option("-Y", "--year", dest="year", type="string", default="", help="Select v2022, v2022EE, v2023, or v2023_BPix  for years")
parser.add_option("--postfix", dest="postfix", type="string", default="", help="Add postfix to output here, to have rebinning for histograms")
parser.add_option("-F", "--forceHadd", dest="forceHadd", action="store_true", default=False, help="Force hadd split files")
parser.add_option("-N", "--noHadd", dest="noHadd", action="store_true", default=False, help="Skip hadd split files")
(options, args) = parser.parse_args()

year = options.year
input = options.infile
forceHadd = options.forceHadd

# starting bin -> 0.01, trick for logX
#rebin_arr = array.array('d', [0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 1.0, 2.0, 5.0, 10.0, 30, 100.0])
rebin_arr = array.array('d', [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.16, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 5.0, 10.0, 30, 100.0])
rebin_pt = array.array('d', [0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 600])
rebin_jetpt = array.array('d', [0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 400])
rebin_mass = array.array('d', [0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500, 1000])

# Set to pre-approval binning if postfix is set, for histogramming
if len(options.postfix) > 0:
    rebin_arr = array.array('d', [0.01, 1.0, 2.0, 5.0, 10.0, 30, 100.0])

#if year not in ['2016pre', '2016post', '2017', '2018']:
if year not in ['v2022', 'v2022EE', 'v2023', 'v2023_BPix', 'v2024']:
    print('Wrong year, check again')
    sys.exit()
if forceHadd: print("Hadd all split MC!!")

yield_name = 'h_ncleanjetspass'
base_path = './'
nom_path = os.path.join(base_path, input, year)
if not os.path.exists(nom_path):
    print("Folder '{}' does not exists.".format(nom_path))
    sys.exit()
else:
    print("Start postprocessing at '{}'.".format(nom_path))

# Set output folders
out_path = os.path.join(base_path, input, year + '_rebin' + options.postfix)
fig_path = os.path.join(base_path, input, '../fig_mu/rebin_' + year + options.postfix)
if not os.path.exists(out_path):
    os.makedirs(out_path)
if not os.path.exists(fig_path):
    os.makedirs(fig_path)
#if not os.path.exists(os.path.join(fig_path, 'qcd')):
#    os.makedirs(os.path.join(fig_path, 'qcd'))
#if not os.path.exists(os.path.join(fig_path, 'dyincl')):
#    os.makedirs(os.path.join(fig_path, 'dyincl'))
#if not os.path.exists(os.path.join(fig_path, 'dyincl', 'qcd')):
#    os.makedirs(os.path.join(fig_path, 'dyincl', 'qcd'))

file_list = [i.replace('.root', '') for i in os.listdir(nom_path) if '.root' in i]
data_list = [i[:i.find('202')] for i in os.listdir(nom_path) if '.root' in i and '202' in i and 'jes' not in i]
data_list = list(set(data_list))

print(data_list)
print(file_list)

# Loop over all files.
for fname in file_list:
    #print(os.path.join(nom_path, fname))
    #if not any(i in fname for i in ['TTTo2L2Nu', 'TTToSemiLeptonic']): continue
    #if not any(i in fname for i in ['hdamp', 'tune']): continue

    #flag for ext. syst with different normalization
    run_on_syst = False

    infile = TFile.Open(os.path.join(nom_path, fname + '.root'), 'READ')
    hlists = [ h.GetName() for h in infile.GetListOfKeys() if '_S' in h.GetName() ]
    #print (hlists)
    hlists.append("hcounter")

    # Collecting Histograms in outfile.
    print("Saving histograms at {}/{}.root".format(out_path, fname))
    outfile = TFile.Open(os.path.join(out_path, fname+'.root'), 'RECREATE')

    nominal_list = []

    do_renorm = True
    if 'LFV' in fname: do_renorm = False

    for hname in hlists:
        if "__" not in hname: nominal_list.append(hname)
        if run_on_syst: continue
        h = infile.Get(hname)
        if 'dnn_pred' in hname:
            h.SetBinContent(2, h.GetBinContent(1) + h.GetBinContent(2))
            h.SetBinError(2, sqrt(pow(h.GetBinError(1), 2) + pow(h.GetBinError(2), 2)))
            h.SetBinContent(1, 0.)
            h.SetBinError(1, 0.)
            h = h.Rebin(len(rebin_arr)-1, h.GetName(), rebin_arr)
        if ('muon1_pt' in hname) or ('electron1_pt' in hname):
            h_rebinned = h.Rebin(len(rebin_pt) - 1, hname+"_rebinned", rebin_pt)
            h.Write()
            h_rebinned.Write()
        if ('tau1_pt' in hname) or ('jet1_pt' in hname):
            h_rebinned = h.Rebin(len(rebin_jetpt) - 1, hname+"_rebinned", rebin_jetpt)
            h.Write()
            h_rebinned.Write()
        if ('leptau_mass' in hname):
            h_rebinned = h.Rebin(len(rebin_mass) - 1, hname+"_rebinned", rebin_mass)
            h.Write()
            h_rebinned.Write()
        if yield_name in hname:
            h1 = h.Clone('h1')
            h1.SetName(hname.replace(yield_name, yield_name + '_yield'))
            if "__" not in h1.GetName():
                nominal_list.append(h1.GetName())
            h1.Write()
        h.Write()

    hcounter = infile.Get('hcounter')

    infile.Close()
    outfile.Close()


for dataname in data_list:
    try:
        subprocess.call(['rm', os.path.join(out_path, dataname + '.root')])
    except: pass
    subprocess.check_call( ["hadd", "-f", os.path.join(out_path, dataname + '.root')] + glob.glob(os.path.join(out_path, dataname) + '201*') )

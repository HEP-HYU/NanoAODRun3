from ROOT import TFile
import ROOT
import sys, os
from subprocess import check_call
from os import listdir, path
import collections
import glob
import multiprocessing
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-I', '--input', dest='input', type=str, default="test")
parser.add_argument("--postfix", dest="postfix", type=str, default="", help="Add postfix to output here, to have rebinning for histograms")
args = parser.parse_args()
input = args.input
channel = "electron"
if "muon" in args.input: channel = "muon"

lumi_dict = {'2016pre': 19502, '2016post': 16812, '2017': 41480, '2018':59832}#16: 36314, run2:137625
lumi_dict = {'2022': 7980.4, '2022EE': 26671.7, '2023': 17794.0, '2023BPix': 9451.0, '2024':109000.0}
file_names = collections.OrderedDict()

if not os.path.exists(os.path.join(input, 'Run3' + args.postfix)):
    try: os.makedirs(os.path.join(input, 'Run3' + args.postfix))
    except: pass

def store_file(it):
    path = it[0]
    file_name = it[1]
    for f in file_name:
        print(os.path.join(path, f[:f.rfind('_')] + '.root'))
        ftmp = TFile.Open(os.path.join(path, f[:f.rfind('_')] + '.root'), 'READ')

        hcounter = ftmp.Get("hcounter")

        hist_names = [x.GetName() for x in ftmp.GetListOfKeys()]
        hist_names = list(dict.fromkeys(hist_names)) #remove duplicates from more than one instances 
        hist_names[:] = [item for item in hist_names if item not in ['hcounter']]
        hist_names.sort()

        ntmp = ftmp.Get("hcounter").GetBinContent(2)

        dest_name = os.path.join(input, 'Run3' + args.postfix, f)
        print("destination :", dest_name)
        dest = TFile.Open(dest_name, 'RECREATE')
        dest.cd()
        hcounter.Write()

        print('Writing scaled histogram to ' + dest_name)

        for hist in hist_names:
            htmp = ftmp.Get(hist)

            if "hcounter" in hist:
                htmp.Write()
            else:
                htmp.Scale((lumi_dict[path.split('/')[-1].split('_')[1]]/170897.1)/ntmp)
                dest.cd()
                htmp.Write()
        dest.Write()
        dest.Close()


if __name__ == '__main__':

    for era in ['2022', '2022EE', '2023', '2023BPix', '2024']:
        dir_path = os.path.join(input, "v12_"+era+'_postprocess' + args.postfix)
        if era == "2024": dir_path = os.path.join(input, "v15_"+era+'_postprocess' + args.postfix)
        dirs = os.listdir(dir_path)
        print("POST process path: " , dir_path)
        dirs[:] = [item.replace('.root', "_"+era+'.root') for item in dirs if any(i in item for i in ['LFV']) if '__' not in item]
        print("EDITED DIRS: " , dirs)
        file_names[dir_path] = dirs

    pool = multiprocessing.Pool(12)
    pool.map(store_file, file_names.items())
    pool.close()
    pool.join()

    chs = []
    if channel == "muon":
        chs = ['TCMuTau-LFV-Scalar', 'TCMuTau-LFV-Tensor', 'TCMuTau-LFV-Vector',
               'TUMuTau-LFV-Scalar', 'TUMuTau-LFV-Tensor', 'TUMuTau-LFV-Vector',
               'TTtoCMuTau-LFV-Scalar', 'TTtoCMuTau-LFV-Tensor', 'TTtoCMuTau-LFV-Vector',
               'TTtoUMuTau-LFV-Scalar', 'TTtoUMuTau-LFV-Tensor', 'TTtoUMuTau-LFV-Vector']
    else:
        chs = ['TCETau-LFV-Scalar', 'TCETau-LFV-Tensor', 'TCETau-LFV-Vector',
               'TUETau-LFV-Scalar', 'TUETau-LFV-Tensor', 'TUETau-LFV-Vector',
               'TTtoCETau-LFV-Scalar', 'TTtoCETau-LFV-Tensor', 'TTtoCETau-LFV-Vector',
               'TTtoUETau-LFV-Scalar', 'TTtoUETau-LFV-Tensor', 'TTtoUETau-LFV-Vector']

    for ch in chs:
        print(os.path.join(input, 'Run3' + args.postfix, 'hist_' + ch + '.root'))
        check_call(['hadd','-f', os.path.join(input, 'Run3' + args.postfix, 'hist_' + ch + '.root')] +  glob.glob(os.path.join(input, 'Run3' + args.postfix, 'hist_' + ch + '_202*.root')))

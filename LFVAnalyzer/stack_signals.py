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
parser.add_argument('-C', '--channel', dest='channel', type=str, default="")
args = parser.parse_args()
input = args.input
channel = args.channel

lumi_dict = {'2016pre': 19502, '2016post': 16812, '2017': 41480, '2018':59832}#16: 36314, run2:137625
lumi_dict = {'2022': 7980.4, '2022EE': 26671.7, '2023': 17794.0, '2023BPix': 9451.0, '2024':109000.0}
file_names = collections.OrderedDict()

syst = ["","jesAbsoluteup","jesAbsolutedown", "jesAbsolute_ERAup", "jesAbsolute_ERAdown",
        "jesBBEC1up", "jesBBEC1down", "jesBBEC1_ERAup", "jesBBEC1_ERAdown",
        "jesFlavorQCDup", "jesFlavorQCDdown", "jesRelativeBalup", "jesRelativeBaldown",
        "jesRelativeSample_ERAup", "jesRelativeSample_ERAdown",
        "jerup", "jerdown", "tesup", "tesdown",
        "hdampup", "hdampdown", "tuneup", "tunedown"]
syst2 = []

for sy in syst:
    if 'ERA' in sy:
        for key in lumi_dict.keys():
            syst2.append('__' + sy.replace('ERA', key[:4]))
    elif sy == "": syst2.append(sy)
    else:          syst2.append('__' + sy)

syst2 = list(set(syst2))

if not os.path.exists(input + "/Run3"):
    try: os.makedirs(input + "/Run3")
    except: pass

def store_file(it):
    path = it[0]
    file_name = it[1]
    for f in file_name:
        print(os.path.join(path, f[:f.rfind('_')] + '.root'))
        ftmp = TFile.Open(os.path.join(path, f[:f.rfind('_')] + '.root'), 'READ')

        hist_names = [x.GetName() for x in ftmp.GetListOfKeys()]
        hist_names = list(dict.fromkeys(hist_names)) #remove duplicates from more than one instances 
        hist_names[:] = [item for item in hist_names if item not in ['hcounter']]
        hist_names.sort()

        ntmp = ftmp.Get("hcounter").GetBinContent(2)

        dest_name = path.split('/')[0] + '/Run3/' + f
        dest = TFile.Open(dest_name, 'RECREATE')

        print('Writing scaled histogram to ' + dest_name)

        for hist in hist_names:
            htmp = ftmp.Get(hist)

            if "hcounter" in hist:
                htmp.Write()
            else:
                htmp.Scale(lumi_dict[path.split('/')[1].split('_')[0]]/ntmp)
                dest.cd()
                htmp.Write()
        dest.Write()
        dest.Close()


if __name__ == '__main__':

    for era in ['v12_2022', 'v12_2022EE', 'v12_2023', 'v12_2023BPix', 'v15_2024']:
        dir_path = os.path.join(input, era+'_postprocess')
        dirs = os.listdir(dir_path)
        print("POST process path: " , dir_path)
        dirs[:] = [item.replace('.root', '_' + era + '.root') for item in dirs if any(i in item for i in ['LFV'])]
        print("EDITED DIRS: " , dirs)
        file_names[dir_path] = dirs

    pool = multiprocessing.Pool(20)
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
        for syst in syst2:
            print(input + '/Run3/hist_' + ch + syst + '_202*.root')
            check_call(['hadd','-f', input + '/Run3/hist_' + ch + syst + '.root'] +  glob.glob(input + '/Run3/hist_' + ch + syst + '_202*.root'))

import warnings
warnings.filterwarnings("ignore")

import os, sys
import uproot
import pandas as pd
import numpy as np

f_dir = sys.argv[1]

#inputvars = ["btagWeight_DeepFlavB_perJet", 'Jet_pt', 'event'] #/data1/common/skimmed_NanoAOD/skim_test2/mc/2016pre/WJetsToLNu_HT1200To2500/280000_0C76D3EA-6565-B24B-B716-B9BB61B55D59.root
inputvars = ["btagWeight_DeepFlavB_perJet", "btagWeight_DeepFlavB_jes_perJet", 'Jet_pt', 'event'] #/data1/common/skimmed_NanoAOD/skim_test2/mc/2016pre/WJetsToLNu_HT1200To2500/120000_6977F95D-BF50-DC4D-9CD6-A627B125DD39.root

infile = uproot.open(f_dir)
tree = infile["Events"]
pd_data = tree.arrays(inputvars,library="pd")
#pd_data = pd_data[pd_data['event']==165452]
pd_data = pd_data[pd_data['event']==734551]
with pd.option_context('display.max_rows', None, 'display.max_columns', None, 'display.max_colwidth', None):
    print(pd_data)

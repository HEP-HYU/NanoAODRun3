import sys, os
import re
from subprocess import call
import argparse

parser = argparse.ArgumentParser(description="Run plotIt for Run 3")
parser.add_argument('-I', '--input', dest='input', type=str, default="process_0715_v4/electron", help="Input directory path")
parser.add_argument('-C', '-ch', '--channel', dest='channel', type=str, default="", choices=["", "electron", "muon"], help="Channel: electron or muon (auto-detected from input if empty)")
parser.add_argument("-D", dest="DNN", action="store_true", default=False, help="Run for DNN histograms")
parser.add_argument("-y", dest="yield_only", action="store_true", default=False, help="Run for yield histograms")
parser.add_argument("--postfix", dest="postfix", type=str, default="", help="Add postfix to output here, to have rebinning for histograms")
args = parser.parse_args()

dest_path = args.input
postfix = args.postfix

# Auto-detect channel if not specified
channel = args.channel
if not channel:
    if 'muon' in dest_path.lower():
        channel = 'muon'
    else:
        channel = 'electron'

config_path = f'../plotIt/configs/Run3_{channel}/'

syst_avoid = ['muonhighscale', 'metUnclust']

eras = {
    '2022':     {'lumi': 7980.4,   'config_tag': '2022',      'dir_candidates': ['v12_2022', 'v2022', '2022']},
    '2022EE':   {'lumi': 26671.7,  'config_tag': '2022EE',    'dir_candidates': ['v12_2022EE', 'v2022EE', '2022EE']},
    '2023':     {'lumi': 17794.0,  'config_tag': '2023',      'dir_candidates': ['v12_2023', 'v2023', '2023']},
    '2023BPix': {'lumi': 9451.0,   'config_tag': '2023_BPix', 'dir_candidates': ['v12_2023BPix', 'v12_2023_BPix', 'v2023_BPix', 'v2023BPix', '2023BPix', '2023_BPix']},
    '2024':     {'lumi': 109000.0, 'config_tag': '2024',      'dir_candidates': ['v15_2024', 'v2024', '2024']}
}

total_lumi = sum(info['lumi'] for info in eras.values()) # 170897.1 pb^-1

# Collect common systematics across era configs
common_syst_list = []
for era_key, era_info in eras.items():
    cfg_file = os.path.join(config_path, f"config_{era_info['config_tag']}.yml")
    if os.path.exists(cfg_file):
        with open(cfg_file) as f:
            lines = f.readlines()
            for line in lines:
                if any(s in line for s in syst_avoid): continue
                if "  - " in line:
                    common_syst_list.append(line)

common_syst_list = sorted(set(common_syst_list), key=common_syst_list.index)
common_syst = 'systematics:\n'
for item in common_syst_list:
    common_syst += item

# Create output figure folder
fig_dir = os.path.join(dest_path, 'figure_run3' + postfix)
if not os.path.exists(fig_dir):
    try: os.makedirs(fig_dir)
    except: pass

# Determine per-era folder structure inside dest_path
resolved_eras = {}
for era_key, era_info in eras.items():
    matched_dir = None
    for cand in era_info['dir_candidates']:
        test_path1 = os.path.join(dest_path, f"{cand}_postprocess{postfix}")
        test_path2 = os.path.join(dest_path, f"{cand}{postfix}")
        test_path3 = os.path.join(dest_path, cand)
        if os.path.exists(test_path1):
            matched_dir = f"{cand}_postprocess{postfix}"
            break
        elif os.path.exists(test_path2):
            matched_dir = f"{cand}{postfix}"
            break
        elif os.path.exists(test_path3):
            matched_dir = cand
            break

    if not matched_dir:
        matched_dir = f"{era_info['dir_candidates'][0]}_postprocess{postfix}"
    resolved_eras[era_key] = matched_dir

has_merged_run3 = os.path.exists(os.path.join(dest_path, 'Run3' + postfix))

string_for_files = ''
for era_key, era_info in eras.items():
    lumi = era_info['lumi']
    era_dir = resolved_eras[era_key]
    files_yml_path = os.path.join(config_path, 'files.yml')

    if os.path.exists(files_yml_path):
        with open(files_yml_path) as f:
            lines = f.readlines()
            skip_signal = False
            for line in lines:
                if skip_signal and 'hist' in line: skip_signal = False
                if 'LFV' in line and has_merged_run3: skip_signal = True
                if line.startswith('#'): skip_signal = True

                if 'hist' in line and not line.startswith('#'):
                    line = line[0] + os.path.join(dest_path, era_dir, line[1:]).lstrip('/')
                    if '\':\n' in line or '\':' in line:
                        hist_file_name = line.split('/')[-1]
                        line = line[0] + os.path.normpath(os.path.join(dest_path, era_dir, hist_file_name)) + '\':\n'

                    if not any(i in line for i in ['LFV', 'SingleMuon', 'Egamma', 'EGamma', 'Muon', 'WJetsToLNu_HT0To100']):
                        line += '  scale: ' + str(float(lumi) / total_lumi) + '\n'
                    elif 'WJetsToLNu_HT0To100' in line:
                        line += '  scale: ' + str(1.0288 * float(lumi) / total_lumi) + '\n'

                if 'group' in line:
                    if 'yields-group' in line:
                        if any(i in line for i in ['ttbar LL', 'ttbar LJ']):
                            line = "  yields-group: '1\\ttbar' \n"
                        elif any(i in line for i in ['ttbar Had', 'tt+X', 'VV', 'Z+Jets', 'W+Jets']):
                            line = "  yields-group: '9Other' \n"
                    else:
                        if any(i in line for i in ['Gttlj', 'Gttll']):
                            line = '  group: Gtt \n'
                        elif any(i in line for i in ['Gttjj', 'GttV', 'GVV', 'GZJets', 'GWJets']):
                            line = '  group: Gother \n'

                if not skip_signal:
                    string_for_files += line

# Write files_Run3.yml
files_run3_yml = os.path.join(config_path, 'files_Run3.yml')
with open(files_run3_yml, 'w+') as fnew:
    if has_merged_run3:
        if channel == 'electron':
            fnew.write(f"""
'{dest_path}/Run3{postfix}/hist_TUETau-LFV-Scalar.root':
  type: signal
  pretty-name: 'LFVSTus'
  cross-section: 0.076313
  generated-events: 1
  group: GLFVSTus
  yields-group: '1ST LFV \\cPqt\\cPqe\\tau Scalar'
  order: 1

'{dest_path}/Run3{postfix}/hist_TUETau-LFV-Vector.root':
  type: signal
  pretty-name: 'LFVSTuv'
  cross-section: 0.352889
  generated-events: 1
  group: GLFVSTuv
  yields-group: '1ST LFV \\cPqt\\cPqe\\tau Vector'
  order: 1

'{dest_path}/Run3{postfix}/hist_TUETau-LFV-Tensor.root':
  type: signal
  pretty-name: 'LFVSTut'
  cross-section: 1.61537
  generated-events: 1
  group: GLFVSTut
  yields-group: '2ST LFV \\cPqt\\cPqe\\tau Tensor'
  order: 1

'{dest_path}/Run3{postfix}/hist_TCETau-LFV-Scalar.root':
  type: signal
  pretty-name: 'LFVSTcs'
  cross-section: 0.00488095
  generated-events: 1
  group: GLFVSTcs
  yields-group: '3ST LFV \\cPqt\\cPqc\\tau Scalar'
  order: 1

'{dest_path}/Run3{postfix}/hist_TCETau-LFV-Vector.root':
  type: signal
  pretty-name: 'LFVSTcv'
  cross-section: 0.0253796
  generated-events: 1
  group: GLFVSTcv
  yields-group: '3ST LFV \\cPqt\\cPqc\\tau Vector'
  order: 3

'{dest_path}/Run3{postfix}/hist_TCETau-LFV-Tensor.root':
  type: signal
  pretty-name: 'LFVSTct'
  cross-section: 0.125558
  generated-events: 1
  group: GLFVSTct
  yields-group: '4ST LFV \\cPqt\\cPqc\\tau Tensor'
  order: 1

'{dest_path}/Run3{postfix}/hist_TTtoUETau-LFV-Scalar.root':
  type: signal
  pretty-name: 'LFVTTus'
  cross-section: 0.00134402
  generated-events: 1
  group: GLFVTTus
  yields-group: '5TT LFV \\cPqt\\cPqu\\e\\tau Scalar'
  order: 1

'{dest_path}/Run3{postfix}/hist_TTtoUETau-LFV-Vector.root':
  type: signal
  pretty-name: 'LFVTTuv'
  cross-section: 0.0106061
  generated-events: 1
  group: GLFVTTuv
  yields-group: '5TT LFV \\cPqt\\cPqu\\e\\tau Vector'
  order: 2

'{dest_path}/Run3{postfix}/hist_TTtoUETau-LFV-Tensor.root':
  type: signal
  pretty-name: 'LFVTTut'
  cross-section: 0.0646215
  generated-events: 1
  group: GLFVTTut
  yields-group: '6TT LFV \\cPqt\\cPqu\\e\\tau Tensor'
  order: 1

'{dest_path}/Run3{postfix}/hist_TTtoCETau-LFV-Scalar.root':
  type: signal
  pretty-name: 'LFVTTcs'
  cross-section: 0.00134201
  generated-events: 1
  group: GLFVTTcs
  yields-group: '7TT LFV \\cPqt\\cPqc\\e\\tau Scalar'
  order: 1

'{dest_path}/Run3{postfix}/hist_TTtoCETau-LFV-Vector.root':
  type: signal
  pretty-name: 'LFVTTcv'
  cross-section: 0.0107468
  generated-events: 1
  group: GLFVTTcv
  yields-group: '7TT LFV \\cPqt\\cPqc\\e\\tau Vector'
  order: 4

'{dest_path}/Run3{postfix}/hist_TTtoCETau-LFV-Tensor.root':
  type: signal
  pretty-name: 'LFVTTct'
  cross-section: 0.0645612
  generated-events: 1
  group: GLFVTTct
  yields-group: '8TT LFV \\cPqt\\cPqc\\e\\tau Tensor'
  order: 1
""")
        else: # muon
            fnew.write(f"""
'{dest_path}/Run3{postfix}/hist_TUMuTau-LFV-Scalar.root':
  type: signal
  pretty-name: 'LFVSTus'
  cross-section: 0.076246
  generated-events: 1
  group: GLFVSTus
  yields-group: '1ST LFV \\cPqt\\cPqu\\mu\\tau Scalar'
  order: 1

'{dest_path}/Run3{postfix}/hist_TUMuTau-LFV-Vector.root':
  type: signal
  pretty-name: 'LFVSTuv'
  cross-section: 0.352822
  generated-events: 1
  group: GLFVSTuv
  yields-group: '1ST LFV \\cPqt\\cPqu\\mu\\tau Vector'
  order: 1

'{dest_path}/Run3{postfix}/hist_TUMuTau-LFV-Tensor.root':
  type: signal
  pretty-name: 'LFVSTut'
  cross-section: 1.60867
  generated-events: 1
  group: GLFVSTut
  yields-group: '2ST LFV \\cPqt\\cPqu\\mu\\tau Tensor'
  order: 1

'{dest_path}/Run3{postfix}/hist_TCMuTau-LFV-Scalar.root':
  type: signal
  pretty-name: 'LFVSTcs'
  cross-section: 0.00488296
  generated-events: 1
  group: GLFVSTcs
  yields-group: '3ST LFV \\cPqt\\cPqc\\mu\\tau Scalar'
  order: 1

'{dest_path}/Run3{postfix}/hist_TCMuTau-LFV-Vector.root':
  type: signal
  pretty-name: 'LFVSTcv'
  cross-section: 0.0251451
  generated-events: 1
  group: GLFVSTcv
  yields-group: '3ST LFV \\cPqt\\cPqc\\mu\\tau Vector'
  order: 3

'{dest_path}/Run3{postfix}/hist_TCMuTau-LFV-Tensor.root':
  type: signal
  pretty-name: 'LFVSTct'
  cross-section: 0.123883
  generated-events: 1
  group: GLFVSTct
  yields-group: '4ST LFV \\cPqt\\cPqc\\mu\\tau Tensor'
  order: 1

'{dest_path}/Run3{postfix}/hist_TTtoUMuTau-LFV-Scalar.root':
  type: signal
  pretty-name: 'LFVTTus'
  cross-section: 0.00298052
  generated-events: 1
  group: GLFVTTus
  yields-group: '5TT LFV \\cPqt\\cPqu\\mu\\tau Scalar'
  order: 1

'{dest_path}/Run3{postfix}/hist_TTtoUMuTau-LFV-Vector.root':
  type: signal
  pretty-name: 'LFVTTuv'
  cross-section: 0.023822
  generated-events: 1
  group: GLFVTTuv
  yields-group: '5TT LFV \\cPqt\\cPqu\\mu\\tau Vector'
  order: 2

'{dest_path}/Run3{postfix}/hist_TTtoUMuTau-LFV-Tensor.root':
  type: signal
  pretty-name: 'LFVTTut'
  cross-section: 0.142932
  generated-events: 1
  group: GLFVTTut
  yields-group: '6TT LFV \\cPqt\\cPqu\\mu\\tau Tensor'
  order: 1

'{dest_path}/Run3{postfix}/hist_TTtoCMuTau-LFV-Scalar.root':
  type: signal
  pretty-name: 'LFVTTcs'
  cross-section: 0.00297597
  generated-events: 1
  group: GLFVTTcs
  yields-group: '7TT LFV \\cPqt\\cPqc\\mu\\tau Scalar'
  order: 1

'{dest_path}/Run3{postfix}/hist_TTtoCMuTau-LFV-Vector.root':
  type: signal
  pretty-name: 'LFVTTcv'
  cross-section: 0.0241416
  generated-events: 1
  group: GLFVTTcv
  yields-group: '7TT LFV \\cPqt\\cPqc\\mu\\tau Vector'
  order: 4

'{dest_path}/Run3{postfix}/hist_TTtoCMuTau-LFV-Tensor.root':
  type: signal
  pretty-name: 'LFVTTct'
  cross-section: 0.142798
  generated-events: 1
  group: GLFVTTct
  yields-group: '8TT LFV \\cPqt\\cPqc\\mu\\tau Tensor'
  order: 1
""")
    fnew.write(string_for_files)

# Build config_Run3.yml from template_Run3.yml
template_yml_path = os.path.join(config_path, 'template_Run3.yml')
config_run3_yml = os.path.join(config_path, 'config_Run3.yml')

if os.path.exists(template_yml_path):
    with open(template_yml_path) as f:
        lines = f.readlines()
    with open(config_run3_yml, 'w+') as f1:
        for line in lines:
            if 'luminosity:' in line and not line.strip().startswith('#'):
                f1.write(f"  luminosity: {total_lumi:.1f}\n")
            else:
                f1.write(line)
        f1.write(common_syst)

        if 'FF' in dest_path:
            if args.yield_only:
                f1.write("\nplots:\n  include: ['histos_yield_S5.yml']\n")
            else:
                f1.write("\nplots:\n  include: ['histos_FFapply.yml', 'histos_yield_S5.yml']\n")
        elif args.DNN:
            f1.write("\nplots:\n  include: ['histos_dnn.yml']\n")
        else:
            if args.yield_only:
                f1.write("\nplots:\n  include: ['histos_yield.yml']\n")
            else:
                if os.path.exists(os.path.join(config_path, 'histos_v1.yml')):
                    f1.write("\nplots:\n  include: ['histos_v1.yml', 'histos_reco.yml']\n")
                else:
                    f1.write("\nplots:\n  include: ['histos_yield.yml', 'histos_control.yml', 'histos_reco.yml']\n")

# Execute plotIt
plotit_bin = '../plotIt/plotIt'
out_fig_arg = os.path.join(dest_path, 'figure_run3' + postfix)
config_arg = os.path.join(config_path, 'config_Run3.yml')

if args.yield_only:
    call([plotit_bin, '-o ' + out_fig_arg, config_arg, '-y', '-a'], shell=False)
else:
    call([plotit_bin, '-o ' + out_fig_arg, config_arg], shell=False)

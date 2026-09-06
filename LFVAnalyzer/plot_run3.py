# use (process plots): python plot_run3.py -I inputdir/channel
# use (DNN plots): python plot_run3.py -I inputdir/channel -D

import sys, os
import re
from subprocess import call
import argparse

parser = argparse.ArgumentParser(description="Run plotIt for Run 3")
parser.add_argument('-I', '--input', dest='input', type=str, default="process_0715_v4/electron", help="Input directory path")
parser.add_argument('-C', '-ch', '--channel', dest='channel', type=str, default="", choices=["", "electron", "muon"], help="Channel: electron or muon (auto-detected from input if empty)")
parser.add_argument("-D", dest="DNN", action="store_true", default=False, help="Run for DNN histograms")
parser.add_argument("-F", '--ff', dest="tauFF", action="store_true", default=False, help="Run for tauFF histograms")
parser.add_argument("-y", dest="yield_only", action="store_true", default=False, help="Run for yield histograms")
parser.add_argument("-q", dest="qcd", action="store_true", default=False, help="include QCD")
parser.add_argument("--postfix", dest="postfix", type=str, default="", help="Add postfix to output here, to have rebinning for histograms")
args = parser.parse_args()

dest_path = args.input.rstrip('/')
postfix = args.postfix

# Auto-detect channel if not specified
channel = args.channel
if not channel:
    if 'muon' in dest_path.lower():
        channel = 'muon'
    else:
        channel = 'electron'

config_path = f'../plotIt/configs/Run3_{channel}/'

syst_avoid = ['jesFlavorPureBottom', 'jesFlavorPureCharm', 'jer', 'jesRelativeBal', 'muonhighscale', 'metUnclust', 'jesAbsolute', 'tes', 'jesFlavorPureQuark',]

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
    matched_rel_path = None
    for cand in era_info['dir_candidates']:
        subpaths_to_test = [
            os.path.join(channel, f"{cand}_postprocess{postfix}"),
            os.path.join(channel, f"{cand}{postfix}"),
            os.path.join(channel, cand),
            f"{cand}_postprocess{postfix}",
            f"{cand}{postfix}",
            cand
        ]
        for sub in subpaths_to_test:
            full_p = os.path.join(dest_path, sub)
            if os.path.exists(full_p):
                matched_rel_path = sub
                break
        if matched_rel_path:
            break

    if not matched_rel_path:
        matched_rel_path = f"{era_info['dir_candidates'][0]}_postprocess{postfix}"
    resolved_eras[era_key] = matched_rel_path

has_merged_run3 = os.path.exists(os.path.join(dest_path, 'Run3' + postfix)) or os.path.exists(os.path.join(dest_path, channel, 'Run3' + postfix))

string_for_files = ''
for era_key, era_info in eras.items():
    lumi = era_info['lumi']
    era_rel_path = resolved_eras[era_key]
    files_yml_path = ''
    if '2024' in era_key:
        files_yml_path = os.path.join(config_path, 'files24.yml')
    else:
        files_yml_path = os.path.join(config_path, 'files.yml')

    if os.path.exists(files_yml_path):
        with open(files_yml_path) as f:
            lines = f.readlines()
            skip_signal = False
            for line in lines:
                raw_line = line
                clean_line = line.strip()

                if skip_signal and 'hist' in clean_line: skip_signal = False
                if 'LFV' in clean_line and has_merged_run3: skip_signal = True
                if clean_line.startswith('#'): skip_signal = True

                if clean_line.startswith("'") and "':" in clean_line and not clean_line.startswith('#'):
                    fname_with_quotes = clean_line[:clean_line.find("':")]
                    fname = fname_with_quotes[1:]
                    basename = os.path.basename(fname)

                    file_entry_path = os.path.normpath(os.path.join(dest_path, era_rel_path, basename))
                    line = f"'{file_entry_path}':\n"

                    if not any(i in line for i in ['LFV', 'SingleMuon', 'Egamma', 'EGamma', 'Muon', 'WJetsToLNu_HT0To100']):
                        line += '  scale: ' + str(float(lumi) / total_lumi) + '\n'
                    elif 'WJetsToLNu_HT0To100' in line:
                        line += '  scale: ' + str(1.0288 * float(lumi) / total_lumi) + '\n'
                else:
                    line = raw_line

                #if 'group' in line:
                #    if 'yields-group' in line:
                #        if any(i in line for i in ['ttbar LL']):
                #            line = "  yields-group: '1\\TTLL' \n"
                #        elif any(i in line for i in ['ttbar LJ']):
                #            line = "  yields-group: '2\\TTLJ' \n"
                #        elif any(i in line for i in ['ttbar Had']):
                #            line = "  yields-group: '3TTHad' \n"
                #        elif any(i in line for i in ['Z+Jets']):
                #            line = "  yields-group: '4ZJets' \n"
                #        elif any(i in line for i in ['W+Jets']):
                #            line = "  yields-group: '6WJets' \n"
                #        elif any(i in line for i in ['tt+X']):
                #            line = "  yields-group: '7ttV' \n"
                #        elif any(i in line for i in ['VV']):
                #            line = "  yields-group: '8VV' \n"
                #        elif any(i in line for i in ['QCD']):
                #            line = "  yields-group: '9QCD' \n"
                #    else:
                #        if any(i in line for i in ['Gttlj', 'Gttll']):
                #            line = '  group: Gtt \n'
                #        elif any(i in line for i in ['Gttjj', 'GttV', 'GVV', 'GZJets', 'GWJets']):
                #            line = '  group: Gother \n'

                if not skip_signal:
                    string_for_files += line

# Write files_Run3.yml
files_run3_yml = os.path.join(config_path, 'files_Run3.yml')
with open(files_run3_yml, 'w+') as fnew:
    if has_merged_run3:
        run3_dir_path = os.path.normpath(os.path.join(dest_path, channel, 'Run3' + postfix)) if os.path.exists(os.path.join(dest_path, channel, 'Run3' + postfix)) else os.path.normpath(os.path.join(dest_path, 'Run3' + postfix))
        if channel == 'electron' and args.yield_only:
            fnew.write(f"""
'{run3_dir_path}/hist_TUETau-LFV-Scalar.root':
  type: signal
  pretty-name: 'LFVSTus'
  cross-section: 0.1139
  generated-events: 1
  group: GLFVSTus
  yields-group: '1ST LFV \\cPqt\\cPqe\\tau Scalar'
  order: 1

'{run3_dir_path}/hist_TUETau-LFV-Vector.root':
  type: signal
  pretty-name: 'LFVSTuv'
  cross-section: 0.5267
  generated-events: 1
  group: GLFVSTuv
  yields-group: '1ST LFV \\cPqt\\cPqe\\tau Vector'
  order: 1

'{run3_dir_path}/hist_TUETau-LFV-Tensor.root':
  type: signal
  pretty-name: 'LFVSTut'
  cross-section: 2.411
  generated-events: 1
  group: GLFVSTut
  yields-group: '2ST LFV \\cPqt\\cPqe\\tau Tensor'
  order: 1

'{run3_dir_path}/hist_TCETau-LFV-Scalar.root':
  type: signal
  pretty-name: 'LFVSTcs'
  cross-section: 0.007285
  generated-events: 1
  group: GLFVSTcs
  yields-group: '3ST LFV \\cPqt\\cPqc\\tau Scalar'
  order: 1

'{run3_dir_path}/hist_TCETau-LFV-Vector.root':
  type: signal
  pretty-name: 'LFVSTcv'
  cross-section: 0.03788
  generated-events: 1
  group: GLFVSTcv
  yields-group: '3ST LFV \\cPqt\\cPqc\\tau Vector'
  order: 3

'{run3_dir_path}/hist_TCETau-LFV-Tensor.root':
  type: signal
  pretty-name: 'LFVSTct'
  cross-section: 0.1874
  generated-events: 1
  group: GLFVSTct
  yields-group: '4ST LFV \\cPqt\\cPqc\\tau Tensor'
  order: 1

'{run3_dir_path}/hist_TTtoUETau-LFV-Scalar.root':
  type: signal
  pretty-name: 'LFVTTus'
  cross-section: 0.002003
  generated-events: 1
  group: GLFVTTus
  yields-group: '5TT LFV \\cPqt\\cPqu\\e\\tau Scalar'
  order: 1

'{run3_dir_path}/hist_TTtoUETau-LFV-Vector.root':
  type: signal
  pretty-name: 'LFVTTuv'
  cross-section: 0.01604
  generated-events: 1
  group: GLFVTTuv
  yields-group: '5TT LFV \\cPqt\\cPqu\\e\\tau Vector'
  order: 2

'{run3_dir_path}/hist_TTtoUETau-LFV-Tensor.root':
  type: signal
  pretty-name: 'LFVTTut'
  cross-section: 0.09636
  generated-events: 1
  group: GLFVTTut
  yields-group: '6TT LFV \\cPqt\\cPqu\\e\\tau Tensor'
  order: 1

'{run3_dir_path}/hist_TTtoCETau-LFV-Scalar.root':
  type: signal
  pretty-name: 'LFVTTcs'
  cross-section: 0.002006
  generated-events: 1
  group: GLFVTTcs
  yields-group: '7TT LFV \\cPqt\\cPqc\\e\\tau Scalar'
  order: 1

'{run3_dir_path}/hist_TTtoCETau-LFV-Vector.root':
  type: signal
  pretty-name: 'LFVTTcv'
  cross-section: 0.01583
  generated-events: 1
  group: GLFVTTcv
  yields-group: '7TT LFV \\cPqt\\cPqc\\e\\tau Vector'
  order: 4

'{run3_dir_path}/hist_TTtoCETau-LFV-Tensor.root':
  type: signal
  pretty-name: 'LFVTTct'
  cross-section: 0.09645
  generated-events: 1
  group: GLFVTTct
  yields-group: '8TT LFV \\cPqt\\cPqc\\e\\tau Tensor'
  order: 1
""")
        elif channel == 'muon' and args.yield_only: # muon
            fnew.write(f"""
'{run3_dir_path}/hist_TUMuTau-LFV-Scalar.root':
  type: signal
  pretty-name: 'LFVSTus'
  cross-section: 0.1138
  generated-events: 1
  group: GLFVSTus
  yields-group: '1ST LFV \\cPqt\\cPqu\\mu\\tau Scalar'
  order: 1

'{run3_dir_path}/hist_TUMuTau-LFV-Vector.root':
  type: signal
  pretty-name: 'LFVSTuv'
  cross-section: 0.5266
  generated-events: 1
  group: GLFVSTuv
  yields-group: '1ST LFV \\cPqt\\cPqu\\mu\\tau Vector'
  order: 1

'{run3_dir_path}/hist_TUMuTau-LFV-Tensor.root':
  type: signal
  pretty-name: 'LFVSTut'
  cross-section: 2.401
  generated-events: 1
  group: GLFVSTut
  yields-group: '2ST LFV \\cPqt\\cPqu\\mu\\tau Tensor'
  order: 1

'{run3_dir_path}/hist_TCMuTau-LFV-Scalar.root':
  type: signal
  pretty-name: 'LFVSTcs'
  cross-section: 0.007288
  generated-events: 1
  group: GLFVSTcs
  yields-group: '3ST LFV \\cPqt\\cPqc\\mu\\tau Scalar'
  order: 1

'{run3_dir_path}/hist_TCMuTau-LFV-Vector.root':
  type: signal
  pretty-name: 'LFVSTcv'
  cross-section: 0.03753
  generated-events: 1
  group: GLFVSTcv
  yields-group: '3ST LFV \\cPqt\\cPqc\\mu\\tau Vector'
  order: 3

'{run3_dir_path}/hist_TCMuTau-LFV-Tensor.root':
  type: signal
  pretty-name: 'LFVSTct'
  cross-section: 0.1849
  generated-events: 1
  group: GLFVSTct
  yields-group: '4ST LFV \\cPqt\\cPqc\\mu\\tau Tensor'
  order: 1

'{run3_dir_path}/hist_TTtoUMuTau-LFV-Scalar.root':
  type: signal
  pretty-name: 'LFVTTus'
  cross-section: 0.0044485
  generated-events: 1
  group: GLFVTTus
  yields-group: '5TT LFV \\cPqt\\cPqu\\mu\\tau Scalar'
  order: 1

'{run3_dir_path}/hist_TTtoUMuTau-LFV-Vector.root':
  type: signal
  pretty-name: 'LFVTTuv'
  cross-section: 0.035555
  generated-events: 1
  group: GLFVTTuv
  yields-group: '5TT LFV \\cPqt\\cPqu\\mu\\tau Vector'
  order: 2

'{run3_dir_path}/hist_TTtoUMuTau-LFV-Tensor.root':
  type: signal
  pretty-name: 'LFVTTut'
  cross-section: 0.21333
  generated-events: 1
  group: GLFVTTut
  yields-group: '6TT LFV \\cPqt\\cPqu\\mu\\tau Tensor'
  order: 1

'{run3_dir_path}/hist_TTtoCMuTau-LFV-Scalar.root':
  type: signal
  pretty-name: 'LFVTTcs'
  cross-section: 0.0044485
  generated-events: 1
  group: GLFVTTcs
  yields-group: '7TT LFV \\cPqt\\cPqc\\mu\\tau Scalar'
  order: 1

'{run3_dir_path}/hist_TTtoCMuTau-LFV-Vector.root':
  type: signal
  pretty-name: 'LFVTTcv'
  cross-section: 0.035555
  generated-events: 1
  group: GLFVTTcv
  yields-group: '7TT LFV \\cPqt\\cPqc\\mu\\tau Vector'
  order: 4

'{run3_dir_path}/hist_TTtoCMuTau-LFV-Tensor.root':
  type: signal
  pretty-name: 'LFVTTct'
  cross-section: 0.21333
  generated-events: 1
  group: GLFVTTct
  yields-group: '8TT LFV \\cPqt\\cPqc\\mu\\tau Tensor'
  order: 1
""")
        elif channel == 'electron' and not args.yield_only:
            fnew.write(f"""
'{run3_dir_path}/hist_TUETau-LFV-Vector.root':
  type: signal
  pretty-name: 'LFVSTuv'
  cross-section: 0.5267
  generated-events: 1
  group: GLFVSTuv
  yields-group: '3ST LFV \\cPqt\\cPqu\\e\\tau Vector'
  order: 3

'{run3_dir_path}/hist_TTtoUETau-LFV-Vector.root':
  type: signal
  pretty-name: 'LFVTTuv'
  cross-section: 0.01604
  generated-events: 1
  group: GLFVTTuv
  yields-group: '7TT LFV \\cPqt\\cPqu\\e\\tau Vector'
  order: 4
""")
        elif channel == 'muon' and not args.yield_only: # muon
            fnew.write(f"""
'{run3_dir_path}/hist_TUMuTau-LFV-Vector.root':
  type: signal
  pretty-name: 'LFVSTuv'
  cross-section: 0.5266
  generated-events: 1
  group: GLFVSTuv
  yields-group: '3ST LFV \\cPqt\\cPqu\\mu\\tau Vector'
  order: 3

'{run3_dir_path}/hist_TTtoUMuTau-LFV-Vector.root':
  type: signal
  pretty-name: 'LFVTTuv'
  cross-section: 0.035555
  generated-events: 1
  group: GLFVTTuv
  yields-group: '7TT LFV \\cPqt\\cPqu\\mu\\tau Vector'
  order: 4
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

        if args.tauFF:
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
                if os.path.exists(os.path.join(config_path, 'histos.yml')):
                    f1.write("\nplots:\n  include: ['histos_yield.yml', 'histos.yml', 'histos_reco.yml']\n")
                else:
                    f1.write("\nplots:\n  include: ['histos_yield.yml', 'histos_control.yml', 'histos_reco.yml']\n")

# Execute plotIt
plotit_bin = '../plotIt/plotIt'
out_fig_arg = os.path.join(dest_path, 'figure_run3' + postfix)
if args.qcd:
    out_fig_arg = os.path.join(dest_path, 'figure_run3_qcd' + postfix)
    os.makedirs(out_fig_arg, exist_ok=True)
config_arg = os.path.join(config_path, 'config_Run3.yml')

if args.yield_only and not args.qcd:
    call([plotit_bin, '-o ' + out_fig_arg, config_arg, '-y'], shell=False)
elif args.yield_only and args.qcd:
    call([plotit_bin, '-o ' + out_fig_arg, config_arg, '-y', '-q'], shell=False)
elif args.qcd:
    call([plotit_bin, '-o ' + out_fig_arg, config_arg, '-q', '-a'], shell=False)
else:
    call([plotit_bin, '-o ' + out_fig_arg, config_arg, '-a'], shell=False)

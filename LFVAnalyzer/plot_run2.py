import sys, os
import re
from subprocess import call
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-I', '--input', dest='input', type=str, default="test")
parser.add_argument("-D", dest="DNN", action="store_true", default=False, help="Run for DNN histograms")
parser.add_argument("-y", dest="yield_only", action="store_true", default=False, help="Run for DNN histograms")
parser.add_argument("--postfix", dest="postfix", type=str, default="", help="Add postfix to output here, to have rebinning for histograms")
args = parser.parse_args()
input = args.input
yield_only = args.yield_only

syst_avoid = ['muonhighscale', 'metUnclust']

do_postfit_scale = False
postfit_scales = {'2018': {'tt': 0.9095, 'other': 0.964, 'singleTop': 0.9528, 'st_lfv_ut': 0.9544}, '2017': {'tt': 0.9513, 'other': 1.0127, 'singleTop': 0.9977, 'st_lfv_ut': 0.9667}, '2016pre': {'tt': 0.9568, 'other': 1.028, 'singleTop': 1.0091, 'st_lfv_ut': 0.9714}, '2016post': {'tt': 0.9569, 'other': 1.0272, 'singleTop': 1.0015, 'st_lfv_ut': 0.968}}

#########################################################
# WARNING: Do NOT use run2 systematic table!!
# 1) no proper correlations
# 2) this add up only nominal signals, not considering unc
# Therefore, yearly unc from weight Unc. are exploding
##########################################################

config_path = '../plotIt/configs/TOP-22-011/'
common_syst = 'systematics:\n'
groups = ['GData', 'Gttll', 'Gttlj', 'Gttjj', 'GttV', 'GZJets', 'GWJets', 'GSingleT', 'GVV', 'GQCD',
          'GLFVSTcv', 'GLFVSTuv', 'GLFVTTcv', 'GLFVTTuv',]
years = {'2016pre': 19502, '2016post': 16812, '2017': 41480, '2018':59832}

common_syst_list = []
for year, lumi in years.items():
    with open(config_path + 'config_' + year + '.yml') as f:
        lines = f.readlines()
        for line in lines:
            if any(s in line for s in syst_avoid): continue
            if "  - " in line: common_syst_list.append(line)

common_syst_list = sorted(set(common_syst_list), key=common_syst_list.index)

dest_path = input
if not os.path.exists(os.path.join(dest_path, 'figure_run2' + args.postfix)):
    try: os.makedirs(os.path.join(dest_path, 'figure_run2' + args.postfix))
    except: pass
if not os.path.exists(os.path.join(dest_path, 'figure_run2_postfit_scale' + args.postfix)) and do_postfit_scale:
    try: os.makedirs(os.path.join(dest_path, 'figure_run2_postfit_scale' + args.postfix))
    except: pass
#if not os.path.exists(os.path.join(dest_path, 'figure_run2/qcd')):
#    try: os.makedirs(os.path.join(dest_path, 'figure_run2/qcd'))
#    except: pass
#if not os.path.exists(os.path.join(dest_path, 'figure_run2/dyincl')):
#    try: os.makedirs(os.path.join(dest_path, 'figure_run2/dyincl'))
#    except: pass
#if not os.path.exists(os.path.join(dest_path, 'figure_run2/dyincl/qcd')):
#    try: os.makedirs(os.path.join(dest_path, 'figure_run2/dyincl/qcd'))
#    except: pass

for item in common_syst_list:
    common_syst += item

string_for_files = ''
for year, lumi in years.items():
    #Firstly, merge file list + scale
    with open(config_path + 'files.yml') as f:
    #with open(config_path + 'files_fake.yml') as f:
        lines = f.readlines()
        skip_signal = False
        for line in lines:
            #if '#' in line[0]: line = line[1:]
            if skip_signal and 'hist' in line: skip_signal = False
            if 'LFV' in line: skip_signal = True
            #if 'hist_QCD' in line: skip_signal = True
            if '#' in line[0]: skip_signal = True
            if 'hist' in line:
                line = line[0] + dest_path + '/' + year + '_postprocess' + args.postfix + '/' + line[1:]
                if not any(i in line for i in ['LFV', 'SingleMuon', 'WJetsToLNu_HT0To100']):
                    if do_postfit_scale: #FIXME: Hard-coded
                        if any(i in line for i in ['TTToSemiLeptonic', 'TTTo2L2Nu']):
                            line += "  scale: " + str(postfit_scales[year]['tt'] * int(lumi)/137625.0) + ' \n'
                        elif any(i in line for i in ['TTToHadronic', 'ttH', 'TTW', 'TTZ', 'WW', 'WZ', 'ZZ', 'DYJets', 'WJets']):
                            line += "  scale: " + str(postfit_scales[year]['other'] * int(lumi)/137625.0) + ' \n'
                        elif any(i in line for i in ['ST_s', 'ST_t']):
                            line += "  scale: " + str(postfit_scales[year]['singleTop'] * int(lumi)/137625.0) + ' \n'
                    else:
                        line += '  scale: ' + str(int(lumi)/137625.0) + '\n'
                elif 'WJetsToLNu_HT0To100' in line:
                    if do_postfit_scale:
                        line += '  scale: ' + str(postfit_scales[year]['other'] * 1.0288 * int(lumi)/137625.0) + '\n'
                    else:
                        line += '  scale: ' + str(1.0288 * int(lumi)/137625.0) + '\n'
            # Merge processes. This is only needed for paper plots - i.e. Run2 plots
            if 'group' in line:
                #First, yield group
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
            #if not skip_signal and not any(i in line for i in ['yields-group']):
            if not skip_signal:
                #if 'group' in line and not any(i in line for i in groups): string_for_files += '  group: Gother \n'
                #else: string_for_files += line
                string_for_files += line


with open(config_path + 'files_Run2.yml', 'w+') as fnew:
    print("""
'{0}/Run2{1}/hist_ST_LFV_TUMuTau_Scalar.root':
  type: signal
  pretty-name: 'LFVSTus'
  cross-section: 0.08494
  generated-events: 1
  group: GLFVSTus
  #yields-group: '1ST LFV$^{{\\cPqt\\cPqu\\mu\\tau}}_{{Scalar}}$'
  yields-group: '1ST LFV $\\cPqt\\cPqu\\mu\\tau$ Scalar'
  order: 1

'{0}/Run2{1}/hist_ST_LFV_TUMuTau_Vector.root':
  type: signal
  pretty-name: 'LFVSTuv'
  cross-section: 0.3924
  generated-events: 1
  group: GLFVSTuv
  #yields-group: '1ST LFV$^{{\\cPqt\\cPqu\\mu\\tau}}_{{Vector}}$'
  yields-group: '1ST LFV $\\cPqt\\cPqu\\mu\\tau$ Vector'
  order: 1

'{0}/Run2{1}/hist_ST_LFV_TUMuTau_Tensor.root':
  type: signal
  pretty-name: 'LFVSTut'
  cross-section: 1.781
  generated-events: 1
  group: GLFVSTut
  #yields-group: '2ST LFV$^{{\\cPqt\\cPqu\\mu\\tau}}_{{Tensor}}$'
  yields-group: '2ST LFV $\\cPqt\\cPqu\\mu\\tau$ Tensor'
  order: 1

'{0}/Run2{1}/hist_ST_LFV_TCMuTau_Scalar.root':
  type: signal
  pretty-name: 'LFVSTcs'
  cross-section: 0.007375
  generated-events: 1
  group: GLFVSTcs
  #yields-group: '3ST LFV$^{{\\cPqt\\cPqc\\mu\\tau}}_{{Scalar}}$'
  yields-group: '3ST LFV $\\cPqt\\cPqc\\mu\\tau$ Scalar'
  order: 1

'{0}/Run2{1}/hist_ST_LFV_TCMuTau_Vector.root':
  type: signal
  pretty-name: 'LFVSTcv'
  cross-section: 0.03673
  generated-events: 1
  group: GLFVSTcv
  #yields-group: '3ST LFV$^{{\\cPqt\\cPqc\\mu\\tau}}_{{Vector}}$'
  yields-group: '3ST LFV $\\cPqt\\cPqc\\mu\\tau$ Vector'
  order: 3

'{0}/Run2{1}/hist_ST_LFV_TCMuTau_Tensor.root':
  type: signal
  pretty-name: 'LFVSTct'
  cross-section: 0.177
  generated-events: 1
  group: GLFVSTct
  #yields-group: '4ST LFV$^{{\\cPqt\\cPqc\\mu\\tau}}_{{Tensor}}$'
  yields-group: '4ST LFV $\\cPqt\\cPqc\\mu\\tau$ Tensor'
  order: 1

'{0}/Run2{1}/hist_TT_LFV_TUMuTau_Scalar.root':
  type: signal
  pretty-name: 'LFVTTus'
  cross-section: 0.00269
  generated-events: 1
  group: GLFVTTus
  #yields-group: '5TT LFV$^{{\\cPqt\\cPqu\\mu\\tau}}_{{Scalar}}$'
  yields-group: '5TT LFV $\\cPqt\\cPqu\\mu\\tau$ Scalar'
  order: 1

'{0}/Run2{1}/hist_TT_LFV_TUMuTau_Vector.root':
  type: signal
  pretty-name: 'LFVTTuv'
  cross-section: 0.0215
  generated-events: 1
  group: GLFVTTuv
  #yields-group: '5TT LFV$^{{\\cPqt\\cPqu\\mu\\tau}}_{{Vector}}$'
  yields-group: '5TT LFV $\\cPqt\\cPqu\\mu\\tau$ Vector'
  order: 2

'{0}/Run2{1}/hist_TT_LFV_TUMuTau_Tensor.root':
  type: signal
  pretty-name: 'LFVTTut'
  cross-section: 0.1290
  generated-events: 1
  group: GLFVTTut
  #yields-group: '6TT LFV$^{{\\cPqt\\cPqu\\mu\\tau}}_{{Tensor}}$'
  yields-group: '6TT LFV $\\cPqt\\cPqu\\mu\\tau$ Tensor'
  order: 1

'{0}/Run2{1}/hist_TT_LFV_TCMuTau_Scalar.root':
  type: signal
  pretty-name: 'LFVTTcs'
  cross-section: 0.00269
  generated-events: 1
  group: GLFVTTcs
  #yields-group: '7TT LFV$^{{\\cPqt\\cPqc\\mu\\tau}}_{{Scalar}}$'
  yields-group: '7TT LFV $\\cPqt\\cPqc\\mu\\tau$ Scalar'
  order: 1

'{0}/Run2{1}/hist_TT_LFV_TCMuTau_Vector.root':
  type: signal
  pretty-name: 'LFVTTcv'
  cross-section: 0.0215
  generated-events: 1
  group: GLFVTTcv
  #yields-group: '7TT LFV$^{{\\cPqt\\cPqc\\mu\\tau}}_{{Vector}}$'
  yields-group: '7TT LFV $\\cPqt\\cPqc\\mu\\tau$ Vector'
  order: 4

'{0}/Run2{1}/hist_TT_LFV_TCMuTau_Tensor.root':
  type: signal
  pretty-name: 'LFVTTct'
  cross-section: 0.1290
  generated-events: 1
  group: GLFVTTct
  #yields-group: '8TT LFV$^{{\\cPqt\\cPqc\\mu\\tau}}_{{Tensor}}$'
  yields-group: '8TT LFV $\\cPqt\\cPqc\\mu\\tau$ Tensor'
  order: 1
  """.format(dest_path, args.postfix), file=fnew)
    fnew.write(string_for_files)

with open(config_path + 'template_Run2.yml') as f:
    lines = f.readlines()
    with open(config_path + 'config_Run2.yml', 'w+') as f1:
        for line in lines: f1.write(line)
        f1.write(common_syst)
        if 'FF' in dest_path:
            if yield_only:
                f1.write("\nplots:\n  include: ['histos_yield_S5.yml']\n")
            else:
                f1.write("\nplots:\n  include: ['histos_FFapply.yml', 'histos_yield_S5.yml']\n")
        elif args.DNN:
            f1.write("\nplots:\n  include: ['histos_dnn.yml']\n")
        else:
            if yield_only:
                f1.write("\nplots:\n  include: ['histos_yield.yml']\n")
            else:
                f1.write("\nplots:\n  include: ['histos_yield.yml', 'histos_control.yml', 'histos_reco.yml']\n")

if yield_only:
    call(['../plotIt/plotIt', '-o ' + dest_path + '/figure_run2' + args.postfix, config_path + 'config_Run2.yml', '-y', '-a'], shell=False)
    #call(['../plotIt/plotIt', '-o ' + dest_path + '/figure_run2/qcd', config_path + 'config_Run2.yml', '-y', '-q', '-a'], shell=False)
    #call(['../plotIt/plotIt', '-o ' + dest_path + '/figure_run2/dyincl', config_path + 'config_Run2.yml', '-y', '-d', '-a'], shell=False)
    #call(['../plotIt/plotIt', '-o ' + dest_path + '/figure_run2/dyincl/qcd', config_path + 'config_Run2.yml', '-y', '-d', '-q', '-a'], shell=False)
elif do_postfit_scale:
    call(['../plotIt/plotIt', '-o ' + dest_path + '/figure_run2_postfit_scale' + args.postfix, config_path + 'config_Run2.yml', '-y'], shell=False)
else:
    call(['../plotIt/plotIt', '-o ' + dest_path + '/figure_run2' + args.postfix, config_path + 'config_Run2.yml'], shell=False)
    #call(['../plotIt/plotIt', '-o ' + dest_path + '/figure_run2/qcd', config_path + 'config_Run2.yml', '-q'], shell=False)
    #call(['../plotIt/plotIt', '-o ' + dest_path + '/figure_run2/dyincl', config_path + 'config_Run2.yml', '-d'], shell=False)
    #call(['../plotIt/plotIt', '-o ' + dest_path + '/figure_run2/dyincl/qcd', config_path + 'config_Run2.yml', '-d', '-q'], shell=False)

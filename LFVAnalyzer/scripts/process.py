import os, sys, re, argparse
from subprocess import call

parser = argparse.ArgumentParser(usage="%prog [options]")
parser.add_argument("-V", "--version", dest="version", type=str, default="", help="Skim version: folder under /data2/common/skimmed_NanoAOD/")
parser.add_argument("-O", "--outdir", dest="outdir", type=str, default="test", help="Output folder in your working directory")
parser.add_argument("-Y", "--year", dest="year", type=str, default="", help="Select 2016pre, 2016post, 2017, or 2018 runs")
parser.add_argument("-C", "--ch",  dest="ch", type=str, default="", help="Select electron or muon")
parser.add_argument("-S", "--syst", dest="syst", type=str, default="theory", help="Systematic: 'data' for Data, 'nosyst' for mc without uncertainties. Default is 'theory'. To run without theory unc for TT samples, put 'all'.")
parser.add_argument("-D", "--dataset", dest="dataset", action="store", nargs="+", default=[], help="Put dataset folder name (eg. TTTo2L2Nu) to process specific one.")
parser.add_argument("-F", "--dataOrMC", dest="dataOrMC", type=str, default="", help="data or mc flag, if you want to process data-only or mc-only")
parser.add_argument("--dry", dest="dry", action="store_true", default=False, help="dryrun: not submitting jobs to slurm")
parser.add_argument("-M", "--mode", dest="mode", type=str, default="", help="Only for fake rate: lss, los, tss, tos")
parser.add_argument("--ff", dest="ff", action="store_true", default=False, help="Apply tau fake factor for final selection")
options = parser.parse_args()

year = options.year
ch   = options.ch
workdir = os.getcwd()
indir = os.path.join('/data2/common/skimmed_NanoAOD/', options.version, ch)
tgdir = os.path.join(workdir, options.outdir, ch, year)
logdir = os.path.join(workdir, options.outdir, ch, year, 'log')
splitList = ["TTto4Q",
#splitList = ["TTto2L2Nu", "TTtoLNu2Q", "TTto4Q",
             #"ST_LFV_TCMuTau_Scalar", "ST_LFV_TCMuTau_Vector", "ST_LFV_TCMuTau_Tensor",
             #"ST_LFV_TUMuTau_Scalar", "ST_LFV_TUMuTau_Vector", "ST_LFV_TUMuTau_Tensor",
             #"TT_LFV_TCMuTau_Scalar", "TT_LFV_TCMuTau_Vector", "TT_LFV_TCMuTau_Tensor",
             #"TT_LFV_TUMuTau_Scalar", "TT_LFV_TUMuTau_Vector", "TT_LFV_TUMuTau_Tensor",
            ]

os.makedirs(tgdir, exist_ok=True)
os.makedirs(logdir, exist_ok=True)

datadir = os.path.join(indir, 'data', year)
data_list = [os.path.join(datadir, s) for s in os.listdir(datadir)]
mcdir = os.path.join(indir, 'mc', year)
mc_list = [os.path.join(mcdir, s) for s in os.listdir(mcdir)]

dataset_list = data_list + mc_list

syst_list = ["", "__tesup", "__tesdown", "__jerup","__jerdown", "__jesAbsoluteup","__jesAbsolutedown",
             "__jesAbsolute_"+year[:4]+"up", "__jesAbsolute_"+year[:4]+"down",
             "__jesBBEC1up", "__jesBBEC1down", "__jesBBEC1_"+year[:4]+"up", "__jesBBEC1_"+year[:4]+"down",
             #"__jesFlavorQCDup", "__jesFlavorQCDdown", "__jesRelativeBalup", "__jesRelativeBaldown",
             "__jesRelativeBalup", "__jesRelativeBaldown",
             "__jesRelativeSample_"+year[:4]+"up", "__jesRelativeSample_"+year[:4]+"down"]
if year == "2018": syst_list.extend(["__jesHEMup", "__jesHEMdown"])
syst_list.extend(["__jesFlavorPureGluonup", "__jesFlavorPureGluondown", "__jesFlavorPureQuarkup", "__jesFlavorPureQuarkdown",
                  "__jesFlavorPureCharmup", "__jesFlavorPureCharmdown", "__jesFlavorPureBottomup", "__jesFlavorPureBottomdown"])
syst_list.extend(["__metUnclustup", "__metUnclustdown"])
syst_list.extend(["__muonhighscaleup", "__muonhighscaledown"])

# tune and hdamp will appear as an individual dataset.
# thus no need to run in loop, but left here for double check
syst_ext = ["__tuneup", "__tunedown", "__hdampup", "__hdampdown",]

# Default syst for MC = "theory" = all SF + theory
# To run central values only, pass "nosyst"
# To calculate theory uncertainties (PDF, ME, PS, if available), put 'theory'
# Data should use systematic flag "data" not to construct event weights
if options.syst == "nosyst": syst_list = [""]

#For split processing

parameters = [] #order: (tgdir, indir, outdir, outfilename, syst)
for ds in dataset_list:

    #if any(n in ds for n in ['QCD_']):
    #    continue

    if len(options.mode) > 0 and any(n in ds for n in ['_LFV_']):
        continue

    # will use data/signal in genuine folder
    if 'fakeTau' in options.outdir and any(n in ds for n in ['_LFV_', 'SingleMuon', 'QCD_']):
        continue
    elif 'genuineTau' in options.outdir and any(n in ds for n in ['QCD_']):
        continue

    #Split TT2L2Nu and TTToSemilep into files
    toSplit = False
    if any(name_ in ds for name_ in splitList) and '__' not in ds and "2024" in year and "muon" in ch: toSplit = True
    rootfilestoprocess = []
    fullnamelist =[]
    if toSplit:
        #if '__' in ds: continue
        print("collecting root files in "+ds)
        flist = os.listdir(ds)
        for fname in flist:
            fullname = os.path.join(ds, fname)
            fullnamelist.append(fullname)
        for fname in fullnamelist:
            if re.match('.*\.root', fname) and os.path.isfile(fname): # if it has .root file extension
                rootfilestoprocess.append(fname)

    if len(options.dataset) > 0 and not any(i in ds for i in options.dataset): continue

    for src in syst_list:

        if options.dataOrMC == 'data':
            if 'data' not in ds[5:]: continue
        elif options.dataOrMC == 'mc':
            if 'mc' not in ds[5:]: continue

        dataset_name = ds.split('/')[-1]
        outdir = tgdir
        outfname = 'hist_' + dataset_name + src + '.root'

        ext_syst = False
        if any(i in dataset_name for i in syst_ext): ext_syst = True

        if 'data' in ds[5:]:
            if src == "":
                parameters.append([year, ch, ds, outdir, outfname, "data"])
            else: continue

        else:
            if src == "" and not ext_syst:
                if options.syst == "all":
                    if toSplit:
                        os.makedirs(tgdir.replace(year, year + '/' + "split"), exist_ok=True)
                        parameters.append((year, ch, rootfilestoprocess, outdir.replace(year, year + '/' + "split"), outfname, "all"))
                    else:
                        parameters.append((year, ch, ds, outdir, outfname, "all"))
                elif options.syst == "theory":
                    if any(i in dataset_name for i in ["TTto", "ST_t", "TT_LFV", "ST_LFV"]):
                        if toSplit:
                            os.makedirs(tgdir.replace(year, year + '/' + "split"), exist_ok=True)
                            parameters.append([year, ch, rootfilestoprocess, outdir.replace(year, year + '/' + "split"), outfname, "theory"])
                        else:
                            parameters.append([year, ch, ds, outdir, outfname, "theory"])
                    else: parameters.append([year, ch, ds, outdir, outfname, "all"])
                elif options.syst == "nosyst":
                    if toSplit:
                        os.makedirs(tgdir.replace(year, year + '/' + "split"), exist_ok=True)
                        parameters.append([year, ch, rootfilestoprocess, outdir.replace(year, year + '/' + "split"), outfname, "nosyst"])
                    else:
                        parameters.append([year, ch, ds, outdir, outfname, "nosyst"])
            elif src == "" and ext_syst:
                if any(i in dataset_name for i in syst_ext) and options.syst == "nosyst": continue
                else: parameters.append([year, ch, ds, outdir, outfname, "nosyst"])
            else:
                #os.makedir(os.path.join(tgdir, dataset_name+src)
                if any(i in dataset_name for i in syst_ext): continue
                parameters.append([year, ch, ds, outdir, outfname, src[2:]])
                # for tests
                #if toSplit:
                #    os.makedirs(tgdir.replace(year, year + '/' + "split"), exist_ok=True)
                #    parameters.append([year, rootfilestoprocess, outdir.replace(year, year + '/' + "split"), outfname, src[2:]])
                #else:
                #    parameters.append([year, ds, outdir, outfname, src[2:]])


runString_list = []
for item in parameters:
    if isinstance(item[2], str):
        runString = "sbatch -J " + item[0] + "_" + item[4].replace("hist_", "").replace(".root", "") +\
                    " scripts/job_slurm_process.sh " + item[0] + " " + item[1] + " " +\
                    item[2] + " " + item[3] + " " + item[4] + " " +\
                    workdir + " " +logdir + " " + item[5]
        if len(options.mode) > 0: runString += " " + options.mode
        elif options.ff:          runString += " --ff"
        runString_list.append(runString)
    else:
        for fidx in range(len(item[2])):
            runString = "sbatch -J " + item[0] + "_" + item[4].replace("hist_", "").replace(".root", "_" + str(fidx)) +\
                        " scripts/job_slurm_process.sh " + item[0] + " " + item[1] + " " + item[2][fidx] + " " +\
                        item[3] + " " + item[4].replace(".root", "_" + str(fidx) + ".root")  + " " +\
                        workdir + " " +logdir + " " + item[5]
            if len(options.mode) > 0: runString += " " + options.mode
            elif options.ff: runString += " --ff"
            runString_list.append(runString)

for runString in runString_list:
    print(runString)
    if not options.dry:
        call([runString], shell=True)

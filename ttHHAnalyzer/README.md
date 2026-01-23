# NanoAOD RDataFrame

With this package, we can do the following to NanoAOD files
- Apply good JSON
- Calculate b-tag weights for MC
- Apply Jet/MET corrections
- Skim events
- Process events with applying selections
- SF and uncertainty calculation

For this, we use the data frame concept to simplify the code,
but it requires to learn new way of using ROOT.
In a data frame concept, the operations on data should be defined in terms of functions,
which means that result of applying a function should produce outputs and no internal states are allowed.
In ROOT, RDataFrame allows us to borrow the concept for analysis.

This package code implements LFV analysis using NanoAOD format.


## I. Introduction

- The advantage of using RDataFrame is that you don't have to
worry about looping, filling of histograms at the right place,
writing out trees. In some sense you can concentrate more on thinking about data analysis.

- Draw back is that you have to think of a new way of doing old things.
This could getting used to. Also complicated algorithms may not
be so easy to implement with RDataFrame

- Purely data frame concept is not ideal for HEP data since
we have to treat real data and simulated data differently.
We must apply different operations depending on the input data, such as scale factors and corrections.
Therefore, some form of global state should be stored. 
In this package, the data frame concept is used together with object oriented concept for this purpose.


## II. Code

- The code consists of a main class NanoAODAnalyzerrdframe. 
There is one header file `.h` and one source file `.cpp` for it.
The class has several methods:
    - Object definitions. `selectElectrons()`, `selectMuons()`, ...
    - Additional derived variables definition. `defineMoreVars()`
    - Histogram definitions. `bookHists()`
    - Selections. `defineCuts()`
    - Plus some utility methods/functions. `gen4vec()`, `readjson()`, `removeOverlaps()`, `helper_1DHistCreator()`, `createHists()`

- **Users should modify**: object definitions, define additional variables, histogram definitions, selections.

- To use the class, you need to pass a pointer to TChain object, output root file name and a JSON file.
  If there is `"genWeight"` branch in the ROOT Tree then the input is assumed to be MC, otherwise data.
  Look at `nanoaoddataframe.cpp` to find how to use within C++.

- Compile
Make sure to source proper evnv
``` bash
  source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc12-opt/setup.sh
```
Then compile
``` bash
  make clean && make -j 4   # max 6 effective.
```
Need ROOT 6.26.06 or later then for correctionlib,
this will compile and create a `libnanoaodrdframe.so` shared library that can be loaded in a ROOT session or macro:
```c++
  gSystem->Load("libnanoadrdframe.so");
```
or within pyROOT (look in `processnanoaod.py`).

Now you install plotIt for plotting histograms. Tha package is maintained from the other repo, to be used in the various analyses.
``` bash
# Start from the base directory
cd LFVRun3
git clone https://github.com/yeonsu108/plotIt.git
cd plotIt
cd external
./build-external.sh
cd ..
make -j 4
```

## III. Running over large dataset

### Usage of scripts
In the `scripts/` folder, there are scripts for skimming or processing the NanoAOD files.
All of the options are set in the script files `scripts/*.sh`.

#### Skim
```bash
# CAUTION::::::Create working directory under /data1 to reduce file i/o on disks
# You are at LFVRun3/nanoaodframe/
# First fetch dataset info from google spreadsheet:
# https://docs.google.com/spreadsheets/d/1KNvsRvXi3sgU45T2qOUztFSX3As4elZB325WaPnSkA8/edit#gid=569299692
# The argument will be the name of campaign in 'Campaign Summary' tab
# This will create file lists in data/dataset

python getDatasetInfo.py v2023_BPix
# output: data/dataset/v2023_BPix/dataset.json

# Get file list from DAS
# Set up your vOMS proxy before running this step
source getDatasetDasList.sh v2023_BPix
# output: data/dataset/v2023_BPix/dataset_{dataset_name}.txt
```
### Command Line Options

This script supports the following command line options:

| Option |               | Argument     | Description                                                                                          |
|--------|---------------|--------------|------------------------------------------------------------------------------------------------------|
| `-V`   | `--version`   | *string*     | **Skim version.** Specify the folder name under `/data1/common/skimmed_NanoAOD/` to use for skimming. |
| `-Y`   | `--year`      | *string*     | **Year of data.** Select one of: `v2022`, `v2023`, or `v2023_BPix`.                                   |
| `-C`   | `--ch`        | *string*     | **Channel.** Choose the analysis channel: `muon` or `electron`. *(Default: muon)*                     |
| `-D`   | `--dataset`   | *string(s)*  | **Dataset(s) to process.** Specify one or more dataset folder names (e.g., `TTTo2L2Nu`).              |
| `-N`   | `--name`      | *string*     | **Output file name.** Process only a specific file (e.g., `280000_7316D0F0-4250-7D44-8244-921B41B9C092`). |
| `-F`   | `--dataOrMC`  | *string*     | **Data/MC flag.** Set to `data` or `mc` to process only data or only Monte Carlo samples.             |
|        | `--dry`       | *(flag)*     | **Dry run.** If set, the script will not submit jobs to Slurm; it only prints the commands.           |
|        | `--local`     | *(flag)*     | **Local run.** If set, the script will process local ROOT files instead of submitting jobs.           |

### Example Usage

```bash
# You can add data or mc flag at the end of command if needed

python scripts/skim.py -V skim_test -Y v2023_BPix -C muon

OR

python scripts/skim.py -V skim_test -F mc -Y v2023_BPix -C muon

# To run specific datasets

python scripts/skim.py -V skim_test -F mc -Y v2023_BPix -C muon -D TTtoLNu2Q

# This script uses slurm on htop.
# Please specify which node to EXCLUDE from the batch job, in scripts/job_slurm_skim.sh

#To find failed jobs knowing that `folder_filename.log`
find log/*/*/* | xargs grep runtime_error
find log/*/*/* | xargs grep fault
find log/*/*/* | xargs grep fatal
find log/*/*/* | xargs grep Traceback
find log/*/*/* | xargs grep ERROR
find log/*/*/* | xargs grep error

#Do this ONLY for systematic root file, unless will submit all variations in addition to nominal one
python scripts/skim.py -V skim_test -Y 2018 --dry | grep 270000_221AB515 | sh

#Using commands above, write resubmit list. Create each list and run submit command. Unless, there may be duplicated list.
find log/*/*/* | xargs grep -l runtime_error | sed 's/.log//' | sed 's/log\/*202*\/.*\///' > ~/resub
find log/*/*/* | xargs grep -l fault | sed 's/.log//' | sed 's/log\/*202*\/.*\///' > ~/resub
find log/*/*/* | xargs grep -l fatal | sed 's/.log//' | sed 's/log\/*202*\/.*\///' > ~/resub
find log/*/*/* | xargs grep -l Traceback | sed 's/.log//' | sed 's/log\/*202*\/.*\///' > ~/resub
find log/*/*/* | xargs grep -l error | sed 's/.log//' | sed 's/log\/*202*\/.*\///' > ~/resub
find log/*/*/* | xargs grep ERROR | grep -v "SimpleJetCorrectionUncertainty" | sed 's/.log//' | sed 's/log\/*202*\/.*\///' > ~/resub

#Submit for each year...
cat ~/resub | xargs -i -l1 python scripts/skim.py -V skim_test -Y v2023_BPix -C muon -N 
```

#### Processing
`scripts/process.py` scripts can automatically run over all ROOT files in an input directory.
``` txt
Usage: python scripts/process.py -V skim_test -O test -Y v2023_BPix -S theory

Options:
  -V, --version         Skim version: folder under /data1/common/skimmed_NanoAOD/
  -O, --outdir          Output folder in your working directory
  -Y YEAR, --year=YEAR  Select v2023_BPix
  -S SYST, --syst=SYST  Systematic: 'data' for Data, 'nosyst' for mc without uncertainties. Default is 'theory'. To run without theory unc for TT samples, put 'all'
  -D, --dataset         Put dataset folder name (eg. -D TTto2L2Nu,QCD_Pt1000_MuEnriched) to process specific dataset.
  -F, --dataOrMC        Flag to choose Data or MC.
```

In some cases, you may want to submit single file per core using slurm.
`scripts/process.py` will do the job
``` txt
Usage: python scripts/process.py -V skim_test -O test -Y v2023_BPix -C muon -S nosyst (-F data/mc)

For failing job detection:
find log/* | xargs grep runtime_error
find log/* | xargs grep fault
find log/* | xargs grep fatal
find log/* | xargs grep Traceback
```
To apply Tau Fake Factor with ABCD, simple analyzer is introduced (TauFakeFactorAnalyzer). Output foulder MUST include 'fake'. The l stands for loose tau, t for tight.
``` txt
python scripts/process.py -V skim_v9_230714_jes -O v9_0714_fake_lss -Y 2018 -S nosyst -M lss
python scripts/process.py -V skim_v9_230714_jes -O v9_0714_fake_los -Y 2018 -S nosyst -M los
python scripts/process.py -V skim_v9_230714_jes -O v9_0714_fake_tss -Y 2018 -S nosyst -M tss
python scripts/process.py -V skim_v9_230714_jes -O v9_0714_fake_tos -Y 2018 -S nosyst -M tos

# Run simple script for FF calculation. In the 'fromROOTFile', one can calculate per tau_pt bin
python fake_factor_calculator.py
# or
python fake_factor_calculator_fromROOTFile.py

python scripts/process.py -V reweight_test -O v9_0714_1010_FF -Y 2018 -S theory --ff
```
Now, apply b SF rescaling, compute uncertainty envelope, etc. Let the `test/2018` is the folder containing histograms.
``` txt
python postprocess.py -I test -Y 2018
```
Drawing histogram by plotIt
``` txt
cd v2023_BPix
mkdir ../figure_v2023BPix
hadd hist_Muon.root hist_Muon*.root
../../../../plotIt/plotIt -o ../figure_v2023BPix/ ../../../../plotIt/configs/Run3_muon/config_2023_BPix.yml -y -s
# Add -q to draw with QCD
# Add --allSig to draw all signals for plot ir yield
cd ../../
python print_syst_table.py -I test -Y v2023_BPix

#For Run2,
python stack_signals_v2.py -I test #use v2 for now, v1 cannot deal with year based uncertainties
python plot_run2.py -I test
```

### Usage of scripts
In the `scripts/` folder, there are scripts for skimming or processing the NanoAOD files.
All of the options are set in the script files `scripts/*.sh`.

# NanoAOD RDataFrame — LFV Top Analysis

This package implements an LFV (Lepton Flavour Violation) search in top quark decays
using NanoAOD files processed with ROOT RDataFrame.

**With this package, we can do the following to NanoAOD files**
- Apply golden JSON for luminosity masking (data)
- Apply JEC/JER, lepton, tau, and b-tag scale factors for MC
- Compute systematic uncertainty variations
- Skim events with loose object preselection and branch reduction
- Process skimmed events with full analysis selections, corrections, and event weighting
- Produce flat ntuples and histograms for statistical analysis and DNN training


## I. Framework Architecture

The data flow is:

```
NanoAOD (CMS central)
     │
     ▼  SkimEvents (src/SkimEvents.cpp)
Skimmed NanoAOD    ←  stored under $LFV_SKIM_DIR/{version}/{ch}/{data|mc}/{year}/
     │
     ▼  TopLFVAnalyzer (src/TopLFVAnalyzer.cpp)
Flat Ntuple + Histograms  ←  stored in the process output directory
     │
     ├──▶  postprocess.py          (b-SF rescaling, histogram merging)
     ├──▶  DNN/eval_multi.py       (DNN training / inference)
     └──▶  plotIt / print_syst_table.py  (plotting, systematic table)
```

### Skim stage (`SkimEvents`)

Runs once per NanoAOD file. Its job is:

- Apply a golden JSON lumi mask (data only)
- Apply JEC (L1FastJet + L2Relative + L3Absolute + L2L3Residual), JER smearing, and propagation to MET
- Store pile-up weights (`puWeight[0/1/2]`), lepton SFs (ID/Iso/Trg), b-tag shape SFs (`btagWeight[0..16]`), tau ID SFs, and top-pT reweighting — all as branches in the skimmed NanoAOD
- Apply a **loose** object and event preselection (≥1 lepton, ≥1 tau, ≥1 jet)
- Drop all unnecessary NanoAOD branches; keep only what the process stage needs

The skim output is still in NanoAOD format (`TTree Events`), just smaller.

### Process stage (`TopLFVAnalyzer`)

Runs per dataset (or per file for shape systematics). Does the full analysis:

- **Object selection** — full analysis cuts from `data/config/analysis_config.json`; for shape systematics (JES/JER/TES/MET Unclust/muon GE scale), the primary object is varied before all downstream quantities are computed
- **Kinematic variables** — 4-vectors, ΔR, M_T, S_T^MET, top reconstruction χ², lepton+τ mass, object-level branches (`Jet1_pt`, `Tau1_pt`, etc.)
- **Event weights** — nominal weight and all systematic variations stored as separate branches (`eventWeight`, `eventWeight__puup`, `eventWeight__btaghfup`, …)
- **b-tag normalization** — per-njet-bin normalization pass before writing weights, following CMS BTV recommendation
- **Tau fake factor** — applied optionally via `--ff` flag
- **Histogram booking** — 1-D distributions at each selection step, with all weight variations
- **Flat ntuple writing** — `Events` TTree with all output branches for DNN and limit setting

The process stage output is a standard ROOT file with `Events` tree + histograms.


## II. Code Structure

```
LFVAnalyzer/
├── src/
│   ├── SkimEvents.{h,cpp}           # Skim stage
│   ├── TopLFVAnalyzer.{h,cpp}       # Process stage (analysis + weights + ntuple)
│   ├── TauFakeFactorAnalyzer.{h,cpp} # Fake-rate measurement mode
│   └── nanoaodrdataframe.cpp        # Entry point for standalone binary
├── scripts/
│   ├── skim.py                      # Batch submission for skim
│   ├── skimonefile.py               # Single-file skim worker
│   ├── process.py                   # Batch submission for process
│   ├── processonefile.py            # Single-file process worker
│   ├── processonedataset.py         # Per-dataset process worker
│   ├── job_slurm_skim.sh            # Slurm wrapper for skim
│   ├── job_slurm_process.sh         # Slurm wrapper for process
│   ├── checkerror.py                # Scan log directory for errors
│   ├── isZoombie.py                 # Scan output for zombie/incomplete ROOT files
│   └── validate_systematics.py      # Systematic propagation validation
├── data/
│   ├── config/analysis_config.json  # Object cuts, tau WPs, b-tag WPs
│   ├── dataset/{year}/              # File lists (see Section IV)
│   ├── GoldenJSON/                  # CMS golden JSON per era
│   ├── BTV/                         # b-tag correctionlib JSONs
│   ├── JME/                         # JEC/JER correctionlib JSONs
│   ├── MuonSF/, ElectronSF/, TauIDSFs/ # Lepton correctionlib JSONs
│   ├── LUM/                         # Pileup reweighting JSONs
│   └── TauFF/                       # Tau fake factors (JSON)
├── postprocess.py                   # b-SF rescaling and histogram merging
├── print_syst_table.py              # Systematic uncertainty table
├── fake_factor_calculator*.py       # Tau fake rate calculation helpers
├── Makefile
└── README.md
```

`../CommonTools/src/` contains the base class `NanoAODAnalyzerrdframe` shared with the skim and process stages.


## III. Setup

### 1. Environment

Source the LCG environment (EL9 node required for Run 3):

```bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_108/x86_64-el9-gcc15-opt/setup.sh
```

> **Note:** LCG_108 is the current default used by both skim and process Slurm scripts.  
> LCG_103 (CentOS7) and LCG_105 (EL9 early) are kept as comments for Run 2 UL jobs.

`correctionlib` is bundled with LCG_108. The Makefile auto-detects its headers and library via `correction config`.

### 2. Compile

```bash
make clean && make -j 4
```

This produces `libnanoadrdframe.so`, loaded at runtime by the Python process scripts via `cppyy`.

Debug build:
```bash
make DEBUG=1 -j 4
```

### 3. Optional: storage path overrides

By default the code looks for NanoAOD and skimmed files at:

| Purpose | Default path | Override env var |
|---------|-------------|-----------------|
| Skimmed NanoAOD output/input | `/data2/common/skimmed_NanoAOD` | `LFV_SKIM_DIR` |
| Raw NanoAOD (local mode) | `/data2/common/NanoAOD` | `LFV_NANO_DIR` |

```bash
export LFV_SKIM_DIR=/my/skim/storage
export LFV_NANO_DIR=/my/nanoaod/storage
```

### 4. Install plotIt (optional, for plotting)

```bash
cd ..              # go to NanoAODRun3/
git clone https://github.com/yeonsu108/plotIt.git
cd plotIt/external && ./build-external.sh
cd .. && make -j 4
```


## IV. Dataset Preparation

### File lists

The skim script reads file lists from `data/dataset/{year}/dataset_{DatasetName}.txt`.  
Each line contains a DAS path (e.g. `/store/mc/Run3Summer22NanoAODv12/...`).

Fetch the list from the analysis spreadsheet:

```bash
python getDatasetInfo.py v12_2023BPix
# output: data/dataset/v12_2023BPix/dataset.json
```

Populate DAS paths (requires a valid VOMS proxy):

```bash
voms-proxy-init -voms cms
source getDatasetDasList.sh v12_2023BPix
# output: data/dataset/v12_2023BPix/dataset_{DatasetName}.txt
```

Supported era keys: `v12_2022`, `v12_2022EE`, `v12_2023`, `v12_2023BPix`, `v15_2024`.

### Analysis configuration

Object cuts, tau DeepTau WPs, and b-tag WP/tagger are defined in `data/config/analysis_config.json`:

```json
{
  "object_selection": {
    "muon_veto_cut":     "Muon_pt>15.0 && abs(Muon_eta)<2.4 && Muon_looseId && ...",
    "electron_veto_cut": "Electron_pt>15.0 && abs(Electron_eta)<2.5 && Electron_cutBased == 1",
    "jet_cut":           "Jet_pt>40.0 && abs(Jet_eta)<2.5 && Jet_passJetIdTightLepVeto==1.0 ...",
    "tau_cut":           "Tau_pt>40.0 && abs(Tau_eta)<2.5 && Tau_idDecayModeNewDMs ..."
  },
  "deeptau_wp": {
    "muon_ch":     { "vsjet": "7", "vsmu": "4", "vse": "2" },
    "electron_ch": { "vsjet": "7", "vsmu": "1", "vse": "6" }
  },
  "btag": {
    "wp_medium": 0.1272,
    "tagger_run3": "particleNet_shape",
    "tagger_2024": "UParTAK4_comb"
  }
}
```

> Modifying cuts in this file affects both skim and process stages.


## V. Running the Skim

### Options

| Option | | Argument | Description |
|--------|--|---------|-------------|
| `-V` | `--version` | *string* | Skim version tag — output folder under `$LFV_SKIM_DIR/` |
| `-Y` | `--year` | *string* | Era: `v12_2022`, `v12_2022EE`, `v12_2023`, `v12_2023BPix`, `v15_2024` |
| `-C` | `--ch` | *string* | Channel: `muon` or `electron` (default: `muon`) |
| `-D` | `--dataset` | *string(s)* | Restrict to specific dataset name(s) |
| `-N` | `--name` | *string* | Restrict to a single output filename (e.g. `280000_7316D0F0`) |
| `-F` | `--dataOrMC` | *string* | `data` or `mc` to skip the other |
| `-P` | `--path` | *string* | XRootD prefix (default: `root://xrootd-cms.infn.it/`; use `fnal` for FNAL) |
|  | `--dry` | flag | Print commands only — do not submit to Slurm |
|  | `--local` | flag | Read from `$LFV_NANO_DIR` instead of XRootD |

### Example usage

```bash
# Skim all datasets for the 2023BPix era, muon channel (batch)
python scripts/skim.py -V skim_v1 -Y v12_2023BPix -C muon

# Skim MC only
python scripts/skim.py -V skim_v1 -F mc -Y v12_2023BPix -C muon

# Skim a specific dataset
python scripts/skim.py -V skim_v1 -Y v12_2023BPix -C muon -D TTtoLNu2Q

# Dry-run to inspect commands
python scripts/skim.py -V skim_v1 -Y v12_2023BPix -C muon --dry

# Read local files instead of XRootD
python scripts/skim.py -V skim_v1 -Y v12_2023BPix -C muon --local
```

### Skim output structure

```
$LFV_SKIM_DIR/{version}/{ch}/
├── data/{year}/{DatasetName}/{dirNum}_{filename}.root
├── mc/{year}/{DatasetName}/{dirNum}_{filename}.root
└── log/{year}/{DatasetName}/{dirNum}_{filename}.log
```

Each output file is a NanoAOD-format ROOT file (`Events` tree) with:
- All original branches needed downstream
- Added weight branches: `puWeight`, `btagWeight`, `muonWeightId/Iso/Trg`, `elecWeightReco/Id/Trg`, `tauWeightIdVsJet/VsEl/VsMu`, `TopPtWeight`
- JEC/JER already applied to jet and MET branches

### Check for failed skim jobs

```bash
# Automatic log scan
python scripts/checkerror.py /target/path
# Results in error_logs.txt and error_but_done_logs.txt

# Manual grep (any error keyword)
find log/*/*/*  | xargs grep -l "runtime_error\|fault\|fatal\|Traceback\|ERROR" > ~/resub

# Check for zombie or incomplete output files
python scripts/isZoombie.py $LFV_SKIM_DIR/skim_v1/muon/mc/v12_2023BPix --isSkim

# Resubmit individual files from resub list
cat ~/resub | xargs -I{} python scripts/skim.py -V skim_v1 -Y v12_2023BPix -C muon -N {}
```

### Adjust Slurm exclusion list

Edit `scripts/job_slurm_skim.sh` to exclude problematic nodes:
```bash
#SBATCH -p high_cpu,cpu,gpu -x gpu-0-1
```


## VI. Running the Process

### Options

| Option | | Argument | Description |
|--------|--|---------|-------------|
| `-V` | `--version` | *string* | Skim version tag to read from `$LFV_SKIM_DIR/` |
| `-O` | `--outdir` | *string* | Output directory under the working directory |
| `-Y` | `--year` | *string* | Era tag (same values as skim) |
| `-C` | `--ch` | *string* | `muon` or `electron` |
| `-S` | `--syst` | *string* | Systematic mode (see below) |
| `-D` | `--dataset` | *string(s)* | Restrict to specific dataset(s) |
| `-F` | `--dataOrMC` | *string* | `data` or `mc` |
| `-M` | `--mode` | *string* | Fake-rate mode: `lss`, `los`, `tss`, `tos` |
|  | `--ff` | flag | Apply tau fake factor (output folder name must contain `FF`) |
|  | `--dry` | flag | Print commands only — do not submit |
|  | `--split` | flag | Submit one Slurm job per input ROOT file instead of per dataset |

### Systematic modes (`-S`)

| Value | Meaning |
|-------|---------|
| `data` | Data: no weights, no systematics, golden JSON applied |
| `nosyst` | MC nominal only — no shape systs submitted, no weight syst branches |
| `all` | MC with all weight syst branches (no theory weights for TTbar) |
| `theory` | MC with all weight syst branches **plus** ME/PS/PDF branches for TTbar and signal |

Shape systematics (JES, JER, TES, MET unclust, muon GE scale) are **always** submitted as separate jobs for each `--syst` mode except `nosyst`. Each shape syst produces a separate output file:

```
hist_{DatasetName}__{syst}.root    # e.g. hist_TTto2L2Nu__jesAbsoluteup.root
```

### Example usage

```bash
# Process all 2023BPix muon MC with full theory uncertainties (batch)
python scripts/process.py -V skim_v1 -O process_v1 -Y v12_2023BPix -C muon -S theory

# Process data only, no systematics
python scripts/process.py -V skim_v1 -O process_v1 -Y v12_2023BPix -C muon -S data -F data

# Process specific dataset, all weight systs but no theory weights
python scripts/process.py -V skim_v1 -O process_v1 -Y v12_2023BPix -C muon -S all -D WJetsToLNu

# Quick nominal-only test
python scripts/process.py -V skim_v1 -O test_out -Y v12_2023BPix -C muon -S nosyst -D TTto2L2Nu

# Run a single file directly (no Slurm)
python scripts/processonefile.py \
    -Y v12_2023BPix -C muon -S all \
    -I /path/to/skimmed.root \
    -O /path/to/output/hist_TTto2L2Nu.root
```

### Process output structure

```
{outdir}/{ch}/{year}/
├── hist_{DatasetName}.root                  # nominal
├── hist_{DatasetName}__jesAbsoluteup.root   # JES shape syst
├── hist_{DatasetName}__tesup.root           # TES shape syst
├── hist_{DatasetName}__metUnclustup.root    # MET unclust shape syst
├── ...
└── log/{DatasetName}_{syst}.log
```

Each ROOT file contains:
- **`Events` TTree** — flat ntuple with all output branches (see Section VIII)
- **1-D histograms** — `h_tau1_pt`, `h_leptau_mass`, `h_chi2`, `h_st_met`, … at each selection step
- **`hcounter`** — sum of generator weights for normalisation
- **`LHEPdfWeightSum`, `PSWeightSum`, `ScaleWeightSum`** — theory weight normalisations (theory mode only)

### Check for failed process jobs

```bash
# Automatic log scan (run from process output directory)
python scripts/checkerror.py /target/path

# Check for zombie/incomplete output files
python scripts/isZoombie.py process_v1/muon/v12_2023BPix/

# Find logs with errors
find log/ | xargs grep -l "runtime_error\|fault\|Traceback" > ~/resub_process

# Validate systematic propagation (requires uproot)
python scripts/validate_systematics.py \
    --indir process_v1/muon/v12_2023BPix/ \
    --sample TTto2L2Nu \
    --year v12_2023BPix \
    --ch muon \
    --nevents 10000
```


## VII. Systematic Uncertainties

### Shape systematics (separate output files)

These shift object kinematics and are run as separate jobs. All downstream derived variables (`chi2`, `leptau_mass`, `st_met`, `lepMET_mt`, etc.) and DNN input features are automatically shifted.

| Systematic | Objects varied | Era |
|-----------|---------------|-----|
| `jesAbsoluteup/down` | Jet pT, mass, MET | All |
| `jesAbsolute_{year}up/down` | Same, year-specific | All |
| `jesBBEC1up/down` | Jet pT, mass, MET | All |
| `jesBBEC1_{year}up/down` | Same, year-specific | All |
| `jesRelativeBalup/down` | Jet pT, mass, MET | All |
| `jesRelativeSample_{year}up/down` | Same | All |
| `jesFlavorPure{Gluon,Quark,Charm,Bottom}up/down` | Jet pT, mass, MET | All |
| `jesHEMup/down` | Jet pT, mass, MET | 2018 only |
| `jerup/down` | Jet pT, mass, MET (smearing) | All |
| `tesup/down` | Tau pT, mass | All |
| `metUnclustup/down` | PuppiMET pT, φ | All |
| `muonhighscaleup/down` | Muon pT, MET | Muon channel |

Shape systematics for `hdamp` and `TuneCP5` are handled via dedicated MC datasets (separate dataset entries), not as varied jobs.

### Weight systematics (branches in the nominal file)

For `--syst all` or `--syst theory`, the following branches are added to each MC event alongside `eventWeight`:

| Branch pattern | Physics source |
|---------------|---------------|
| `eventWeight__puup/down` | Pileup reweighting |
| `eventWeight__topptup/down` | Top pT reweighting |
| `eventWeight__muidup/down` | Muon ID SF |
| `eventWeight__muisoup/down` | Muon isolation SF |
| `eventWeight__mutrgup/down` | Muon trigger SF |
| `eventWeight__btagcorrup/down` | b-tag correlated uncertainty |
| `eventWeight__btaguncorrup/down` | b-tag uncorrelated uncertainty |
| `eventWeight__btagstatup/down` | b-tag statistical uncertainty |
| `eventWeight__btagtype3up/down` | b-tag type3 uncertainty |
| `eventWeight__btagbfragup/down` | b-tag b-fragmentation uncertainty |
| `eventWeight__tauidjetUncert{0,1}up/down` | Tau VSjet stat |
| `eventWeight__tauidjetSyst*up/down` | Tau VSjet syst (era/DM) |
| `eventWeight__tauidelup/down` | Tau VSe SF |
| `eventWeight__tauidmuup/down` | Tau VSmu SF |
| `eventWeight__mescaleup/down` | ME scale (theory mode) |
| `eventWeight__renscaleup/down` | Renorm. scale (theory mode) |
| `eventWeight__facscaleup/down` | Factor. scale (theory mode) |
| `eventWeight__isrup/down` | ISR (theory mode) |
| `eventWeight__fsrup/down` | FSR (theory mode) |
| `eventWeight__pdf{i}` | PDF eigenvectors (theory mode) |

**b-tag SF:** Fixed Working Point (Medium) SFs are applied directly to WP-passing (tagged) jets via `btagWeight[0]`, with 10 variation components (`btagWeight[1..10]`). Normalization is not required for fixed-WP SFs.

### DNN input features

The following branches are used by `DNN/train_multi.py`:

```python
# Muon channel (28 features)
["Muon1_pt", "Muon1_eta",
 "Tau1_pt", "Tau1_mass", "Tau1_eta",
 "Jet1_pt", "Jet1_mass", "Jet1_eta", "Jet1_btagPNetB",
 "Jet2_pt", "Jet2_mass", "Jet2_eta", "Jet2_btagPNetB",
 "Jet3_pt", "Jet3_mass", "Jet3_eta", "Jet3_btagPNetB",
 "chi2", "chi2_SMW_mass", "chi2_SMTop_mass",
 "chi2_wqq_dEta", "chi2_wqq_dPhi", "chi2_wqq_dR",
 "leptau_mass", "leptau_dEta", "leptau_dPhi", "leptau_dR",
 "PuppiMET_pt"]

# Electron channel: replace Muon1_pt/eta with Electron1_pt/eta
```


## IX. Postprocessing and Plotting

### b-SF rescaling and histogram merging

```bash
python postprocess.py -I process_v1 -Y v12_2023BPix
```

Options:
- `-I` : process output directory (same as `-O` argument used in `process.py`)
- `-Y` : year tag; supported: `v12_2022`, `v12_2022EE`, `v12_2023`, `v12_2023BPix`, `v15_2024`
- `--postfix` : suffix for output folder (e.g. for alternative binning)
- `-F` / `--forceHadd` : force re-hadding split files

Output goes to `{indir}/{year}_postprocess/`.

### Plotting with plotIt

```bash
cd process_v1/muon/v12_2023BPix
mkdir ../v2023BPix
hadd hist_Muon.root hist_Muon*.root
../../../../plotIt/plotIt -o ../figure_v2023BPix/ \
    ../../../../plotIt/configs/Run3_muon/config_2023BPix.yml -y -s
```

### Print systematic uncertainty table

```bash
python print_syst_table.py -I process_v1 -Y v12_2023BPix
```

### Tau fake factor workflow

```bash
# Step 1: run in fake-rate measurement mode (output folder must contain 'fake')
python scripts/process.py -V skim_v1 -O process_v1_fake_lss -Y v12_2023BPix -C muon -S nosyst -M lss
python scripts/process.py -V skim_v1 -O process_v1_fake_los -Y v12_2023BPix -C muon -S nosyst -M los
python scripts/process.py -V skim_v1 -O process_v1_fake_tss -Y v12_2023BPix -C muon -S nosyst -M tss
python scripts/process.py -V skim_v1 -O process_v1_fake_tos -Y v12_2023BPix -C muon -S nosyst -M tos

# Step 2: compute fake factors
python fake_factor_calculator_fromROOTFile.py   # per-pt-bin version

# Step 3: apply in final selection (output folder must contain 'FF')
python scripts/process.py -V skim_v1 -O process_v1_FF -Y v12_2023BPix -C muon -S theory --ff
```

### DNN training and inference

```bash
cd ../DNN
# Train
python train_multi.py ...

# Inference / evaluation
python eval_multi.py \
    -O DNN_eval \
    -I top_lfv_<timestamp> \
    -C muon \
    -P /path/to/process_v1/muon/ \
    --xsecfile xsec.yaml \
    --alpha 0.1
```


## X. Validation

### Systematic propagation check

```bash
python scripts/validate_systematics.py \
    --indir process_v1/muon/v2023_BPix/ \
    --sample TTto2L2Nu \
    --year v12_2023BPix \
    --ch muon \
    --nevents 10000 \
    --outfile syst_report.json
```

The script runs 7 categories of automated checks and exits with code 1 if any fail. See `scripts/validate_systematics.py` for full documentation.

### Output file integrity

```bash
# Check for zombie or incomplete skim files
python scripts/isZoombie.py $LFV_SKIM_DIR/skim_v1/muon/mc/v2023_BPix --isSkim

# Check for zombie or incomplete process files
python scripts/isZoombie.py process_v1/muon/v2023_BPix/
```

A valid skim file has: `Events`, `hcounter`, `hcounter_S1`, `LHEPdfWeightSum`, `PSWeightSum`, `ScaleWeightSum`.
A valid process file has: `Events`, `h_tau1_pt_S5` (i.e. histograms at the final selection step).


## XI. Troubleshooting

### Environment / compile

| Problem | Fix |
|---------|-----|
| `TTree.h` not found | Source LCG: `source /cvmfs/sft.cern.ch/.../LCG_108/.../setup.sh` |
| `correction` not found after sourcing | Re-source LCG; `correctionlib` is only available from LCG_104+ |
| `rootcint` fails | Check `rootcling` is on PATH; use the same ROOT version for compile and runtime |
| `libnanoadrdframe.so` not found at runtime | Run from the `LFVAnalyzer/` directory; the `.so` must be in `$PWD` |

### Skim stage

| Problem | Fix |
|---------|-----|
| XRootD timeout / stall | Use `--local` if files are available locally; or try `--path fnal` |
| Missing `data/dataset/{year}/` files | Run `getDatasetInfo.py` and `getDatasetDasList.sh` first |
| `Error in TFile::Open` on data | Check VOMS proxy: `voms-proxy-init -voms cms` |
| All skim jobs exit immediately | Check Slurm node availability; check `#SBATCH -x` exclusion list |

### Process stage

| Problem | Fix |
|---------|-----|
| `No such column: btagWeight` | Skim was run without b-tag SFs; re-skim with current code |
| `hcounter not found` | Input skim file did not complete; re-skim |
| `There is NO EVENT to process` | Input skim file is empty after preselection; expected for some datasets |
| Process jobs finish very quickly with 0 events | Check log for "Zombie" or that the skimmed files exist in `$LFV_SKIM_DIR` |
| Shape syst files missing | `--syst nosyst` was used; re-run without `nosyst` |
| `postprocess.py` year not recognised | Only `v2022`, `v2022EE`, `v2023`, `v2023_BPix` are supported; 2024 not yet |

### Log scanning

```bash
# Automated (recommended)
python scripts/checkerror.py          # from the log directory

# Manual
find log/ | xargs grep -l "runtime_error\|fault\|fatal\|Traceback\|ERROR" > ~/errors.txt

# Errors that are benign (known false positives)
# "SimpleJetCorrectionUncertainty" messages in logs are informational, not errors
find log/ | xargs grep ERROR | grep -v SimpleJetCorrectionUncertainty > ~/real_errors.txt
```

# NanoAOD RDataFrame — LFV Top Analysis

This package implements an LFV (Lepton Flavour Violation) search in top quark decays
using NanoAOD files processed with ROOT RDataFrame.

**With this package, we can do the following to NanoAOD files**
- Apply golden JSON for luminosity masking (data, enforced at the Process stage)
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
     ├──▶  postprocess.py          (theory scale/PS/PDF renormalization, external syst rescaling, histogram merging)
     ├──▶  DNN/eval_multi.py       (DNN training / inference)
     └──▶  plotIt / print_syst_table.py  (plotting, systematic table)
```

### Skim stage (`SkimEvents`)

Runs once per NanoAOD file. Its job is:

- **No Golden JSON at this stage** — the lumi mask is *not* applied in the skim (`skim.py` does not pass `-J`); the Golden JSON is applied and enforced in the Process stage (see below), so the skim keeps raw lumisections and does not need to be redone when a Golden JSON is updated
- Apply JEC (L1FastJet + L2Relative + L3Absolute + L2L3Residual), JER smearing, and a **full Type-1 PuppiMET/PFMET recomputation from `RawPuppiMET`/`RawMET`** (including `CorrT1METJet`; MET XY-shift corrections for 2022–2023BPix are applied afterwards)
- Apply the nominal **muon ScaReKIT** and **electron EGM scale & smearing before object selection** (`Muon_pt_uncor`, `Electron_pt_uncor` are preserved; only the *selected* leptons shift PuppiMET); lepton SFs are evaluated after the selection
- Store pile-up weights (`puWeight[0/1/2]`; 2024 uses `Collisions24_CDEFGHI_goldenJSON`), lepton SFs (ID/Iso/Trg), tau energy scale (`Tau_ES_nom/up/down`) & ID SFs, top-pT reweighting, and pre-selection theory weight sums (`ScaleWeightSum`, `PSWeightSum`, `LHEPdfWeightSum`, accumulated in a single event loop) — all in the skimmed NanoAOD
- Apply a **loose** object and event preselection (≥1 lepton, ≥1 tau, ≥1 jet); jets are kept down to $p_T>15$ GeV, $|\eta|<5.2$ so that systematic MET propagation sees all jets
- Drop all unnecessary NanoAOD branches; keep only what the process stage needs (including `Muon_nTrackerLayers`, `Electron_seedGain`, `Electron_r9`, `Electron_deltaEtaSC`)

The skim output is still in NanoAOD format (`TTree Events`), just smaller.

### Process stage (`TopLFVAnalyzer`)

Runs per dataset (or per file for shape systematics). Does the full analysis:

- **Object selection** — full analysis cuts from `data/config/analysis_config.json`; for shape systematics (JES / JER / TES / MET Unclust / Muon ScaReKIT `muonscale`,`muonsmear` / Electron S&S `elescale`,`elesmear`), the primary object is varied at the 4-vector level and propagated to Type-1 PuppiMET before all downstream quantities are computed
- **Golden JSON (Data)** — `processonedataset.py` / `processonefile.py` resolve the Golden JSON from `--year` (fallback: input path) when `-S data` is used; if the era cannot be resolved or the file is missing they print `[ERROR]` and exit. `readjson()` throws if the file cannot be opened, and `TopLFVAnalyzer::defineObjectSelection()` throws if Data is processed without a loaded Golden JSON
- **Kinematic variables** — 4-vectors, ΔR, M_T, S_T^MET, top reconstruction χ², lepton+τ+b mass, scalar `Jet_HT`, object-level branches (`Jet1_pt`, `bJet1_pt`, `Tau1_pt`, etc.)
- **Event weights** — nominal weight and all systematic variations stored as separate branches (`eventWeight`, `eventWeight__puup`, `eventWeight__btagSFbc_correlatedup`, `eventWeight__btagSFbc_<era>up`, `eventWeight__btagSFlight_correlatedup`, `eventWeight__btagSFlight_<era>up`, …)
- **b-tag Fixed-WP SF (Method 1a)** — event-by-event reweighting over tagged and untagged jets using 2D MC efficiency maps (`data/BTV/<ch>/<era>/btag_eff_*.root`) and `btagging.json.gz`, with heavy-flavour (`bc`) and light-flavour (`light`) correlated/uncorrelated nuisances (`btagWeight[0..12]`, `btagWeight_loose[0..12]`)
- **Tau fake factor** — applied optionally via `--ff` flag (and measured via `TauFakeFactorAnalyzer` with Medium b-tag WP and `btagWeight_loose[0]`)
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
    "muon_veto_cut":     "Muon_pt>15.0 && abs(Muon_eta)<2.4 && Muon_looseId && Muon_pfRelIso04_all<0.25",
    "electron_veto_cut": "Electron_pt>15.0 && abs(Electron_eta)<2.5 && !(abs(Electron_eta + Electron_deltaEtaSC)>1.4442 && abs(Electron_eta + Electron_deltaEtaSC)<1.566) && Electron_cutBased >= 1",
    "jet_cut":           "Jet_pt>40.0 && abs(Jet_eta)<2.5 && Jet_passJetIdTightLepVeto==1.0 ...",
    "tau_cut":           "Tau_pt>40.0 && abs(Tau_eta)<2.5 && abs(Tau_dz)<0.2 && Tau_idDecayModeNewDMs ..."
  },
  "deeptau_wp": {
    "muon_ch":     { "vsjet": "7", "vsmu": "4", "vse": "2" },
    "electron_ch": { "vsjet": "7", "vsmu": "1", "vse": "6" }
  },
  "btag": {
    "wp_medium": { "2022_Summer22": 0.2450, "2022_Summer22EE": 0.2605,
                   "2023_Summer23": 0.1917, "2023_Summer23BPix": 0.1919,
                   "2024_Summer24": 0.1272 },
    "tagger_run3": "particleNet_comb",  "tagger_light_run3": "particleNet_light",
    "tagger_2024": "UParTAK4_comb",     "tagger_light_2024": "UParTAK4_light"
  }
}
```

> Modifying cuts in this file affects both skim and process stages. The b-tag Medium WP (`btag.wp_medium`) is the **single source** for `selectJets`, `btagSF_fixedWP`, `TopLFVAnalyzer` and `BTagEfficiencyAnalyzer` (strict `score > wp`); `applyBSFs` cross-checks it against `*_wp_values("M")` in `btagging.json.gz` and aborts with `[ERROR]` on a mismatch. DeepTau WP integers are 1-based (`VSjet/VSe`: 1=VVVLoose … 7=VTight, 8=VVTight; `VSmu`: 1=VLoose … 4=Tight).


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
|  | `--local` | flag | Read from `$LFV_NANO_DIR/{data,mc}/<nano_ver>/…` instead of XRootD, where `<nano_ver>` (`v12` or `v15`) is the prefix of `--year` (e.g. `v12_2023BPix` → `v12`). *(The `v12` sub-folder layout is assumed, not verified.)* |

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
| `data` | Data: no weights, no systematics; Golden JSON is **required** (resolved from `--year`, `[ERROR]`+exit if unresolved/missing). Data vs MC is decided by which directory the dataset was listed from (`data/` vs `mc/`), not by a string match on the path |
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

# Validate systematic propagation (requires uproot; pass plain era to --year for JES checks)
python scripts/validate_systematics.py \
    --indir process_v1/muon/v12_2023BPix/ \
    --sample TTto2L2Nu \
    --year 2023BPix \
    --ch muon \
    --nevents 10000
```


## VII. Systematic Uncertainties

### Shape systematics (separate output files)

These shift object kinematics and are run as separate jobs. All downstream derived variables (`chi2`, `leptau_mass`, `st_met`, `lepMET_mt`, etc.) and DNN input features are automatically shifted.

| Systematic | Objects varied | Era / Channel |
|-----------|---------------|-----|
| `jesAbsoluteup/down` | Jet pT, mass, Type-1 PuppiMET | All (Correlated) |
| `jesAbsolute_{era}up/down` | Jet pT, mass, Type-1 PuppiMET | Per era (`2022`, `2022EE`, `2023`, `2023BPix`, `2024` — sub-era names, `jes_year = year`) |
| `jesBBEC1up/down` | Jet pT, mass, Type-1 PuppiMET | All (Correlated) |
| `jesBBEC1_{era}up/down` | Jet pT, mass, Type-1 PuppiMET | Per era (sub-era names) |
| `jesEC2up/down`, `jesEC2_{era}up/down` | Jet pT, mass, Type-1 PuppiMET | All / Per era |
| `jesHFup/down`, `jesHF_{era}up/down` | Jet pT, mass, Type-1 PuppiMET | All / Per era |
| `jesRelativeBalup/down` | Jet pT, mass, Type-1 PuppiMET | All (Correlated) |
| `jesRelativeSample_{era}up/down` | Jet pT, mass, Type-1 PuppiMET | Per era (sub-era names) |
| `jesFlavorPure{Gluon,Quark,Charm,Bottom}up/down` | Jet pT, mass (by parton flavour), Type-1 PuppiMET | All (Correlated) |
| `jerup/down` | Jet pT, mass ($c_{\text{var}}/c_{\text{nom}}$ ratio), Type-1 PuppiMET | All (Uncorrelated per era) |
| `tesup/down` | Tau pT, mass (`Tau_ES_up/down / Tau_ES_nom`) | All (Uncorrelated per era) |
| `metUnclustup/down` | PuppiMET pT, φ (2D shift preserving XY-correction) | All (Uncorrelated per era) |
| `muonscaleup/down`, `muonsmearup/down` | Muon pT, mass (ScaReKIT), Type-1 PuppiMET | Muon channel (`scale` Corr, `smear` Uncorr) |
| `elescaleup/down`, `elesmearup/down` | Electron pT, mass (EGM S&S), Type-1 PuppiMET | Electron channel (`scale` Corr, `smear` Uncorr) |

Shape systematics for `hdamp` and `TuneCP5` (`tune`) are handled via dedicated SM $t\bar{t}$ MC datasets (`-S nosyst` + `postprocess.py` `hcounter` rescaling).

### Weight systematics (branches in the nominal file)

For `--syst all` or `--syst theory`, the following branches are added to each MC event alongside `eventWeight`:

| Branch pattern | Physics source |
|---------------|---------------|
| `eventWeight__puup/down` | Pileup reweighting (Correlated across Run 3) |
| `eventWeight__topptup/down` | Top pT reweighting (Correlated across Run 3) |
| `eventWeight__mu{id,iso,trg}_{stat,syst}up/down` | Muon ID / Iso / Trigger SF (`stat` Uncorr per era, `syst` Corr). The unsplit `mu{id,iso,trg}up/down` are **no longer booked** (Combine still discards them if found in old files) |
| `eventWeight__elec{id,reco}_{stat,syst}up/down`, `eventWeight__electrgup/down` | Electron ID / Reco SF (`stat` Uncorr, `syst` Corr) and Trigger SF (total only, per-era Uncorr). The unsplit `elec{id,reco}up/down` are **no longer booked** |
| `eventWeight__btagSFbc_correlatedup/down` | Heavy-flavour ($b,c$) b-tag SF correlated uncertainty across Run 3 (`btagWeight[5/6]`) |
| `eventWeight__btagSFbc_{era}up/down` | Heavy-flavour ($b,c$) b-tag SF uncorrelated uncertainty for `{era}` $\in \{\text{2022, 2022EE, 2023, 2023BPix, 2024}\}$ (`btagWeight[7/8]`) |
| `eventWeight__btagSFlight_correlatedup/down` | Light-flavour ($u,d,s,g$) b-tag SF correlated uncertainty across Run 3 (`btagWeight[9/10]`) |
| `eventWeight__btagSFlight_{era}up/down` | Light-flavour ($u,d,s,g$) b-tag SF uncorrelated uncertainty for `{era}` (`btagWeight[11/12]`) |
| `eventWeight__tauidjet*up/down` | Tau VSjet ID SF — **decomposed sources only**: 2022–2023: `Stat1/2dm*`, `Systalleras`, `Syst{era}`, `SystTESdm*`; 2024: `Stat1dm{0,1,10,11}2024`. The inclusive `tauidjet`, 0%-variation `Highpt*` and legacy aliases (`Uncert0/1`, `Systdm*`) are no longer booked |
| `eventWeight__tauidelup/down`, `eventWeight__tauidmuup/down` | Tau VSe and VSmu ID SFs |
| `eventWeight__tauFF{stat,syst}up/down` | Data-driven $\text{jet}\to\tau_h$ Fake Factor stat/syst (`--ff` mode) |
| `eventWeight__{mescale,renscale,facscale}up/down` | Inclusive 7-point ME $\mu_R/\mu_F$ scale variations (theory mode; renormalized by `ScaleWeightSum` in `postprocess.py`) |
| `eventWeight__{isr,fsr}up/down` | Parton Shower ISR / FSR scale variations (theory mode; renormalized by `PSWeightSum` in `postprocess.py`) |
| `eventWeight__pdf{1..102}`, `eventWeight__pdfalphasup/down` | PDF eigenvectors and $\alpha_S$ variations (theory mode; renormalized by `LHEPdfWeightSum` in `postprocess.py`) |

**b-tag SF:** Fixed Working Point (Medium) SFs are applied via BTV Method 1a (`btagWeight[0]`), with heavy-flavour (`bc`) and light-flavour (`light`) separated into correlated (`btagSFbc_correlated`, `btagSFlight_correlated`) and era-uncorrelated (`btagSFbc_<era>`, `btagSFlight_<era>`) components following the official CMS BTV POG Wiki naming convention. Per-njet normalization is not required for Method 1a.

### DNN input features

The following 24 kinematic branches (defined in `DNN/utils/feature_config.py`) are used by `DNN/train_multi.py` and `DNN/eval_multi.py`:

```python
# Muon channel (24 features)
[
    "Muon1_pt",
    "Tau1_pt", "Tau1_mass",
    "Jet1_pt", "Jet1_mass",
    "Jet2_pt", "Jet2_mass",
    "Jet3_pt", "Jet3_mass",
    "bJet1_pt", "bJet1_mass",
    "chi2_SMW_mass", "chi2_SMTop_mass", "chi2_SMTop_pt", "chi2_wqq_dR",
    "leptau_mass", "leptau_dR",
    "lepbjet_dR", "taubjet_dR", "lepMET_dPhi", "tauMET_dPhi",
    "PuppiMET_pt", "Jet_HT", "leptaubjet1_mass",
]

# Electron channel: replace "Muon1_pt" with "Electron1_pt"
```


## IX. Postprocessing and Plotting

### b-SF rescaling and histogram merging

```bash
python postprocess.py -I process_v1/muon -Y v12_2023BPix
```

Options:
- `-I` : channel subdirectory inside process output (containing `{year}/`, e.g. `process_v1/muon`)
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
python print_syst_table.py -I process_v1/muon -Y v12_2023BPix -C muon
```

### Tau fake factor workflow

```bash
# Step 1: run in fake-rate measurement mode (output folder must contain 'fake')
python scripts/process.py -V skim_v1 -O process_v1_fake_lss -Y v12_2023BPix -C muon -S nosyst -M lss
python scripts/process.py -V skim_v1 -O process_v1_fake_los -Y v12_2023BPix -C muon -S nosyst -M los
python scripts/process.py -V skim_v1 -O process_v1_fake_tss -Y v12_2023BPix -C muon -S nosyst -M tss
python scripts/process.py -V skim_v1 -O process_v1_fake_tos -Y v12_2023BPix -C muon -S nosyst -M tos

# Step 2: compute fake factors
python fake_factor_calculator_fromROOTFile.py   # per-decay-mode version

# Step 3: apply in final selection (output folder must contain 'FF')
python scripts/process.py -V skim_v1 -O process_v1_FF -Y v12_2023BPix -C muon -S theory --ff
```

### DNN training and inference

```bash
cd ../DNN
# Train
python train_multi.py ...

# Inference / evaluation (`-P` is the base dir above `<ch>/<year>/`, e.g. `process_v1/`)
python eval_multi.py \
    -O DNN_eval \
    -I top_lfv_<timestamp> \
    -C muon \
    -P /path/to/process_v1/ \
    --xsecfile xsec.yaml \
    --alpha 0.1
```


## X. Validation

### 1. Systematic propagation check (`validate_systematics.py`)

Automated 7-category validation suite verifying branch existence, event count invariance, up/down symmetry, non-zero/finite weights, physical bounds, frozen variables, and Z-scores:

```bash
# Muon channel (e.g. 2022EE)
python3 scripts/validate_systematics.py \
    --indir process_v1/muon/v12_2022EE/ \
    --sample TTto2L2Nu \
    --year 2022EE \
    --ch muon \
    --nevents 10000 \
    --outfile report_2022EE_muon.json

# Electron channel (e.g. 2023BPix)
python3 scripts/validate_systematics.py \
    --indir process_v1/electron/v12_2023BPix/ \
    --sample TTto2L2Nu \
    --year 2023BPix \
    --ch electron \
    --nevents 10000 \
    --outfile report_2023BPix_electron.json
```

The script runs 7 categories of automated checks and exits with code 1 if any fail.

### 2. Systematic yield table generation (`print_syst_table.py`)

Batch runner for `plotIt` to calculate yield impact ($\pm \Delta N/N [\%]$) for all systematics and produce `total_syst.tex`:

```bash
# S5 yields with Tau Fake Factor (`-Y` requires v12_.../v15_... because config lookup uses year[4:])
python3 print_syst_table.py -I process_v1/muon -Y v12_2022EE -C muon -F

# DNN discriminant yields (Electron channel)
python3 print_syst_table.py -I DNN_eval/electron -Y v12_2023BPix -C electron -D
```

### 3. External systematic injection (`inject_ext_syst.py`)

Injects dedicated shape variation files (`hdamp`, `tune`, `tes`, `muonscale`, `elescale`, etc.) into nominal files:

```bash
# Dry run preview
python3 scripts/inject_ext_syst.py -d process_v1/muon/v12_2022EE/ --dry-run

# Inject all external systematics
python3 scripts/inject_ext_syst.py -d process_v1/muon/v12_2022EE/
```

*(Note: its default `--systs` list omits `jesEC2*`, `jesHF*`, and the `2022EE`/`2023BPix` sub-era names; pass `--systs ...` explicitly when injecting those.)*

### Output file integrity

```bash
# Check for zombie or incomplete skim files
python scripts/isZoombie.py $LFV_SKIM_DIR/skim_v1/muon/mc/v12_2023BPix --isSkim

# Check for zombie or incomplete process files
python scripts/isZoombie.py process_v1/muon/v12_2023BPix/
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
| `postprocess.py` year not recognised | Use the era tags listed in Section IX (`-Y v12_2022`, `v12_2022EE`, `v12_2023`, `v12_2023BPix`, `v15_2024`) |
| `[ERROR] ... Muon_pt_uncor` / missing `Electron_seedGain`, `Electron_r9`, `Muon_nTrackerLayers` | The skim predates the S&S-before-selection change; re-skim with the current code (Section XII) |
| `[ERROR] syst is data but cannot determine Golden JSON` / `Golden JSON file does not exist` | `--year` (and input path) contain no 2022/2023/2024 tag, or `data/GoldenJSON/Cert_Collisions*.json` is missing |
| `[ERROR] readjson: Problem reading json file` / `Data processing requires a valid Golden JSON` | Golden JSON path passed to the C++ analyzer cannot be opened (or none was given) for Data |
| `[ERROR] get_efficiency: ...` / `BTV Method 1a efficiency map is missing` | `data/BTV/<ch>/<era>/btag_eff*.root` missing or has `eff <= 0` / `eff >= 1` in a bin; (re)derive the maps (Section XII, step 3) |
| `[ERROR] applyWeights: Missing genWeight branch in MC sample!` | Input MC file has no `genWeight`; check the dataset |
| `[ERROR] applyWeights: Failed to load pileup correction ...` | Pileup JSON/key not found for the era (2024 expects `2024_Summer24/puWeights_CDEFGHI.json.gz`, key `Collisions24_CDEFGHI_goldenJSON`) |
| `[ERROR] [tauFF] ...` / `[ERROR] ... unsupported decay mode` | Missing `data/TauFF/tau_fake_factors_run3.json`, unknown channel, or a tau with `dm ∉ {0,1,10,11}` |
| `[ERROR] JES correction key '...' not found` | `-S jes*` name does not map to a key in `data/JME/<era>/jet_jerc.json.gz` |
| `RDataFrame: ... JIT ... no member named 'Cast' in namespace 'ROOT::VecOps'` | Old build: `ROOT::VecOps::Cast` is not available in ROOT 6.36 / LCG_108; rebuild with the current sources (JIT strings now use `* 1` / `* 1.0f`) |

### Log scanning

```bash
# Automated (recommended; input_path must contain log/)
python scripts/checkerror.py process_v1/muon/v12_2023BPix

# Manual
find log/ | xargs grep -l "runtime_error\|fault\|fatal\|Traceback\|ERROR" > ~/errors.txt

# Errors that are benign (known false positives)
# "SimpleJetCorrectionUncertainty" messages in logs are informational, not errors
find log/ | xargs grep ERROR | grep -v SimpleJetCorrectionUncertainty > ~/real_errors.txt
```


## XII. Recommended Execution Order (verified CLI)

> All options below were checked against the `argparse`/`optparse` definitions in the scripts (2026-11-12). Replace `<…>` placeholders. Run every step from `LFVAnalyzer/` after `source /cvmfs/sft.cern.ch/lcg/views/LCG_108/x86_64-el9-gcc15-opt/setup.sh`.

```
make → skim → b-tag efficiency maps → process → postprocess → DNN eval → plotIt / print_syst_table → Combine
```

The b-tag efficiency step **must come after the skim** because `run_all_btag_efficiency.py` reads the skimmed MC. After the 2026 Object Review changes, **all eras and both channels must be re-skimmed** (old skims lack `Muon_pt_uncor`, `Electron_pt_uncor`, `Electron_seedGain`, … and have a different jet/MET/PU content).

**Era tag conventions (`<SKIM_ERA>` vs `<PLAIN_ERA>`):**
- `<SKIM_ERA>` $\in \{\texttt{v12\_2022}, \texttt{v12\_2022EE}, \texttt{v12\_2023}, \texttt{v12\_2023BPix}, \texttt{v15\_2024}\}$ — required by `skim.py`, `postprocess.py`, `print_syst_table.py`, and `DNN/eval_multi.py`.
- `<PLAIN_ERA>` $\in \{\texttt{2022}, \texttt{2022EE}, \texttt{2023}, \texttt{2023BPix}, \texttt{2024}\}$ — expected by `run_all_btag_efficiency.py`, `plotIt` `config_<PLAIN_ERA>.yml`, and `prepareShapesAndCards.py -dataYear`.

| # | Step | Command | Notes |
|---|------|---------|-------|
| 0 | Build | `make clean && make -j 4` | Check that `libnanoadrdframe.so` is created; the build must be warning-free for `impl_Selection.cpp`/`utility.cpp` |
| 1a | Skim — single-file test | `python scripts/skimonefile.py -I <NANOAOD_FILE_OR_XROOTD_URL> -O <OUT>.root -Y <SKIM_ERA> -C muon` | Check `hcounter`, `LHEPdfWeightSum`, `PSWeightSum`, `ScaleWeightSum` in the output; Data skim needs no Golden JSON |
| 1b | Skim — batch | `python scripts/skim.py -V <SKIM_VERSION> -Y <SKIM_ERA> -C muon [--local] [--dry]` | Repeat for all 5 `<SKIM_ERA>` and `-C electron` |
| 1c | Skim — log scan | `python scripts/checkerror.py <SKIM_DIR>/<SKIM_VERSION>/<ch>` | `input_path` is required (reads `<input_path>/log`); every `[ERROR]` in this README's troubleshooting table is caught |
| 2 | b-tag efficiency maps | `python scripts/run_all_btag_efficiency.py -C both -Y 2022 2022EE 2023 2023BPix 2024 -B "<SKIM_DIR>/<SKIM_VERSION>/{channel}/mc/" -j 4` (first with `--dry-run`) | **`-B` is required** (default path is hard-coded to `skim_0812_LFV`); outputs `data/BTV/<ch>/<era>/btag_eff_*.root` |
| 3a | Process — Data | `python scripts/process.py -V <SKIM_VERSION> -O <PROC_OUT> -Y <SKIM_ERA> -C muon -S data -F data` | Golden JSON resolved automatically; `[ERROR]` + exit if not found |
| 3b | Process — MC | `python scripts/process.py -V <SKIM_VERSION> -O <PROC_OUT> -Y <SKIM_ERA> -C muon -S theory -F mc` | `-S nosyst/all/theory`; add `--dry` first. **Check JES naming bug first (`jes_year = year` conflict with `<SKIM_ERA>`, see `Object_Review_TODO.md` §8.3 item 12)** |
| 3c | Process — log scan | `python scripts/checkerror.py <PROC_OUT>/<ch>/<SKIM_ERA>` | |
| 3d | Tau fake factor (optional) | `python scripts/process.py ... -S nosyst -M lss|los|tss|tos` then `--ff` for the final selection | See Section IX |
| 4 | Postprocess | `python postprocess.py -I <PROC_OUT>/<ch> -Y <SKIM_ERA>` | optparse: `-I -Y [--postfix] [-F] [-N]` (no `-V`/`-C`); strictly requires `<SKIM_ERA>` and reads `./<I>/<Y>` |
| 5 | DNN evaluation | `python ../DNN/eval_multi.py -I <TRAIN_DIR> -O <EVAL_DIR> -C muon -P <PROC_OUT> [--xsecfile <YAML>] [--alpha 0.1]` | no `--year`/`--version`; **`-P` must omit `/<ch>`** because `eval_multi.py` appends `<ch>/<SKIM_ERA>/` |
| 6a | Yield / systematic table | `python print_syst_table.py -I <DIR>/<ch> -Y <SKIM_ERA> -C muon [-D] [-F]` | `-D` DNN histograms, `-F` tauFF; uses `year[4:]` for `config_<PLAIN_ERA>.yml`, but also hits `jes_year = year` (§8.3 item 12) |
| 6b | Control plots | `../plotIt/plotIt -o <FIG_DIR> ../plotIt/configs/Run3_muon/config_<PLAIN_ERA>.yml -y -s` (combined: `config_Run3.yml`) | luminosity uncertainty 2024 = 1.6 %, Run 3 = 1.2 % |
| 7 | Combine inputs | `python ../Run3LFVLimits/prepareShapesAndCards.py -p <ROOT_DIR> -C muon -o <CARD_DIR> -xsecfile <files*.yml> -dataYear <PLAIN_ERA>` | single-dash `-xsecfile`/`-dataYear`; `<ROOT_DIR>` is `<EVAL_DIR>/<ch>/<SKIM_ERA>_postprocess*/`; `2023BPix` gets `lumi_1`/`lumi_2`; unknown `dataYear` → `[ERROR]` |
| 7b | Validation | `python scripts/validate_systematics.py --indir <PROC_OUT>/<ch>/<SKIM_ERA>/ --sample <DS> --year <PLAIN_ERA> --ch muon --nevents 10000` | see Section X (note: `--year` sets `jes*_<year>`) |

### Open points to confirm before a full production
- `-Y` naming conflict across `process.py` / `print_syst_table.py` (`jes_year = year`) vs required `<SKIM_ERA>`, and Golden-JSON era detection order — tracked in `recommendations/Object_Review_TODO.md` §8.3 (items 1 and 12).
- The 2026-11-11 build/JIT fixes have not yet been confirmed by a successful `make` + `skimonefile.py` run (§8.1).

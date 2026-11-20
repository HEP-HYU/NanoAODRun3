## LFV Run 3 Top Quark Analysis (`NanoAODRun3`)

### Repository Overview
This repository contains the end-to-end CMS Run 3 (`2022`, `2022EE`, `2023`, `2023BPix`, `2024`, $\sqrt{s} = 13.6\text{ TeV}$) analysis framework for the search for Lepton Flavour Violation (LFV) in top quark interactions ($\mu\tau$ and $e\tau$ channels, AN-25-111):

- **`CommonTools/`**: Shared C++ RDataFrame base classes (`NanoAODAnalyzerrdframe`) and object correction/selection modules:
  - `impl_JetMET.cpp` / `impl_Selection.cpp`: AK4PFPuppi JEC/JER, 11 Regrouped + 4 Pure-Flavour JES & JER 4-vector variations, Type-1 PuppiMET propagation, BTV Fixed-WP Method 1a (`btagSFbc_*`, `btagSFlight_*`), DeepTau2018v2p5 ID SFs.
  - `impl_Leptons.cpp` / `MuonScaReHelper.h` / `ElectronSSHelper.h`: Muon ScaReKIT momentum scale & smearing, Electron EGM scale & smearing, and stat/syst-separated lepton SFs.
  - `impl_Tau.cpp` / `impl_Weights.cpp` / `impl_Reweight.cpp`: Unified Tau Energy Scale (`Tau_ES_nom/up/down`), Pileup reweighting, Top $p_T$ reweighting, and pre-selection LHE/PS weight sums.
- **`LFVAnalyzer/`**: Main analysis package (`SkimEvents`, `TopLFVAnalyzer`, `TauFakeFactorAnalyzer`, `BTagEfficiencyAnalyzer`), batch submission scripts (`scripts/skim.py`, `scripts/process.py`), theory renormalization & merging (`postprocess.py`), yield summary table generator (`print_syst_table.py`), and systematic validator (`scripts/validate_systematics.py`).
- **`DNN/`**: Multi-class and binary DNN training (`train_multi.py`) and systematic shape evaluation (`eval_multi.py`, `eval_multi_bi.py`) on 24 kinematic features.
- **`plotIt/`**: Data/MC control plot and systematic uncertainty calculator with Run 3 era and combined configurations (`configs/Run3_muon/`, `configs/Run3_electron/`).
- **`Run3LFVLimits/`**: CMS Combine / CombineHarvester shape preparation, nuisance correlation mapping, datacard generation (`prepareShapesAndCards.py`), LOWESS smoothing/symmetrization (`symmetrize.py`), shape comparison plots (`plot_syst.py`), and limit/impact/GoF/postfit runners.

### Key Documentation & Systematics References
- Detailed package usage: [`LFVAnalyzer/README.md`](LFVAnalyzer/README.md), [`DNN/README.md`](DNN/README.md), [`Run3LFVLimits/README.md`](Run3LFVLimits/README.md)
- Systematics status & POG checklist: [`../recommendations/Systematics_Review_Status.md`](../recommendations/Systematics_Review_Status.md)
- End-to-end pipeline consistency & bug audit: [`../recommendations/Systematics_Pipeline_Consistency_Report.md`](../recommendations/Systematics_Pipeline_Consistency_Report.md)

### Status of the CMS Object Review round (as of 2026-11-12)
- All six POG/XPOG checklists (EGM, MUO, TAU, JME, BTV, XPOG & Lumi/PU) were audited; **34 of 37 tracked items are completed** and 3 are review notes / pending external input (`EGM-07` combined HLT SF, `MUO-06`, `JME-06`). Item-by-item record: [`../recommendations/Object_Review_TODO.md`](../recommendations/Object_Review_TODO.md).
- Behavioural changes worth knowing when running the code:
  - Muon ScaReKIT / electron S&S are applied **before** object selection at skim level → **all eras and channels must be re-skimmed** and the b-tag efficiency maps re-derived afterwards (`LFVAnalyzer/README.md` Section XII for the verified command order).
  - Data requires a Golden JSON at the **Process** stage (the skim stays unmasked); failures print `[ERROR]` and abort.
  - Missing SFs / branches / correction keys no longer fall back silently: they print `[ERROR]` (detected by `LFVAnalyzer/scripts/checkerror.py`) and throw.
  - Unsplit `muid/muiso/mutrg/elecid/elecreco` and the inclusive `tauidjet` weight variations are no longer booked (only stat/syst and decomposed sources).
- Open points (no code changed): [`../recommendations/Object_Review_TODO.md`](../recommendations/Object_Review_TODO.md) §8.3 and [`../recommendations/Systematics_Review_Status.md`](../recommendations/Systematics_Review_Status.md) §6.3 — notably the `process.py -Y` naming vs JES syst names and the Golden-JSON era-detection order.
- Build status: the 2026-11-11 compile/JIT fixes (`correction::Correction::Ref`, removal of `ROOT::VecOps::Cast` from JIT strings) are not yet confirmed by a successful cluster run.

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

### Recent Updates & Remaining Work (as of 2026-11-25)
- **Completed milestones (`2026-11-23` ~ `11-25`)**:
  - Cluster C++ build/JIT verified; Run 3 skimming completed; b-tag efficiency maps regenerated and uploaded (`02f87ad`); `process.py`/`print_syst_table.py`/`inject_ext_syst.py` JES `-Y` sub-era naming fixed (`e4293ae`); lumi & xsec configs synced (`7431eac`, `f1abfdc`); muon full Run 3 TauFF updated (`4b6f688`) and 2024 `postprocess`/`tauFF` verified.
  - `LFVAnalyzer/plot_era.py` & `plot_run3.py`: `-f` auto-excludes signals and clears `systematics: []` (preserving `luminosity-error`), skips missing low-$p_T$ QCD ($\le 30\text{ GeV}$) with `[WARNING]`, `plot_era.py -f` supports version-prefix auto-expansion + merged data auto-`hadd`, `plot_run3.py` fixes `root: '#'` path handling, and `-F` dynamically extracts all 52 `_S5` plots from `histos.yml` + `histos_reco.yml`.
  - `plotIt/configs/Run3_{muon,electron}/{histos.yml,histos_FFapply.yml}`: Added `show-onlyoverflow: true` on `h_leptau_mass_blind_*` so `-1.0` blinded data in underflow bin 0 is not folded into bin 1.
  - `LFVAnalyzer/print_syst_table.py`: Forces 2-decimal precision (`yields-table-numerical-precision-yields: 2`) in `total_syst.tex`.
  - **TOP-23-008 lepton phase-space extrapolation**: Added `muExtra` (0.5%) and `elExtra` (1.0%) across all 12 Run 3 `plotIt` configs, `print_syst_table.py`, and `Run3LFVLimits/prepareShapesAndCards.py` (`lnN`, correlated across eras).
- **Remaining tasks**:
  - Re-run `plot_era.py -F`, `plot_run3.py -F`, and `print_syst_table.py` (no re-process needed) to pick up the 52 S5 plots, blind-mass bin-1 fix, 2-decimal table, and `muExtra`/`elExtra`.
  - Finalize electron-channel TauFF (if pending), finish all 5 eras × 2 channels through `postprocess`/`plotIt`/`print_syst_table.py`, then run `DNN/eval_multi.py` → `Run3LFVLimits/prepareShapesAndCards.py` → Combine limits/impacts/GoF/postfit.
  - Decide/resolve open analysis items: `Luminosity` row in `total_syst_template.tex`, 100-eigenvector `isPDFenv` in `postprocess.py`, Run 3 TauFF ABCD independence plot (`misid_prob_os_vs_ss` / 2D correlation) & optional $|\eta|$-split / $|\tau_\eta|<2.3$ checks, Golden JSON era-detection order (`processonedataset.py`/`processonefile.py`), `EGM-07`/`MUO-06`/`JME-06` review forms, and syncing/committing the 3 separate repos (`NanoAODRun3`, `plotIt`, `Run3LFVLimits`).

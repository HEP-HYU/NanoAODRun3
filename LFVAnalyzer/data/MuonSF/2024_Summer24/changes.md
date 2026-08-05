## Changes: 2026-06-18 (Fix infinity to inf, update HLT SFs)

Merge Request: [!5](https://gitlab.cern.ch/cms-analysis-corrections/MUO/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/-/merge_requests/5)

Updated the `muon_Z.json.gz` content.

1. Fixed the infinite treatment. (Was using Infinity, but correct one is "inf")
2. Updated the existing muon HLT SFs to latest version, measured with muon's pT vs eta (previous version measured with muon pT vs absolute eta)
3. Added 2 new HLT SFs with medium ID and medium PF isolation as denominator(`NUM_IsoMu24_DEN_CutBasedIdMedium_and_PFIsoTight`, `NUM_IsoMu24_or_Mu50_or_CascadeMu100_or_HighPtTkMu100_DEN_CutBasedIdMedium_and_PFIsoTight`)
4. Updated the `"version": 1` to `"version": 2` for all SFs in `muon_Z.json.gz` file

After the update, verified the file with `correction validate --version 2 <file name>`

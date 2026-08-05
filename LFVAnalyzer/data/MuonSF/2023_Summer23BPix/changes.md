## Changes: 2026-06-18 (Fix infinity to inf, update HLT SFs)

Merge Request: [!2](https://gitlab.cern.ch/cms-analysis-corrections/MUO/Run3-23DSep23-Summer23BPix-NanoAODv12/-/merge_requests/2)

Updated the `muon_Z.json.gz` content.

1. Fixed the infinite treatment. (Was using Infinity, but correct one is "inf")
2. Updated the muon HLT SFs to latest version, measured with muon's pT vs eta (previous version measured with muon pT vs absolute eta)
3. Updated the `"version": 1` to `"version": 2` for all SFs in `muon_Z.json.gz` file

After the update, verified the file with `correction validate --version 2 <file name>`

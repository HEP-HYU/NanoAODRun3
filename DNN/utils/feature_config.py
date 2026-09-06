# feature_config.py
# Single source of truth for DNN input variables.
# Used by: train_multi.py, eval_multi.py, feature_importance.py, hyperparameter_scan.py
#
# Why .py (not .json)?
#   The muon/electron channel branch is a Python expression; keeping it here avoids
#   duplicating the if/else in every script that imports the list.

# Variables shared across both channels
INPUTVARS_COMMON = [
    "Tau1_pt", "Tau1_mass",
    "Jet1_pt", "Jet1_mass",
    "Jet2_pt", "Jet2_mass",
    "Jet3_pt", "Jet3_mass",
    "bJet1_pt", "bJet1_mass",
    "chi2_SMW_mass", "chi2_SMTop_mass", "chi2_SMTop_pt", "chi2_wqq_dR",
    "leptau_mass", "leptau_dR",
    "lepbjet_dR", "taubjet_dR", "lepMET_dPhi", "tauMET_dPhi",
    "PuppiMET_pt", "Jet_HT", "LFV_top_mass",  # m(lep+tau+b)
]

INPUTVARS_MUON     = ["Muon1_pt",    ]
INPUTVARS_ELECTRON = ["Electron1_pt",]


def get_inputvars(ch: str) -> list:
    """Return the ordered input variable list for the given channel ('muon' or 'electron')."""
    if ch == "muon":
        return INPUTVARS_MUON + INPUTVARS_COMMON
    else:
        return INPUTVARS_ELECTRON + INPUTVARS_COMMON


# Signal/background sample lists — also centralised here to avoid per-script drift
SIGLIST_ST = {
    "muon":     ["TCMuTau-LFV-Scalar", "TCMuTau-LFV-Vector", "TCMuTau-LFV-Tensor",
                 "TUMuTau-LFV-Scalar", "TUMuTau-LFV-Vector", "TUMuTau-LFV-Tensor"],
    "electron": ["TCETau-LFV-Scalar",  "TCETau-LFV-Vector",  "TCETau-LFV-Tensor",
                 "TUETau-LFV-Scalar",  "TUETau-LFV-Vector",  "TUETau-LFV-Tensor"],
}

SIGLIST_TT = {
    "muon":     ["TTtoCMuTau-LFV-Scalar", "TTtoCMuTau-LFV-Vector", "TTtoCMuTau-LFV-Tensor",
                 "TTtoUMuTau-LFV-Scalar", "TTtoUMuTau-LFV-Vector", "TTtoUMuTau-LFV-Tensor"],
    "electron": ["TTtoCETau-LFV-Scalar",  "TTtoCETau-LFV-Vector",  "TTtoCETau-LFV-Tensor",
                 "TTtoUETau-LFV-Scalar",  "TTtoUETau-LFV-Vector",  "TTtoUETau-LFV-Tensor"],
}

YEARS = ["v12_2022", "v12_2022EE", "v12_2023", "v12_2023BPix", "v15_2024"]

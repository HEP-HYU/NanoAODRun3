import os
import sys
import pickle
import argparse
import uproot
import pandas as pd
import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt

# Arguments
parser = argparse.ArgumentParser()
parser.add_argument("-C", "--ch", dest="ch", type=str, default="muon")
parser.add_argument("-P", "--project-dir", dest="project_dir", type=str, default="/Users/su/Desktop/antigravity_lfvcode/NanoAODRun3/process_0513_v7/")
parser.add_argument("-M", "--model-dir", dest="model_dir", type=str, required=True, help="Model directory")
parser.add_argument("-O", "--outdir", dest="outdir", type=str, default="physics_validation")
args = parser.parse_args()
ch = args.ch
project_dir = args.project_dir
model_dir = args.model_dir
outdir = args.outdir

os.makedirs(outdir, exist_ok=True)

# Load Model and Scaler
print("Loading model and scaler...")
model = tf.keras.models.load_model(os.path.join(model_dir, "best_model.keras"))
with open(os.path.join(model_dir, "scaler.pkl"), "rb") as f:
    scaler = pickle.load(f)

# Input variables
inputvars = [
    "Tau1_pt","Tau1_mass","Tau1_eta",
    "Jet1_pt","Jet1_mass","Jet1_eta","Jet1_btagPNetB",
    "Jet2_pt","Jet2_mass","Jet2_eta","Jet2_btagPNetB",
    "Jet3_pt","Jet3_mass","Jet3_eta","Jet3_btagPNetB",
    "chi2","chi2_SMW_mass","chi2_SMTop_mass",
    "chi2_wqq_dEta","chi2_wqq_dPhi","chi2_wqq_dR",
    "leptau_mass","leptau_dEta","leptau_dPhi","leptau_dR",
    "PuppiMET_pt"
]
if ch == "muon":
    inputvars = ["Muon1_pt", "Muon1_eta"] + inputvars
else:
    inputvars = ["Electron1_pt", "Electron1_eta"] + inputvars

# Process names
st_signals = {
    "TCMuTau-Scalar": "TCMuTau-LFV-Scalar",
    "TCMuTau-Vector": "TCMuTau-LFV-Vector",
    "TCMuTau-Tensor": "TCMuTau-LFV-Tensor",
    "TUMuTau-Scalar": "TUMuTau-LFV-Scalar",
    "TUMuTau-Vector": "TUMuTau-LFV-Vector",
    "TUMuTau-Tensor": "TUMuTau-LFV-Tensor",
} if ch == "muon" else {
    "TCETau-Scalar": "TCETau-LFV-Scalar",
    "TCETau-Vector": "TCETau-LFV-Vector",
    "TCETau-Tensor": "TCETau-LFV-Tensor",
    "TUETau-Scalar": "TUETau-LFV-Scalar",
    "TUETau-Vector": "TUETau-LFV-Vector",
    "TUETau-Tensor": "TUETau-LFV-Tensor",
}

tt_signals = {
    "TTtoCMuTau-Scalar": "TTtoCMuTau-LFV-Scalar",
    "TTtoCMuTau-Vector": "TTtoCMuTau-LFV-Vector",
    "TTtoCMuTau-Tensor": "TTtoCMuTau-LFV-Tensor",
    "TTtoUMuTau-Scalar": "TTtoUMuTau-LFV-Scalar",
    "TTtoUMuTau-Vector": "TTtoUMuTau-LFV-Vector",
    "TTtoUMuTau-Tensor": "TTtoUMuTau-LFV-Tensor",
} if ch == "muon" else {
    "TTtoCETau-Scalar": "TTtoCETau-LFV-Scalar",
    "TTtoCETau-Vector": "TTtoCETau-LFV-Vector",
    "TTtoCETau-Tensor": "TTtoCETau-LFV-Tensor",
    "TTtoUETau-Scalar": "TTtoUETau-LFV-Scalar",
    "TTtoUETau-Vector": "TTtoUETau-LFV-Vector",
    "TTtoUETau-Tensor": "TTtoUETau-LFV-Tensor",
}

year = "v15_2024"
project_dir_y = os.path.join(project_dir, ch, year)

# Evaluate each ST signal separately
plt.figure(figsize=(10, 6))
for label, rootname in st_signals.items():
    filepath = os.path.join(project_dir_y, f"hist_{rootname}.root")
    tree = uproot.open(filepath)["Events"]
    df = tree.arrays(inputvars, library="pd")
    df = abs(df)
    x = np.array(df.filter(items=inputvars), dtype=np.float32)
    x_scaled = scaler.transform(x)
    pred = model.predict(x_scaled, batch_size=1024, verbose=0)
    
    # Plot probability of ST LFV category
    plt.hist(pred[:, 2], bins=40, histtype="step", density=True, label=label, linewidth=1.5)

plt.xlabel("Probability P(ST LFV)")
plt.ylabel("Normalized Density")
plt.title(f"ST LFV Signal Process Score Dependence ({ch} channel)")
plt.legend(loc="upper left")
plt.yscale("log")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(outdir, "st_signals_dependence.png"))
plt.close()

# Evaluate each TT signal separately
plt.figure(figsize=(10, 6))
for label, rootname in tt_signals.items():
    filepath = os.path.join(project_dir_y, f"hist_{rootname}.root")
    tree = uproot.open(filepath)["Events"]
    df = tree.arrays(inputvars, library="pd")
    df = abs(df)
    x = np.array(df.filter(items=inputvars), dtype=np.float32)
    x_scaled = scaler.transform(x)
    pred = model.predict(x_scaled, batch_size=1024, verbose=0)
    
    # Plot probability of TT LFV category
    plt.hist(pred[:, 1], bins=40, histtype="step", density=True, label=label, linewidth=1.5)

plt.xlabel("Probability P(TT LFV)")
plt.ylabel("Normalized Density")
plt.title(f"TT LFV Signal Process Score Dependence ({ch} channel)")
plt.legend(loc="upper left")
plt.yscale("log")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(outdir, "tt_signals_dependence.png"))
plt.close()

print(f"Physics validation score dependence plots saved to {outdir}")

import os
import sys
import pickle
import argparse
import uproot
import pandas as pd
import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt
from tensorflow.keras.utils import to_categorical

# Arguments
parser = argparse.ArgumentParser()
parser.add_argument("-C", "--ch", dest="ch", type=str, default="muon")
parser.add_argument("-P", "--project-dir", dest="project_dir", type=str, default="/Users/su/Desktop/antigravity_lfvcode/NanoAODRun3/process_0513_v7/")
parser.add_argument("-M", "--model-dir", dest="model_dir", type=str, required=True, help="Directory containing the model and scaler")
parser.add_argument("--seed", dest="seed", type=int, default=42)
args = parser.parse_args()
ch = args.ch
model_dir = args.model_dir
project_dir = args.project_dir
np.random.seed(args.seed)

# Load Model and Scaler
print("Loading model and scaler...")
model = tf.keras.models.load_model(os.path.join(model_dir, "best_model.keras"))
with open(os.path.join(model_dir, "scaler.pkl"), "rb") as f:
    scaler = pickle.load(f)

# Load Data (same logic as train_multi.py)
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

siglist_st = ["TCMuTau-LFV-Scalar", "TCMuTau-LFV-Vector", "TCMuTau-LFV-Tensor", "TUMuTau-LFV-Scalar", "TUMuTau-LFV-Vector", "TUMuTau-LFV-Tensor"] if ch == "muon" else ["TCETau-LFV-Scalar", "TCETau-LFV-Vector", "TCETau-LFV-Tensor", "TUETau-LFV-Scalar", "TUETau-LFV-Vector", "TUETau-LFV-Tensor"]
siglist_tt = ["TTtoCMuTau-LFV-Scalar", "TTtoCMuTau-LFV-Vector", "TTtoCMuTau-LFV-Tensor", "TTtoUMuTau-LFV-Scalar", "TTtoUMuTau-LFV-Vector", "TTtoUMuTau-LFV-Tensor"] if ch == "muon" else ["TTtoCETau-LFV-Scalar", "TTtoCETau-LFV-Vector", "TTtoCETau-LFV-Tensor", "TTtoUETau-LFV-Scalar", "TTtoUETau-LFV-Vector", "TTtoUETau-LFV-Tensor"]

years = ["v15_2024"]
df_sig_st_list = []
df_sig_tt_list = []
df_bkg_tt_list = []

for year in years:
    project_dir_y = os.path.join(project_dir, ch, year)
    
    for sig_tree in siglist_st:
        tree = uproot.open(os.path.join(project_dir_y, f"hist_{sig_tree}.root"))["Events"]
        df = tree.arrays(inputvars, library="pd")
        df_sig_st_list.append(df)
        
    for sig_tree in siglist_tt:
        tree = uproot.open(os.path.join(project_dir_y, f"hist_{sig_tree}.root"))["Events"]
        df = tree.arrays(inputvars, library="pd")
        df_sig_tt_list.append(df)
        
    bkg1 = uproot.open(os.path.join(project_dir_y, "hist_TTto2L2Nu.root"))["Events"].arrays(inputvars, library="pd")
    bkg2 = uproot.open(os.path.join(project_dir_y, "hist_TTtoLNu2Q.root"))["Events"].arrays(inputvars, library="pd")
    bkg1 = bkg1.sample(n=min(len(bkg1), 7*len(bkg2)), random_state=args.seed)
    df_bkg_tt_list.extend([bkg1, bkg2])

df_sig_st = pd.concat(df_sig_st_list)
df_sig_tt = pd.concat(df_sig_tt_list)
df_bkg = pd.concat(df_bkg_tt_list)

nsig = min(len(df_sig_st), len(df_sig_tt), len(df_bkg))
df_sig_st = df_sig_st.sample(n=nsig, random_state=args.seed)
df_sig_tt = df_sig_tt.sample(n=nsig, random_state=args.seed)
df_bkg = df_bkg.sample(n=nsig, random_state=args.seed)

df_sig_st["category"] = 2
df_sig_tt["category"] = 1
df_bkg["category"] = 0

pd_data = pd.concat([df_sig_tt, df_sig_st, df_bkg])
pd_data = abs(pd_data) # Training is on absolute values
pd_data = pd_data.sample(frac=1, random_state=args.seed).reset_index(drop=True).astype('float32')

x_total = np.array(pd_data.filter(items=inputvars))
y_total = np.array(pd_data.filter(items=['category']))
y_cat = to_categorical(y_total)

from sklearn.model_selection import train_test_split
_, x_val, _, y_val = train_test_split(x_total, y_cat, test_size=0.3, random_state=args.seed)

# Evaluate Baseline Performance
x_val_scaled = scaler.transform(x_val)
baseline_loss, baseline_acc = model.evaluate(x_val_scaled, y_val, verbose=0)
print(f"Baseline Validation Loss: {baseline_loss:.4f}, Accuracy: {baseline_acc:.4f}")

# Permutation Feature Importance
importances = []
for i, var_name in enumerate(inputvars):
    x_val_perm = x_val.copy()
    # Shuffle the specific feature values across events
    np.random.shuffle(x_val_perm[:, i])
    x_val_perm_scaled = scaler.transform(x_val_perm)
    perm_loss, perm_acc = model.evaluate(x_val_perm_scaled, y_val, verbose=0)
    
    # Measure impact as drop in accuracy
    acc_drop = baseline_acc - perm_acc
    importances.append((var_name, acc_drop, perm_acc))
    print(f"Shuffling {var_name:20s} | Accuracy: {perm_acc:.4f} (Drop: {acc_drop:.4f})")

# Sort and Save Report
importances.sort(key=lambda x: x[1], reverse=True)
df_imp = pd.DataFrame(importances, columns=["Feature", "AccuracyDrop", "ShuffledAccuracy"])
df_imp.to_csv(os.path.join(model_dir, "feature_importance.csv"), index=False)

# Plot Feature Importance
plt.figure(figsize=(12, 8))
plt.barh(df_imp["Feature"][::-1], df_imp["AccuracyDrop"][::-1], color='steelblue')
plt.xlabel("Accuracy Drop (Importance)")
plt.ylabel("Input Feature")
plt.title(f"Permutation Feature Importance ({ch} channel)")
plt.tight_layout()
plt.savefig(os.path.join(model_dir, "feature_importance.png"))
plt.close()

print(f"Feature importance results saved to {model_dir}")

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
from sklearn.utils import resample

# Arguments
parser = argparse.ArgumentParser()
parser.add_argument("-C", "--ch", dest="ch", type=str, default="muon")
parser.add_argument("-P", "--project-dir", dest="project_dir", type=str, default="/Users/su/Desktop/antigravity_lfvcode/NanoAODRun3/process_0513_v7/")
parser.add_argument("-M", "--model-dir", dest="model_dir", type=str, required=True, help="Model directory")
parser.add_argument("-O", "--outdir", dest="outdir", type=str, default="robustness_test")
parser.add_argument("--seed", dest="seed", type=int, default=42)
args = parser.parse_args()
ch = args.ch
project_dir = args.project_dir
model_dir = args.model_dir
outdir = args.outdir
np.random.seed(args.seed)

os.makedirs(outdir, exist_ok=True)

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
print(f"Baseline Validation Accuracy: {baseline_acc:.4f}")

# 1. Feature Stability Test (+-10% shift)
stability_results = []
for i, var_name in enumerate(inputvars):
    # Shift up
    x_val_shift_up = x_val.copy()
    x_val_shift_up[:, i] *= 1.10
    x_val_scaled_up = scaler.transform(x_val_shift_up)
    _, acc_up = model.evaluate(x_val_scaled_up, y_val, verbose=0)
    
    # Shift down
    x_val_shift_dn = x_val.copy()
    x_val_shift_dn[:, i] *= 0.90
    x_val_scaled_dn = scaler.transform(x_val_shift_dn)
    _, acc_dn = model.evaluate(x_val_scaled_dn, y_val, verbose=0)
    
    max_shift_impact = max(abs(baseline_acc - acc_up), abs(baseline_acc - acc_dn))
    stability_results.append((var_name, acc_up, acc_dn, max_shift_impact))
    print(f"Feature stability {var_name:20s} | Shift +10% Acc: {acc_up:.4f} | Shift -10% Acc: {acc_dn:.4f} | Max Impact: {max_shift_impact:.4f}")

df_stab = pd.DataFrame(stability_results, columns=["Feature", "Accuracy_ShiftUp", "Accuracy_ShiftDown", "MaxImpact"])
df_stab.to_csv(os.path.join(outdir, "feature_stability.csv"), index=False)

# Plot Feature Stability
df_stab_sorted = df_stab.sort_values(by="MaxImpact", ascending=True)
plt.figure(figsize=(12, 8))
plt.barh(df_stab_sorted["Feature"], df_stab_sorted["MaxImpact"], color='orange')
plt.xlabel("Maximum Accuracy Shift on 10% Variation")
plt.ylabel("Input Feature")
plt.title(f"Feature Stability Test - 10% Variation Impact ({ch} channel)")
plt.tight_layout()
plt.savefig(os.path.join(outdir, "feature_stability.png"))
plt.close()

# 2. Bootstrap Validation (100 resamples)
print("Running bootstrap validation...")
n_iterations = 100
bootstrap_accs = []
for i in range(n_iterations):
    x_res, y_res = resample(x_val_scaled, y_val, random_state=args.seed + i)
    _, acc = model.evaluate(x_res, y_res, verbose=0)
    bootstrap_accs.append(acc)

mean_acc = np.mean(bootstrap_accs)
std_acc = np.std(bootstrap_accs)
conf_interval = np.percentile(bootstrap_accs, [2.5, 97.5])

print(f"Bootstrap Validation Accuracy: {mean_acc:.4f} +/- {std_acc:.4f}")
print(f"95% Confidence Interval: [{conf_interval[0]:.4f}, {conf_interval[1]:.4f}]")

with open(os.path.join(outdir, "robustness_report.txt"), "w") as f:
    f.write(f"--- Robustness & Bootstrap Report ({ch} channel) ---\n")
    f.write(f"Mean Bootstrap Accuracy = {mean_acc:.4f}\n")
    f.write(f"Standard Error (SD) = {std_acc:.4f}\n")
    f.write(f"95% Confidence Interval = [{conf_interval[0]:.4f}, {conf_interval[1]:.4f}]\n\n")
    f.write("Feature stability summary table:\n")
    f.write(df_stab.to_string(index=False))

print(f"Robustness results saved to {outdir}")

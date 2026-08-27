# use: python feature_importance.py -C muon -P /path/to/process_dir/ -M top_lfv_date/muon_nom/
# Permutation feature importance against the trained model.

import os
os.environ["CUDA_VISIBLE_DEVICES"] = "3"
import sys
import pickle
import argparse
import uproot
import pandas as pd
import numpy as np
import tensorflow as tf
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from tensorflow.keras.utils import to_categorical
from utils.feature_config import get_inputvars, SIGLIST_ST, SIGLIST_TT, YEARS

# Arguments
parser = argparse.ArgumentParser()
parser.add_argument("-C", "--ch",          dest="ch",          type=str, default="muon")
parser.add_argument("-P", "--project-dir", dest="project_dir", type=str,
                    default="/home/itseyes/github/anti_NanoAODRun3/LFVAnalyzer/process_0715_v9/",
                    help="Processed NanoAOD histogram base directory")
parser.add_argument("-M", "--model-dir",   dest="model_dir",   type=str, required=True,
                    help="Directory containing best_model.keras and scaler.pkl")
parser.add_argument("--seed",              dest="seed",        type=int, default=42)
args = parser.parse_args()

ch          = args.ch
model_dir   = args.model_dir
project_dir = args.project_dir
np.random.seed(args.seed)
tf.random.set_seed(args.seed)

# ── Input variables (from central config) ────────────────────────────────────
inputvars = get_inputvars(ch)
siglist_st = SIGLIST_ST[ch]
siglist_tt = SIGLIST_TT[ch]
years      = YEARS

# ── Load model & scaler ───────────────────────────────────────────────────────
print("Loading model and scaler...")
model_path = os.path.join(model_dir, "best_model.keras")
if not os.path.exists(model_path):
    model_path = os.path.join(model_dir, "best_model.h5")
model = tf.keras.models.load_model(model_path)
with open(os.path.join(model_dir, "scaler.pkl"), "rb") as f:
    scaler = pickle.load(f)

# ── Load data (same pipeline as train_multi.py) ───────────────────────────────
print("Loading data files...")
df_sig_st_list, df_sig_tt_list, df_bkg_tt_list = [], [], []

for eny, year in enumerate(years):
    project_dir_y = os.path.join(project_dir, ch, year)

    for sig_tree in siglist_st:
        tree = uproot.open(os.path.join(project_dir_y, f"hist_{sig_tree}.root"))["Events"]
        df   = tree.arrays(inputvars, library="pd")
        df["year"] = eny
        df_sig_st_list.append(df)

    for sig_tree in siglist_tt:
        tree = uproot.open(os.path.join(project_dir_y, f"hist_{sig_tree}.root"))["Events"]
        df   = tree.arrays(inputvars, library="pd")
        df["year"] = eny
        df_sig_tt_list.append(df)

    bkg1 = uproot.open(os.path.join(project_dir_y, "hist_TTto2L2Nu.root"))["Events"].arrays(inputvars, library="pd")
    bkg2 = uproot.open(os.path.join(project_dir_y, "hist_TTtoLNu2Q.root"))["Events"].arrays(inputvars, library="pd")
    bkg1["year"] = eny
    bkg2["year"] = eny
    bkg1 = bkg1.sample(n=min(len(bkg1), 7 * len(bkg2)), random_state=args.seed)
    df_bkg_tt_list.extend([bkg1, bkg2])

df_sig_st = pd.concat(df_sig_st_list)
df_sig_tt = pd.concat(df_sig_tt_list)
df_bkg    = pd.concat(df_bkg_tt_list)

nsig = min(len(df_sig_st), len(df_sig_tt), len(df_bkg))
df_sig_st = df_sig_st.sample(n=nsig, random_state=args.seed)
df_sig_tt = df_sig_tt.sample(n=nsig, random_state=args.seed)
df_bkg    = df_bkg.sample(n=nsig,    random_state=args.seed)

df_sig_st["category"] = 2
df_sig_tt["category"] = 1
df_bkg["category"]    = 0

pd_data = pd.concat([df_sig_tt, df_sig_st, df_bkg])
# NOTE: abs() removed — same rationale as train_multi.py.
# StandardScaler handles all ranges; abs() destroys eta/dPhi sign information.
pd_data = pd_data.sample(frac=1, random_state=args.seed).reset_index(drop=True).astype("float32")

x_total = np.array(pd_data.filter(items=inputvars))
y_total = np.array(pd_data.filter(items=["category"]))
y_cat   = to_categorical(y_total)

from sklearn.model_selection import train_test_split
_, x_val, _, y_val = train_test_split(x_total, y_cat, test_size=0.3, random_state=args.seed)

# ── Evaluate baseline ─────────────────────────────────────────────────────────
x_val_scaled = scaler.transform(x_val)
baseline_loss, baseline_acc = model.evaluate(x_val_scaled, y_val, verbose=0)
print(f"Baseline Validation Loss: {baseline_loss:.4f}, Accuracy: {baseline_acc:.4f}")

# ── Permutation feature importance ────────────────────────────────────────────
importances = []
for i, var_name in enumerate(inputvars):
    x_val_perm = x_val.copy()
    np.random.shuffle(x_val_perm[:, i])          # shuffle single feature
    x_val_perm_scaled = scaler.transform(x_val_perm)
    perm_loss, perm_acc = model.evaluate(x_val_perm_scaled, y_val, verbose=0)

    acc_drop  = baseline_acc - perm_acc
    loss_rise = perm_loss - baseline_loss
    importances.append((var_name, acc_drop, loss_rise, perm_acc))
    print(f"Shuffling {var_name:25s} | Acc: {perm_acc:.4f} (Drop: {acc_drop:+.4f}) | Loss rise: {loss_rise:+.4f}")

# ── Save & plot ───────────────────────────────────────────────────────────────
importances.sort(key=lambda x: x[1], reverse=True)
df_imp = pd.DataFrame(importances, columns=["Feature", "AccuracyDrop", "LossRise", "ShuffledAccuracy"])
df_imp.to_csv(os.path.join(model_dir, "feature_importance.csv"), index=False)

fig, axes = plt.subplots(1, 2, figsize=(16, max(6, len(inputvars) * 0.35)))

# Accuracy drop
axes[0].barh(df_imp["Feature"][::-1], df_imp["AccuracyDrop"][::-1], color="steelblue")
axes[0].axvline(0, color="black", linewidth=0.8, linestyle="--")
axes[0].set_xlabel("Accuracy Drop (higher = more important)")
axes[0].set_title(f"Permutation Importance — Accuracy\n({ch} channel)")

# Loss rise
axes[1].barh(df_imp["Feature"][::-1], df_imp["LossRise"][::-1], color="salmon")
axes[1].axvline(0, color="black", linewidth=0.8, linestyle="--")
axes[1].set_xlabel("Loss Rise (higher = more important)")
axes[1].set_title(f"Permutation Importance — Loss\n({ch} channel)")

plt.tight_layout()
plt.savefig(os.path.join(model_dir, "feature_importance.png"), dpi=150)
plt.close()

print(f"\nFeature importance results saved to {model_dir}")
print("\nTop 10 most important features:")
print(df_imp[["Feature", "AccuracyDrop", "LossRise"]].head(10).to_string(index=False))
print("\nBottom 5 (candidates for removal if AccuracyDrop ≈ 0):")
print(df_imp[["Feature", "AccuracyDrop", "LossRise"]].tail(5).to_string(index=False))

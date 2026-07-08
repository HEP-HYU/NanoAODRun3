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
from sklearn.metrics import roc_curve, auc, confusion_matrix

# Arguments
parser = argparse.ArgumentParser()
parser.add_argument("-C", "--ch", dest="ch", type=str, default="muon")
parser.add_argument("-P", "--project-dir", dest="project_dir", type=str, default="/Users/su/Desktop/antigravity_lfvcode/NanoAODRun3/process_0513_v7/")
parser.add_argument("-B", "--baseline-dir", dest="baseline_dir", type=str, required=True, help="Baseline model directory")
parser.add_argument("-I", "--improved-dir", dest="improved_dir", type=str, required=True, help="Improved model directory")
parser.add_argument("-O", "--outdir", dest="outdir", type=str, default="model_comparison")
parser.add_argument("--seed", dest="seed", type=int, default=42)
args = parser.parse_args()
ch = args.ch
project_dir = args.project_dir
baseline_dir = args.baseline_dir
improved_dir = args.improved_dir
outdir = args.outdir
np.random.seed(args.seed)

os.makedirs(outdir, exist_ok=True)

# Load Models and Scalers
print("Loading baseline model and scaler...")
model_b = tf.keras.models.load_model(os.path.join(baseline_dir, "best_model.keras"))
with open(os.path.join(baseline_dir, "scaler.pkl"), "rb") as f:
    scaler_b = pickle.load(f)

print("Loading improved model and scaler...")
model_i = tf.keras.models.load_model(os.path.join(improved_dir, "best_model.keras"))
with open(os.path.join(improved_dir, "scaler.pkl"), "rb") as f:
    scaler_i = pickle.load(f)

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
_, x_val, _, y_val_cat = train_test_split(x_total, y_cat, test_size=0.3, random_state=args.seed)
y_val = np.argmax(y_val_cat, axis=1)

# Evaluate and get predictions
print("Evaluating models...")
x_val_b = scaler_b.transform(x_val)
x_val_i = scaler_i.transform(x_val)

pred_b = model_b.predict(x_val_b)
pred_i = model_i.predict(x_val_i)

baseline_loss, baseline_acc = model_b.evaluate(x_val_b, y_val_cat, verbose=0)
improved_loss, improved_acc = model_i.evaluate(x_val_i, y_val_cat, verbose=0)

# 1. Output performance summary text
with open(os.path.join(outdir, "performance_summary.txt"), "w") as f:
    f.write(f"--- Model Comparison Summary ({ch} channel) ---\n")
    f.write(f"Baseline model path: {baseline_dir}\n")
    f.write(f"Improved model path: {improved_dir}\n\n")
    f.write(f"Baseline: Loss = {baseline_loss:.4f}, Accuracy = {baseline_acc:.4f}\n")
    f.write(f"Improved: Loss = {improved_loss:.4f}, Accuracy = {improved_acc:.4f}\n")
    f.write(f"Accuracy Difference: {improved_acc - baseline_acc:+.4f}\n")

print(f"Baseline Accuracy: {baseline_acc:.4f}")
print(f"Improved Accuracy: {improved_acc:.4f} (Difference: {improved_acc - baseline_acc:+.4f})")

# 2. Draw ROC Curves
plt.figure(figsize=(10, 8))
classes = ["Background", "Signal TT LFV", "Signal ST LFV"]
colors_b = ["red", "blue", "green"]
colors_i = ["darkred", "darkblue", "darkgreen"]

for c in range(3):
    # Baseline ROC
    fpr_b, tpr_b, _ = roc_curve((y_val == c).astype(int), pred_b[:, c])
    auc_b = auc(fpr_b, tpr_b)
    plt.plot(fpr_b, tpr_b, linestyle="--", color=colors_b[c], label=f"Baseline {classes[c]} (AUC = {auc_b:.4f})")

    # Improved ROC
    fpr_i, tpr_i, _ = roc_curve((y_val == c).astype(int), pred_i[:, c])
    auc_i = auc(fpr_i, tpr_i)
    plt.plot(fpr_i, tpr_i, linestyle="-", color=colors_i[c], label=f"Improved {classes[c]} (AUC = {auc_i:.4f})")

plt.plot([0, 1], [0, 1], 'k--', alpha=0.5)
plt.xlim(0.0, 1.0)
plt.ylim(0.0, 1.05)
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title(f"ROC Curves comparison ({ch} channel)")
plt.legend(loc="lower right")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(outdir, "roc_comparison.png"))
plt.close()

# 3. Output confusion matrices
def draw_cm(y_true, y_pred_labels, title, filename):
    cm = confusion_matrix(y_true, y_pred_labels)
    cm_norm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(cm_norm, interpolation='nearest', cmap=plt.cm.Blues, vmin=0, vmax=1)
    fig.colorbar(im, ax=ax)
    ax.set(xticks=np.arange(3), yticks=np.arange(3), xticklabels=classes, yticklabels=classes,
           title=title, ylabel='True Label', xlabel='Predicted Label')
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    for i in range(3):
        for j in range(3):
            ax.text(j, i, f"{cm_norm[i, j]:.2f}\n({cm[i, j]})", ha="center", va="center",
                    color="white" if cm_norm[i, j] > 0.5 else "black")
    fig.tight_layout()
    plt.savefig(os.path.join(outdir, filename))
    plt.close()

draw_cm(y_val, np.argmax(pred_b, axis=1), "Baseline Confusion Matrix", "cm_baseline.png")
draw_cm(y_val, np.argmax(pred_i, axis=1), "Improved Confusion Matrix", "cm_improved.png")

# 4. Compare DNN Output Probability Distributions
plt.figure(figsize=(12, 5))
# For ST LFV probability distribution
plt.subplot(1, 2, 1)
plt.hist(pred_b[y_val == 2, 2], bins=40, histtype="step", color="blue", linestyle="--", label="ST LFV (Baseline)")
plt.hist(pred_i[y_val == 2, 2], bins=40, histtype="step", color="blue", linestyle="-", label="ST LFV (Improved)")
plt.hist(pred_b[y_val == 0, 2], bins=40, histtype="step", color="red", linestyle="--", label="Bkg (Baseline)")
plt.hist(pred_i[y_val == 0, 2], bins=40, histtype="step", color="red", linestyle="-", label="Bkg (Improved)")
plt.xlabel("P(ST LFV)")
plt.ylabel("Events")
plt.yscale("log")
plt.title("ST LFV Score Distribution")
plt.legend()

# For TT LFV probability distribution
plt.subplot(1, 2, 2)
plt.hist(pred_b[y_val == 1, 1], bins=40, histtype="step", color="green", linestyle="--", label="TT LFV (Baseline)")
plt.hist(pred_i[y_val == 1, 1], bins=40, histtype="step", color="green", linestyle="-", label="TT LFV (Improved)")
plt.hist(pred_b[y_val == 0, 1], bins=40, histtype="step", color="red", linestyle="--", label="Bkg (Baseline)")
plt.hist(pred_i[y_val == 0, 1], bins=40, histtype="step", color="red", linestyle="-", label="Bkg (Improved)")
plt.xlabel("P(TT LFV)")
plt.ylabel("Events")
plt.yscale("log")
plt.title("TT LFV Score Distribution")
plt.legend()

plt.tight_layout()
plt.savefig(os.path.join(outdir, "score_distributions.png"))
plt.close()

print(f"Comparison results saved in {outdir}")

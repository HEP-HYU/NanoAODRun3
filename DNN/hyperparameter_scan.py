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
from tensorflow.keras.callbacks import EarlyStopping

# Arguments
parser = argparse.ArgumentParser()
parser.add_argument("-C", "--ch", dest="ch", type=str, default="muon")
parser.add_argument("-P", "--project-dir", dest="project_dir", type=str, default="/Users/su/Desktop/antigravity_lfvcode/NanoAODRun3/process_0513_v7/")
parser.add_argument("-O", "--outdir", dest="outdir", type=str, default="scan_results")
parser.add_argument("--seed", dest="seed", type=int, default=42)
args = parser.parse_args()
ch = args.ch
project_dir = args.project_dir
outdir = args.outdir
np.random.seed(args.seed)
tf.random.set_seed(args.seed)

os.makedirs(outdir, exist_ok=True)

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
from sklearn.preprocessing import StandardScaler
x_train, x_val, y_train, y_val = train_test_split(x_total, y_cat, test_size=0.3, random_state=args.seed)

scaler = StandardScaler()
x_train = scaler.fit_transform(x_train)
x_val = scaler.transform(x_val)

# Define configurations to scan
# format: (name, [layer widths], lr, dropout, l2_reg, activation)
configs = [
    ("Baseline (ReLU, default)", [50, 50, 50], 1e-3, 0.0, 1e-4, 'relu'),
    ("ELU default", [50, 50, 50], 1e-3, 0.0, 1e-4, 'elu'),
    ("ELU deeper/wider", [100, 100, 100], 1e-3, 0.0, 1e-4, 'elu'),
    ("ELU with Dropout", [100, 100, 100], 1e-3, 0.2, 1e-4, 'elu'),
    ("ELU high L2 reg", [100, 100, 100], 1e-3, 0.0, 1e-3, 'elu'),
    ("ELU low learning rate", [100, 100, 100], 2e-4, 0.0, 1e-4, 'elu'),
]

results = []

for name, layers, lr, dropout, l2_reg, activation in configs:
    print(f"\nTraining configuration: {name} ...")
    model = tf.keras.models.Sequential()
    model.add(tf.keras.layers.Flatten(input_shape=(x_train.shape[1],)))
    
    for idx, width in enumerate(layers):
        model.add(tf.keras.layers.BatchNormalization())
        if l2_reg > 0:
            model.add(tf.keras.layers.Dense(width, activation=activation, kernel_regularizer=tf.keras.regularizers.l2(l2_reg), kernel_initializer='he_uniform'))
        else:
            model.add(tf.keras.layers.Dense(width, activation=activation, kernel_initializer='he_uniform'))
        if dropout > 0:
            model.add(tf.keras.layers.Dropout(dropout))
            
    model.add(tf.keras.layers.Dense(3, activation='softmax'))
    
    # Optional learning rate schedule
    model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=lr, clipvalue=0.5), loss="categorical_crossentropy", metrics=["accuracy"])
    
    es = EarlyStopping(monitor='val_loss', mode='min', patience=10, verbose=0, restore_best_weights=True)
    
    history = model.fit(x_train, y_train, batch_size=1024, epochs=50, validation_data=(x_val, y_val), callbacks=[es], verbose=0)
    
    best_epoch = np.argmin(history.history['val_loss'])
    val_loss = history.history['val_loss'][best_epoch]
    val_acc = history.history['val_accuracy'][best_epoch]
    train_loss = history.history['loss'][best_epoch]
    train_acc = history.history['accuracy'][best_epoch]
    
    results.append({
        "Configuration": name,
        "Train Loss": train_loss,
        "Train Acc": train_acc,
        "Val Loss": val_loss,
        "Val Acc": val_acc,
        "Best Epoch": best_epoch + 1
    })
    print(f"Result -> Best Epoch: {best_epoch+1} | Val Loss: {val_loss:.4f} | Val Acc: {val_acc:.4f}")

# Make dataframe and print/save
df_results = pd.DataFrame(results)
df_results.to_csv(os.path.join(outdir, "scan_results.csv"), index=False)

print("\nHyperparameter Scan Summary Table:")
print(df_results.to_string(index=False))

# Plot performance comparison
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
bars = ax1.barh(df_results["Configuration"], df_results["Val Acc"], color='skyblue')
ax1.set_xlabel("Validation Accuracy")
ax1.set_title("Validation Accuracy by Configuration")
ax1.set_xlim(0.8, 1.0)
for bar in bars:
    width = bar.get_width()
    ax1.text(width + 0.005, bar.get_y() + bar.get_height()/2, f'{width:.4f}', ha='left', va='center')

bars2 = ax2.barh(df_results["Configuration"], df_results["Val Loss"], color='salmon')
ax2.set_xlabel("Validation Loss")
ax2.set_title("Validation Loss by Configuration (Lower is better)")
for bar in bars2:
    width = bar.get_width()
    ax2.text(width + 0.005, bar.get_y() + bar.get_height()/2, f'{width:.4f}', ha='left', va='center')

plt.tight_layout()
plt.savefig(os.path.join(outdir, "scan_results.png"))
plt.close()

print(f"Hyperparameter scan results and plot saved to {outdir}")

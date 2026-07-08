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
print("Loading data files...")
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

# Define grid
depths = list(range(1, 9))
widths = [64, 128, 256, 512, 1024]
l2_reg = 1e-4
activation = 'elu'

results = []

val_acc_grid = np.zeros((len(depths), len(widths)))
overtrain_gap_grid = np.zeros((len(depths), len(widths)))

# 1. Systematic MLP Grid Search
for d_idx, depth in enumerate(depths):
    for w_idx, width in enumerate(widths):
        print(f"\nTraining MLP: Depth = {depth} layers, Width = {width} nodes ...")
        model = tf.keras.models.Sequential()
        model.add(tf.keras.layers.Flatten(input_shape=(x_train.shape[1],)))
        
        for _ in range(depth):
            model.add(tf.keras.layers.BatchNormalization())
            model.add(tf.keras.layers.Dense(width, activation=activation, kernel_regularizer=tf.keras.regularizers.l2(l2_reg), kernel_initializer='he_uniform'))
            
        model.add(tf.keras.layers.Dense(3, activation='softmax'))
        
        # Compile with CosineDecay scheduler
        batch_size = 1024
        steps_per_epoch = int(np.ceil(len(x_train) / batch_size))
        decay_steps = steps_per_epoch * 30
        lr_schedule = tf.keras.optimizers.schedules.CosineDecay(
            initial_learning_rate=1e-3,
            decay_steps=decay_steps,
            alpha=1e-2
        )
        model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=lr_schedule, clipvalue=0.5), loss="categorical_crossentropy", metrics=["accuracy"])
        
        es = EarlyStopping(monitor='val_loss', mode='min', patience=5, verbose=0, restore_best_weights=True)
        
        # Train for max 30 epochs
        history = model.fit(x_train, y_train, batch_size=batch_size, epochs=30, validation_data=(x_val, y_val), callbacks=[es], verbose=0)
        
        best_epoch = np.argmin(history.history['val_loss'])
        val_loss = history.history['val_loss'][best_epoch]
        val_acc = history.history['val_accuracy'][best_epoch]
        train_loss = history.history['loss'][best_epoch]
        train_acc = history.history['accuracy'][best_epoch]
        loss_gap = val_loss - train_loss
        acc_gap = train_acc - val_acc
        
        results.append({
            "Architecture": "MLP",
            "Depth": depth,
            "Width": width,
            "Train Loss": train_loss,
            "Train Acc": train_acc,
            "Val Loss": val_loss,
            "Val Acc": val_acc,
            "Loss Gap": loss_gap,
            "Acc Gap": acc_gap,
            "Best Epoch": best_epoch + 1
        })
        
        val_acc_grid[d_idx, w_idx] = val_acc
        overtrain_gap_grid[d_idx, w_idx] = loss_gap
        print(f"MLP -> Best Epoch: {best_epoch+1} | Val Loss: {val_loss:.4f} | Val Acc: {val_acc:.4f} | Loss Gap: {loss_gap:.4f}")

# 2. Test Tabular ResNet (Depth=3, Width=256)
print("\nTraining Tabular ResNet: Depth = 3 blocks, Width = 256 nodes ...")
inputs = tf.keras.Input(shape=(x_train.shape[1],))
x = tf.keras.layers.Flatten()(inputs)
# Initial projection layer
x = tf.keras.layers.Dense(256, activation=activation, kernel_initializer='he_uniform')(x)
for _ in range(3):
    y = tf.keras.layers.BatchNormalization()(x)
    y = tf.keras.layers.Dense(256, activation=activation, kernel_regularizer=tf.keras.regularizers.l2(l2_reg), kernel_initializer='he_uniform')(y)
    y = tf.keras.layers.BatchNormalization()(y)
    y = tf.keras.layers.Dense(256, activation=activation, kernel_regularizer=tf.keras.regularizers.l2(l2_reg), kernel_initializer='he_uniform')(y)
    x = tf.keras.layers.Add()([x, y])

outputs = tf.keras.layers.Dense(3, activation='softmax')(x)
model_resnet = tf.keras.Model(inputs, outputs)
model_resnet.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=1e-3, clipvalue=0.5), loss="categorical_crossentropy", metrics=["accuracy"])
es_res = EarlyStopping(monitor='val_loss', mode='min', patience=5, verbose=0, restore_best_weights=True)
hist_resnet = model_resnet.fit(x_train, y_train, batch_size=1024, epochs=30, validation_data=(x_val, y_val), callbacks=[es_res], verbose=0)
best_epoch_res = np.argmin(hist_resnet.history['val_loss'])
results.append({
    "Architecture": "ResNet",
    "Depth": 3,
    "Width": 256,
    "Train Loss": hist_resnet.history['loss'][best_epoch_res],
    "Train Acc": hist_resnet.history['accuracy'][best_epoch_res],
    "Val Loss": hist_resnet.history['val_loss'][best_epoch_res],
    "Val Acc": hist_resnet.history['val_accuracy'][best_epoch_res],
    "Loss Gap": hist_resnet.history['val_loss'][best_epoch_res] - hist_resnet.history['loss'][best_epoch_res],
    "Acc Gap": hist_resnet.history['accuracy'][best_epoch_res] - hist_resnet.history['val_accuracy'][best_epoch_res],
    "Best Epoch": best_epoch_res + 1
})
print(f"ResNet -> Val Loss: {results[-1]['Val Loss']:.4f} | Val Acc: {results[-1]['Val Acc']:.4f} | Loss Gap: {results[-1]['Loss Gap']:.4f}")

# 3. Test Conv1D (Filters=64, Width=256)
print("\nTraining Conv1D Network ...")
inputs_c = tf.keras.Input(shape=(x_train.shape[1],))
# Reshape for convolution: (batch, features, channels) -> (batch, 28, 1)
x_c = tf.keras.layers.Reshape((x_train.shape[1], 1))(inputs_c)
x_c = tf.keras.layers.Conv1D(filters=64, kernel_size=3, activation=activation, kernel_initializer='he_uniform')(x_c)
x_c = tf.keras.layers.MaxPooling1D(pool_size=2)(x_c)
x_c = tf.keras.layers.Flatten()(x_c)
x_c = tf.keras.layers.Dense(256, activation=activation, kernel_initializer='he_uniform')(x_c)
outputs_c = tf.keras.layers.Dense(3, activation='softmax')(x_c)
model_conv = tf.keras.Model(inputs_c, outputs_c)
model_conv.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=1e-3, clipvalue=0.5), loss="categorical_crossentropy", metrics=["accuracy"])
es_conv = EarlyStopping(monitor='val_loss', mode='min', patience=5, verbose=0, restore_best_weights=True)
hist_conv = model_conv.fit(x_train, y_train, batch_size=1024, epochs=30, validation_data=(x_val, y_val), callbacks=[es_conv], verbose=0)
best_epoch_conv = np.argmin(hist_conv.history['val_loss'])
results.append({
    "Architecture": "Conv1D",
    "Depth": 1,
    "Width": 256,
    "Train Loss": hist_conv.history['loss'][best_epoch_conv],
    "Train Acc": hist_conv.history['accuracy'][best_epoch_conv],
    "Val Loss": hist_conv.history['val_loss'][best_epoch_conv],
    "Val Acc": hist_conv.history['val_accuracy'][best_epoch_conv],
    "Loss Gap": hist_conv.history['val_loss'][best_epoch_conv] - hist_conv.history['loss'][best_epoch_conv],
    "Acc Gap": hist_conv.history['accuracy'][best_epoch_conv] - hist_conv.history['val_accuracy'][best_epoch_conv],
    "Best Epoch": best_epoch_conv + 1
})
print(f"Conv1D -> Val Loss: {results[-1]['Val Loss']:.4f} | Val Acc: {results[-1]['Val Acc']:.4f} | Loss Gap: {results[-1]['Loss Gap']:.4f}")

# Save results
df_results = pd.DataFrame(results)
df_results.to_csv(os.path.join(outdir, "scan_results.csv"), index=False)

# Function to plot heatmaps
def plot_heatmap(data, title, filename, cmap, fmt=".4f"):
    fig, ax = plt.subplots(figsize=(9, 8))
    im = ax.imshow(data, cmap=cmap, origin='lower')
    fig.colorbar(im, ax=ax)
    
    ax.set_xticks(np.arange(len(widths)))
    ax.set_yticks(np.arange(len(depths)))
    ax.set_xticklabels(widths)
    ax.set_yticklabels(depths)
    
    ax.set_xlabel("Nodes per Layer (Width)")
    ax.set_ylabel("Number of Layers (Depth)")
    ax.set_title(title)
    
    for i in range(len(depths)):
        for j in range(len(widths)):
            ax.text(j, i, format(data[i, j], fmt), ha="center", va="center",
                    color="white" if data[i, j] > np.median(data) else "black")
            
    fig.tight_layout()
    plt.savefig(os.path.join(outdir, filename))
    plt.close()

# Save MLP Heatmaps
plot_heatmap(val_acc_grid, "Validation Accuracy Heatmap (MLP)", "val_accuracy_heatmap.png", "viridis")
plot_heatmap(overtrain_gap_grid, "Overtraining Gap (Val Loss - Train Loss)", "overtraining_loss_gap_heatmap.png", "plasma")

# Compare Best MLP vs ResNet vs Conv1D
best_mlp = df_results[df_results["Architecture"] == "MLP"].sort_values(by="Val Acc", ascending=False).iloc[0]
comparison_df = pd.DataFrame([
    {"Model": "Best MLP\n(Depth=4, Width=256)", "Accuracy": best_mlp["Val Acc"], "Loss Gap": best_mlp["Loss Gap"]},
    {"Model": "ResNet\n(Depth=3, Width=256)", "Accuracy": results[-2]["Val Acc"], "Loss Gap": results[-2]["Loss Gap"]},
    {"Model": "Conv1D\n(Filters=64)", "Accuracy": results[-1]["Val Acc"], "Loss Gap": results[-1]["Loss Gap"]}
])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
bars1 = ax1.bar(comparison_df["Model"], comparison_df["Accuracy"], color=['skyblue', 'lightgreen', 'salmon'])
ax1.set_ylabel("Validation Accuracy")
ax1.set_title("Validation Accuracy Comparison")
ax1.set_ylim(0.8, 0.9)
for bar in bars1:
    yval = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2.0, yval + 0.002, f"{yval:.4f}", ha='center', va='bottom')

bars2 = ax2.bar(comparison_df["Model"], comparison_df["Loss Gap"], color=['skyblue', 'lightgreen', 'salmon'])
ax2.set_ylabel("Overtraining Loss Gap")
ax2.set_title("Overtraining Gap (Val Loss - Train Loss)\nLower is better")
for bar in bars2:
    yval = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2.0, yval + 0.002 if yval >= 0 else yval - 0.002, f"{yval:.4f}", ha='center', va='bottom')

plt.tight_layout()
plt.savefig(os.path.join(outdir, "architecture_comparison.png"))
plt.close()

print("Scan complete! Output saved to scan_results/")

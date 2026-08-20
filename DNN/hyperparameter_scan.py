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
from utils.feature_config import get_inputvars, SIGLIST_ST, SIGLIST_TT, YEARS

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

# ── Input variables & sample lists (central config) ──────────────────────────
print("Loading data files...")
inputvars  = get_inputvars(ch)
siglist_st = SIGLIST_ST[ch]
siglist_tt = SIGLIST_TT[ch]
years      = YEARS
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
df_bkg = pd.concat(df_bkg_tt_list)

nsig = min(len(df_sig_st), len(df_sig_tt), len(df_bkg))
df_sig_st = df_sig_st.sample(n=nsig, random_state=args.seed)
df_sig_tt = df_sig_tt.sample(n=nsig, random_state=args.seed)
df_bkg = df_bkg.sample(n=nsig, random_state=args.seed)

df_sig_st["category"] = 2
df_sig_tt["category"] = 1
df_bkg["category"] = 0

pd_data = pd.concat([df_sig_tt, df_sig_st, df_bkg])
# NOTE: abs() removed — same rationale as train_multi.py.
# StandardScaler handles all ranges; abs() destroys eta/dPhi sign information.
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

# ── Scan grid ────────────────────────────────────────────────────────────────
# Current production model: Depth=3, Width=256, Dropout=0.3
# Scan covers smaller/larger configurations to find optimal size.
depths        = list(range(1, 6))          # 1-5 layers (keep scan fast)
widths        = [64, 128, 256, 512]        # representative widths
dropout_rates = [0.0, 0.2, 0.3, 0.4]     # 0.0 = no dropout (baseline)
l2_reg        = 1e-4
activation    = 'elu'

results = []

# Grids are keyed by dropout rate (for heatmap at default dropout=0.3)
val_acc_grid      = np.zeros((len(depths), len(widths)))
overtrain_gap_grid = np.zeros((len(depths), len(widths)))
DEFAULT_DROPOUT   = 0.3   # grid is recorded for this dropout value

# ── 1. MLP grid scan (depth × width × dropout) ───────────────────────────────
for d_idx, depth in enumerate(depths):
    for w_idx, width in enumerate(widths):
        for dropout_rate in dropout_rates:
            tag = f"Depth={depth} Width={width} Dropout={dropout_rate}"
            is_current_model = (depth == 3 and width == 256 and dropout_rate == DEFAULT_DROPOUT)
            print(f"\nTraining MLP: {tag}{' [current model]' if is_current_model else ''} ...")

            model = tf.keras.models.Sequential()
            model.add(tf.keras.layers.Flatten(input_shape=(x_train.shape[1],)))

            for _ in range(depth):
                model.add(tf.keras.layers.BatchNormalization())
                model.add(tf.keras.layers.Dense(width, activation=activation,
                                                kernel_regularizer=tf.keras.regularizers.l2(l2_reg),
                                                kernel_initializer='he_uniform'))
                if dropout_rate > 0:
                    model.add(tf.keras.layers.Dropout(dropout_rate))

            model.add(tf.keras.layers.Dense(3, activation='softmax'))

            batch_size = 1024
            steps_per_epoch = int(np.ceil(len(x_train) / batch_size))
            decay_steps = steps_per_epoch * 30
            lr_schedule = tf.keras.optimizers.schedules.CosineDecay(
                initial_learning_rate=1e-3, decay_steps=decay_steps, alpha=1e-2
            )
            model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=lr_schedule, clipvalue=0.5),
                          loss="categorical_crossentropy", metrics=["accuracy"])

            es = EarlyStopping(monitor='val_loss', mode='min', patience=5, verbose=0, restore_best_weights=True)
            history = model.fit(x_train, y_train, batch_size=batch_size, epochs=50,
                                validation_data=(x_val, y_val), callbacks=[es], verbose=0)

            best_epoch = np.argmin(history.history['val_loss'])
            val_loss   = history.history['val_loss'][best_epoch]
            val_acc    = history.history['val_accuracy'][best_epoch]
            train_loss = history.history['loss'][best_epoch]
            train_acc  = history.history['accuracy'][best_epoch]
            loss_gap   = val_loss - train_loss
            acc_gap    = train_acc - val_acc

            results.append({
                "Architecture":   "MLP",
                "Depth":          depth,
                "Width":          width,
                "Dropout":        dropout_rate,
                "CurrentModel":   is_current_model,
                "Train Loss":     train_loss,
                "Train Acc":      train_acc,
                "Val Loss":       val_loss,
                "Val Acc":        val_acc,
                "Loss Gap":       loss_gap,
                "Acc Gap":        acc_gap,
                "Best Epoch":     best_epoch + 1
            })

            # Record in heatmap only for DEFAULT_DROPOUT
            if dropout_rate == DEFAULT_DROPOUT:
                val_acc_grid[d_idx, w_idx]       = val_acc
                overtrain_gap_grid[d_idx, w_idx] = loss_gap

            print(f"  -> Best Epoch: {best_epoch+1:3d} | Val Acc: {val_acc:.4f} | Loss Gap: {loss_gap:+.4f}")

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

# ── Save & summarise results ──────────────────────────────────────────────────
df_results = pd.DataFrame(results)
df_results.to_csv(os.path.join(outdir, "scan_results.csv"), index=False)

# Print current model row for easy reference
current = df_results[df_results["CurrentModel"] == True]
if not current.empty:
    print("\n=== Current model (Depth=3, Width=256, Dropout=0.3) ===")
    print(current[["Depth","Width","Dropout","Val Acc","Val Loss","Loss Gap","Best Epoch"]].to_string(index=False))

# Best overall MLP
best_mlp = df_results[df_results["Architecture"] == "MLP"].sort_values(by="Val Acc", ascending=False).iloc[0]
print(f"\nBest MLP: Depth={int(best_mlp['Depth'])}, Width={int(best_mlp['Width'])}, Dropout={best_mlp['Dropout']:.1f}")
print(f"          Val Acc={best_mlp['Val Acc']:.4f}, Loss Gap={best_mlp['Loss Gap']:+.4f}")

# ── Heatmap for DEFAULT_DROPOUT value ────────────────────────────────────────
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
    ax.set_title(title + f" [Dropout={DEFAULT_DROPOUT}]")
    for i in range(len(depths)):
        for j in range(len(widths)):
            ax.text(j, i, format(data[i, j], fmt), ha="center", va="center",
                    color="white" if data[i, j] > np.median(data) else "black")
    fig.tight_layout()
    plt.savefig(os.path.join(outdir, filename))
    plt.close()

plot_heatmap(val_acc_grid,       "Validation Accuracy Heatmap (MLP)",          "val_accuracy_heatmap.png",       "viridis")
plot_heatmap(overtrain_gap_grid, "Overtraining Gap (Val Loss - Train Loss)",    "overtraining_loss_gap_heatmap.png", "plasma")

# ── Dropout effect plot ───────────────────────────────────────────────────────
# Fix Depth=3, Width=256, vary dropout
mlp_only = df_results[(df_results["Architecture"] == "MLP") &
                       (df_results["Depth"] == 3) &
                       (df_results["Width"] == 256)].sort_values("Dropout")
if not mlp_only.empty:
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(mlp_only["Dropout"], mlp_only["Val Acc"],  marker='o', label="Val Acc")
    ax.plot(mlp_only["Dropout"], mlp_only["Loss Gap"], marker='s', linestyle='--', label="Loss Gap (overtrain)")
    ax.axvline(DEFAULT_DROPOUT, color='red', linestyle=':', label=f"Current ({DEFAULT_DROPOUT})")
    ax.set_xlabel("Dropout Rate")
    ax.set_title("Depth=3, Width=256: Dropout effect")
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "dropout_effect.png"))
    plt.close()

print("\nScan complete! Output saved to", outdir)

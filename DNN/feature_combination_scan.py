# use: python feature_combination_scan.py --mode sfs -C muon -P /home/itseyes/github/NanoAODRun3/LFVAnalyzer/process_0513_v7/ -D sfs_results.json --epochs 1000 --patience 10
import os
import sys
import argparse
import uproot
import json
import pandas as pd
import numpy as np
import tensorflow as tf
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from tensorflow.keras.utils import to_categorical
from tensorflow.keras.callbacks import EarlyStopping
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from concurrent.futures import ThreadPoolExecutor, as_completed

# Arguments
parser = argparse.ArgumentParser()
parser.add_argument("-C", "--ch", dest="ch", type=str, default="muon", help="muon or electron")
parser.add_argument("-P", "--project-dir", dest="project_dir", type=str, default="/home/itseyes/github/NanoAODRun3/LFVAnalyzer/process_0513_v7/")
parser.add_argument("-O", "--outfile", dest="outfile", type=str, default="feature_scan_results.json")
parser.add_argument("--mode", dest="mode", type=str, default="sfs", choices=["sfs", "random"], help="sfs (Sequential Forward Selection) or random (Random Search)")
parser.add_argument("--n-random", dest="n_random", type=int, default=200, help="Number of random combinations to test in random mode")
parser.add_argument("--epochs", dest="epochs", type=int, default=1000, help="Max number of epochs per model run (rely on early stopping)")
parser.add_argument("--patience", dest="patience", type=int, default=10, help="Early stopping patience")
parser.add_argument("--workers", dest="workers", type=int, default=4, help="Number of parallel threads for model evaluation")
parser.add_argument("--seed", dest="seed", type=int, default=42)
args = parser.parse_args()

np.random.seed(args.seed)
tf.random.set_seed(args.seed)

# Limit GPU memory growth to support concurrent threads on GPU
gpus = tf.config.experimental.list_physical_devices('GPU')
if gpus:
    try:
        for gpu in gpus:
            tf.config.experimental.set_memory_growth(gpu, True)
    except RuntimeError as e:
        print(f"GPU growth configuration error: {e}")

# 1. Load Data
print("Loading data files...")
base_vars = [
    "Tau1_pt","Tau1_mass","Tau1_eta",
    "Jet1_pt","Jet1_mass","Jet1_eta","Jet1_btagPNetB",
    "Jet2_pt","Jet2_mass","Jet2_eta","Jet2_btagPNetB",
    "Jet3_pt","Jet3_mass","Jet3_eta","Jet3_btagPNetB",
    "chi2","chi2_SMW_mass","chi2_SMTop_mass",
    "chi2_wqq_dEta","chi2_wqq_dPhi","chi2_wqq_dR",
    "leptau_mass","leptau_dEta","leptau_dPhi","leptau_dR",
    "PuppiMET_pt"
]
if args.ch == "muon":
    all_features = ["Muon1_pt", "Muon1_eta"] + base_vars
else:
    all_features = ["Electron1_pt", "Electron1_eta"] + base_vars

siglist_st = ["TCMuTau-LFV-Scalar", "TCMuTau-LFV-Vector", "TCMuTau-LFV-Tensor", "TUMuTau-LFV-Scalar", "TUMuTau-LFV-Vector", "TUMuTau-LFV-Tensor"] if args.ch == "muon" else ["TCETau-LFV-Scalar", "TCETau-LFV-Vector", "TCETau-LFV-Tensor", "TUETau-LFV-Scalar", "TUETau-LFV-Vector", "TUETau-LFV-Tensor"]
siglist_tt = ["TTtoCMuTau-LFV-Scalar", "TTtoCMuTau-LFV-Vector", "TTtoCMuTau-LFV-Tensor", "TTtoUMuTau-LFV-Scalar", "TTtoUMuTau-LFV-Vector", "TTtoUMuTau-LFV-Tensor"] if args.ch == "muon" else ["TTtoCETau-LFV-Scalar", "TTtoCETau-LFV-Vector", "TTtoCETau-LFV-Tensor", "TTtoUETau-LFV-Scalar", "TTtoUETau-LFV-Vector", "TTtoUETau-LFV-Tensor"]

years = ["v15_2024"]
df_sig_st_list = []
df_sig_tt_list = []
df_bkg_tt_list = []

for year in years:
    project_dir_y = os.path.join(args.project_dir, args.ch, year)

    for sig_tree in siglist_st:
        tree = uproot.open(os.path.join(project_dir_y, f"hist_{sig_tree}.root"))["Events"]
        df = tree.arrays(all_features, library="pd")
        df_sig_st_list.append(df)

    for sig_tree in siglist_tt:
        tree = uproot.open(os.path.join(project_dir_y, f"hist_{sig_tree}.root"))["Events"]
        df = tree.arrays(all_features, library="pd")
        df_sig_tt_list.append(df)

    bkg1 = uproot.open(os.path.join(project_dir_y, "hist_TTto2L2Nu.root"))["Events"].arrays(all_features, library="pd")
    bkg2 = uproot.open(os.path.join(project_dir_y, "hist_TTtoLNu2Q.root"))["Events"].arrays(all_features, library="pd")
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
pd_data = abs(pd_data) # Preprocessing: absolute values
pd_data = pd_data.sample(frac=1, random_state=args.seed).reset_index(drop=True).astype('float32')

print(f"Data loaded. Total events: {len(pd_data)} | Max Features: {len(all_features)}")

# Split and Scale
y_cat = to_categorical(np.array(pd_data['category']))
train_idx, val_idx = train_test_split(np.arange(len(pd_data)), test_size=0.3, random_state=args.seed)

# Define fast model training function
def evaluate_feature_subset(selected_features):
    if len(selected_features) == 0:
        return 0.0, 999.0

    x_sub = np.array(pd_data.filter(items=selected_features))

    x_tr = x_sub[train_idx]
    x_va = x_sub[val_idx]
    y_tr = y_cat[train_idx]
    y_va = y_cat[val_idx]

    # Scale subset
    scaler = StandardScaler()
    x_tr = scaler.fit_transform(x_tr)
    x_va = scaler.transform(x_va)
    model = tf.keras.models.Sequential([
        tf.keras.layers.Flatten(input_shape=(x_tr.shape[1],)),
        tf.keras.layers.BatchNormalization(),
        tf.keras.layers.Dense(256, activation='elu', kernel_regularizer=tf.keras.regularizers.l2(1e-4), kernel_initializer='he_uniform'),
        tf.keras.layers.BatchNormalization(),
        tf.keras.layers.Dense(256, activation='elu', kernel_regularizer=tf.keras.regularizers.l2(1e-4), kernel_initializer='he_uniform'),
        tf.keras.layers.BatchNormalization(),
        tf.keras.layers.Dense(256, activation='elu', kernel_regularizer=tf.keras.regularizers.l2(1e-4), kernel_initializer='he_uniform'),
        tf.keras.layers.Dense(3, activation='softmax')
    ])

    model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=1e-3), loss="categorical_crossentropy", metrics=["accuracy"])
    es = EarlyStopping(monitor='val_loss', mode='min', patience=args.patience, verbose=0, restore_best_weights=True)

    # Train until early stop
    history = model.fit(x_tr, y_tr, batch_size=2048, epochs=args.epochs, validation_data=(x_va, y_va), callbacks=[es], verbose=0)

    best_epoch = np.argmin(history.history['val_loss'])
    val_acc = history.history['val_accuracy'][best_epoch]
    val_loss = history.history['val_loss'][best_epoch]

    return float(val_acc), float(val_loss)

# 2. Search Strategies
results_dict = {}

if args.mode == "sfs":
    print("\nStarting Sequential Forward Selection (SFS)...")
    current_features = []
    best_overall_acc = 0.0

    step = 1
    accuracies = []
    feature_counts = []

    while len(current_features) < len(all_features):
        candidates = [f for f in all_features if f not in current_features]
        step_results = []

        # Parallel evaluation of candidates in the current SFS step
        print(f"\n--- SFS Step {step} | Evaluating {len(candidates)} candidates in parallel (workers={args.workers}) ---")
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = {
                executor.submit(evaluate_feature_subset, current_features + [cand]): cand 
                for cand in candidates
            }

            for future in as_completed(futures):
                cand = futures[future]
                try:
                    acc, loss = future.result()
                    step_results.append((cand, acc, loss))
                    print(f"Tested candidate: {cand} -> Val Acc: {acc:.4f} | Val Loss: {loss:.4f}")
                except Exception as e:
                    print(f"Error evaluating candidate {cand}: {e}")

        # Pick candidate that gives highest validation accuracy
        step_results.sort(key=lambda x: x[1], reverse=True)
        best_cand, best_cand_acc, best_cand_loss = step_results[0]

        if best_cand_acc > best_overall_acc:
            current_features.append(best_cand)
            best_overall_acc = best_cand_acc
            accuracies.append(best_overall_acc)
            feature_counts.append(len(current_features))
            print(f"\n[Step {step}] Added '{best_cand}' | Current Best Acc: {best_overall_acc:.4f}")

            results_dict[f"step_{step}"] = {
                "features": list(current_features),
                "added_feature": best_cand,
                "accuracy": best_cand_acc,
                "loss": best_cand_loss
            }

            # Real-time Autosave at the end of each step
            with open(args.outfile, "w") as f:
                json.dump(results_dict, f, indent=4)
            print(f"Autosaved progress to {args.outfile}")

            step += 1
        else:
            print(f"\nNo single feature addition improved the best validation accuracy ({best_overall_acc:.4f}). Stopping search.")
            break

    # Visualize SFS Curve
    if len(accuracies) > 0:
        plt.figure(figsize=(10, 6))
        plt.plot(feature_counts, accuracies, marker='o', linewidth=2, color='royalblue')
        plt.xticks(feature_counts, [results_dict[f"step_{i}"]["added_feature"] for i in range(1, step)], rotation=45, ha='right')
        plt.xlabel("Features Added")
        plt.ylabel("Validation Accuracy")
        plt.title("Sequential Forward Selection (SFS) Performance Curve")
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.tight_layout()
        plt.savefig("sfs_accuracy_curve.png")
        plt.close()
        print("SFS Curve saved as sfs_accuracy_curve.png")

elif args.mode == "random":
    print(f"\nStarting Random Combination Search ({args.n_random} iterations)...")
    
    accuracies = []
    num_features = []
    top_combinations = []
    
    for i in range(args.n_random):
        k = np.random.randint(5, len(all_features) + 1)
        subset = list(np.random.choice(all_features, size=k, replace=False))
        
        acc, loss = evaluate_feature_subset(subset)
        print(f"Iteration {i+1}/{args.n_random} | Features: {len(subset)} | Val Acc: {acc:.4f} | Loss: {loss:.4f}")
        
        accuracies.append(acc)
        num_features.append(len(subset))
        
        results_dict[f"run_{i+1}"] = {
            "features": subset,
            "accuracy": acc,
            "loss": loss
        }
        top_combinations.append((subset, acc))

    # Save random search plots
    # 1. Scatter plot
    plt.figure(figsize=(9, 6))
    plt.scatter(num_features, accuracies, alpha=0.7, color='indigo')
    plt.xlabel("Number of Features in Combination")
    plt.ylabel("Validation Accuracy")
    plt.title("Random Feature Combination Performance")
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig("random_features_scatter.png")
    plt.close()
    print("Scatter plot saved as random_features_scatter.png")
    
    # 2. Feature frequency in top 10%
    top_combinations.sort(key=lambda x: x[1], reverse=True)
    top_n = max(1, int(0.1 * args.n_random))
    top_comb_sub = top_combinations[:top_n]
    
    feature_counts_dict = {f: 0 for f in all_features}
    for comb, _ in top_comb_sub:
        for f in comb:
            feature_counts_dict[f] += 1
            
    sorted_freq = sorted(feature_counts_dict.items(), key=lambda x: x[1], reverse=True)
    features_sorted, counts = zip(*sorted_freq)
    
    plt.figure(figsize=(12, 8))
    plt.barh(features_sorted[::-1], counts[::-1], color='teal')
    plt.xlabel(f"Occurrences in Top {top_n} Combinations")
    plt.title(f"Feature Selection Frequency in Top 10% Performers (out of {args.n_random} runs)")
    plt.tight_layout()
    plt.savefig("feature_selection_frequency.png")
    plt.close()
    print("Selection frequency chart saved as feature_selection_frequency.png")

# Save results
with open(args.outfile, "w") as f:
    json.dump(results_dict, f, indent=4)

print(f"\nScan complete! Results saved to {args.outfile}")

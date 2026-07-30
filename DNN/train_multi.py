# use: python train_multi.py -C muon -P /home/itseyes/github/NanoAODRun3/LFVAnalyzer/process_0513_v7/ --seed 42
# use: python train_multi.py -C electron
# outdir: top_lfv_date_time

import os
import sys

# os.environ["CUDA_VISIBLE_DEVICES"] = "3"

import uproot
import pandas as pd
import numpy as np
import tensorflow as tf
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
import pickle
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint
from keras.utils import to_categorical
from utils.plots import *
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score, GridSearchCV
from datetime import datetime

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-C", "--ch", dest="ch", type=str, default="muon")
parser.add_argument("-P", "--project-dir", dest="project_dir", type=str, default="/home/itseyes/github/anti_NanoAODRun3/LFVAnalyzer/process_0715_v4/", help="Processed NanoAOD histogram base directory")
parser.add_argument("--seed", dest="seed", type=int, default=42)
args = parser.parse_args()
ch = args.ch
np.random.seed(args.seed)
tf.random.set_seed(args.seed)

def min_max_scaling(series):
    return (series - series.min()) / (series.max() - series.min())


root_dir = os.getcwd().replace("DNN","") # Upper directory
# MODIFY !!!
syst = "nom"
label = "top_lfv"
project_dir = args.project_dir

class_names = ["bkg","sigTT", "sigST"]

print("Start multi LFV Training")
epochs = 1000
inputvars_st = [ # "Muon1_pt","Muon1_eta",
        "Tau1_pt","Tau1_mass","Tau1_eta",
        "Jet1_pt","Jet1_mass","Jet1_eta","Jet1_btagPNetB",
        "Jet2_pt","Jet2_mass","Jet2_eta","Jet2_btagPNetB",
        "Jet3_pt","Jet3_mass","Jet3_eta","Jet3_btagPNetB",
        "chi2","chi2_SMW_mass","chi2_SMTop_mass",
        "chi2_wqq_dEta","chi2_wqq_dPhi","chi2_wqq_dR",
        "leptau_mass","leptau_dEta","leptau_dPhi","leptau_dR",
        "PuppiMET_pt", "PuppiMET_phi"
        ]
if ch == "muon": inputvars_st = ["Muon1_pt", "Muon1_eta"] + inputvars_st
else: inputvars_st = ["Electron1_pt", "Electron1_eta"] + inputvars_st

processed = datetime.now().strftime("%m%d_%H%M")

#"MET_pt" : helps to the expected limits, do not remove from the input vars.
sbratio = 1 # sig:bkg = 1:1

train_outdir = label+"_"+processed+"/"+ch+"_"+syst+"/"
os.makedirs(train_outdir, exist_ok=True)
print ("output dir: ", train_outdir)

siglist_st = []
siglist_tt = []
if ch == "muon":
    siglist_st = ["TCMuTau-LFV-Scalar", "TCMuTau-LFV-Vector", "TCMuTau-LFV-Tensor", "TUMuTau-LFV-Scalar", "TUMuTau-LFV-Vector", "TUMuTau-LFV-Tensor"]
    siglist_tt = ["TTtoCMuTau-LFV-Scalar", "TTtoCMuTau-LFV-Vector", "TTtoCMuTau-LFV-Tensor", "TTtoUMuTau-LFV-Scalar", "TTtoUMuTau-LFV-Vector", "TTtoUMuTau-LFV-Tensor"]
else:
    siglist_st = ["TCETau-LFV-Scalar", "TCETau-LFV-Vector", "TCETau-LFV-Tensor", "TUETau-LFV-Scalar", "TUETau-LFV-Vector", "TUETau-LFV-Tensor"]
    siglist_tt = ["TTtoCETau-LFV-Scalar", "TTtoCETau-LFV-Vector", "TTtoCETau-LFV-Tensor", "TTtoUETau-LFV-Scalar", "TTtoUETau-LFV-Vector", "TTtoUETau-LFV-Tensor"]

years = ["v12_2022", "v12_2022EE", "v12_2023", "v12_2023BPix", "v15_2024"]

df_sig_st_list = []
df_sig_tt_list = []
df_bkg_tt_list = []
#We use all the years together
for eny,year in enumerate(years):
   project_dir_y = project_dir+"/"+ch+"/"+year+"/"
   print("N MC events in eny: " , eny, " year: ", year )

   #Concatinate ST signals
   for sig_tree in siglist_st:
      sig_tree_ = uproot.open(project_dir_y+"hist_"+sig_tree+".root")["Events"]
      df_sig_ = sig_tree_.arrays(inputvars_st,library="pd")
      print("signal ST ",sig_tree," ",  len(df_sig_))
      df_sig_["year"] = eny
      df_sig_st_list.append(df_sig_)

   #Concatinate TT signals
   for sig_tree in siglist_tt:
      sig_tree_ = uproot.open(project_dir_y+"hist_"+sig_tree+".root")["Events"]
      df_sig_ = sig_tree_.arrays(inputvars_st,library="pd")
      print("signal TT ", sig_tree," ", len(df_sig_))
      df_sig_["year"] = eny
      df_sig_tt_list.append(df_sig_)

   bkg1_filedir_tt = project_dir_y+"hist_TTto2L2Nu.root"
   bkg2_filedir_tt = project_dir_y+"hist_TTtoLNu2Q.root"
   bkg1_tree_tt = uproot.open(bkg1_filedir_tt)["Events"]
   bkg2_tree_tt = uproot.open(bkg2_filedir_tt)["Events"]
   df_bkg1_tt = bkg1_tree_tt.arrays(inputvars_st,library="pd")
   df_bkg2_tt = bkg2_tree_tt.arrays(inputvars_st,library="pd")
   df_bkg1_tt["year"] = eny
   df_bkg2_tt["year"] = eny
   n_bkg2 = len(df_bkg2_tt)
   df_bkg1_tt = df_bkg1_tt.sample(n = 7*n_bkg2)
   df_bkg_tt_list.append(df_bkg1_tt)
   df_bkg_tt_list.append(df_bkg2_tt)
   print("background tt ll :", len(df_bkg1_tt ))
   print("background tt lj :", len(df_bkg2_tt ))

df_sig_st = pd.concat(df_sig_st_list)
df_sig_tt = pd.concat(df_sig_tt_list)
df_bkg = pd.concat(df_bkg_tt_list)

ntotsig_tt = len(df_sig_tt)
ntotsig_st = len(df_sig_st)
ntotbkg = len(df_bkg)
print(df_bkg.replace(np.nan, 0))
print(df_sig_st.replace(np.nan, 0))
print(df_sig_tt.replace(np.nan, 0))
print("sig tt:", ntotsig_tt)
print("sig st:",ntotsig_st)
print("tot bkg:",ntotbkg)
nsig_st = ntotsig_st
nsig_tt = ntotsig_tt
nbkg = ntotbkg
nsig = min(nsig_st,nsig_tt)
print("nsig: ", nsig)

if nsig >= nbkg:
    if sbratio == 1: nsig = nbkg
    elif sbratio == 2: nsig = int(nbkg/2)
else:
    if sbratio == 1: nbkg = nsig
    elif nsig>=int(nbkg/2):
       if sbratio == 2: nsig = int(nbkg/2)
    else:
       "Error::Check the number of events!"
       sys.exit()

print("Take LFV : "+str(nsig)+" events")
print("Take TT  : "+str(nbkg)+" events")
df_sig_st = df_sig_st.sample(n=nsig, random_state=args.seed)
df_sig_tt = df_sig_tt.sample(n=nsig, random_state=args.seed)
df_bkg = df_bkg.sample(n=nbkg, random_state=args.seed)
df_sig_st["category"] = 2
df_sig_tt["category"] = 1
df_bkg["category"] = 0

print ("ST LFV: ", len(df_sig_st), " TT LFV: ", len(df_sig_tt), " ttbar: ", len(df_bkg))
pd_data = pd.concat([df_sig_tt,df_sig_st,df_bkg])
pd_data = abs(pd_data)
colnames = pd_data.columns
print(pd_data.head())
print("Col names:",colnames)

#for col in colnames:
#    if col == "category": continue
#    pd_data[col] = min_max_scaling(pd_data[col])

#print(pd_data.head())


pd_sig_st = pd_data[pd_data['category'] == 2]
pd_sig_tt = pd_data[pd_data['category'] == 1]
pd_bkg = pd_data[pd_data['category'] == 0]

print("N MC events after shuffle ====  "  )
for eny, y_name in enumerate(years):
    print("N MC events year %s : ST=%d, TT=%d, BKG=%d" % (y_name, len(pd_sig_st[pd_sig_st['year'] == eny]), len(pd_sig_tt[pd_sig_tt['year'] == eny]), len(pd_bkg[pd_bkg['year'] == eny])))
print("Plotting corr_matrix total")
plot_corrMatrix(pd_data,train_outdir,"total")
print("Plotting corr_matrix ST")
plot_corrMatrix(pd_sig_st,train_outdir,"sig_st")
print("Plotting corr_matrix TT")
plot_corrMatrix(pd_sig_tt,train_outdir,"sig_tt")
print("Plotting corr_matrix bkg")
plot_corrMatrix(pd_bkg,train_outdir,"bkg")
plot_corrMatrix(abs(pd_sig_st.subtract(pd_bkg)),train_outdir,"st-bkg")
pd_data = pd_data.sample(frac=1, random_state=args.seed).reset_index(drop=True)

pd_data = pd_data.astype('float32')
print(pd_data.head().to_string())

x_total = np.array(pd_data.filter(items = inputvars_st))
y_total = np.array(pd_data.filter(items = ['category']))

y_cat = to_categorical(y_total)
print("Y cat: ", y_cat)
# Splitting between training set and cross-validation set
numbertr = len(y_total)
trainlen = int(0.7*numbertr) # Fraction used for training

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
x_train, x_val, y_train, y_val = train_test_split(x_total, y_cat, test_size=0.3, random_state=args.seed)
print(len(x_train),len(x_val),len(y_train),len(y_val))

scaler = StandardScaler()
x_train = scaler.fit_transform(x_train)
x_val = scaler.transform(x_val)

with open(os.path.join(train_outdir, 'scaler.pkl'), 'wb') as f:
    pickle.dump(scaler, f)

patience_epoch = 30
# Early Stopping with Validation Loss for Best Model
es = EarlyStopping(monitor='val_loss', mode='min', verbose=1, patience=patience_epoch)
mc = ModelCheckpoint(train_outdir+'/best_model.h5', monitor='val_loss', mode='min', save_best_only=True)
print("xtrain shape:",x_train.shape)
###################################################
#                      Model                      #
###################################################

model = tf.keras.models.Sequential()

###############    Input Layer      ###############
model.add(tf.keras.layers.Flatten(input_shape = (x_train.shape[1],)))
activation_function='elu'
weight_initializer = 'he_uniform'
l2_factor = 1e-4

###############   Hidden Layer 1    ###############
model.add(tf.keras.layers.BatchNormalization())
model.add(tf.keras.layers.Dense(256, activation=activation_function, kernel_regularizer=tf.keras.regularizers.l2(l2_factor), kernel_initializer=weight_initializer))

###############   Hidden Layer 2    ###############
model.add(tf.keras.layers.BatchNormalization())
model.add(tf.keras.layers.Dense(256, activation=activation_function, kernel_regularizer=tf.keras.regularizers.l2(l2_factor), kernel_initializer=weight_initializer))

###############   Hidden Layer 3    ###############
model.add(tf.keras.layers.BatchNormalization())
model.add(tf.keras.layers.Dense(256, activation=activation_function, kernel_regularizer=tf.keras.regularizers.l2(l2_factor), kernel_initializer=weight_initializer))

###############    Output Layer     ###############
model.add(tf.keras.layers.Dense(3, activation="softmax"))
batch_size = 1024

steps_per_epoch = int(np.ceil(len(x_train) / batch_size))
decay_steps = steps_per_epoch * 150 # Decay over 150 epochs
lr_schedule = tf.keras.optimizers.schedules.CosineDecay(
    initial_learning_rate=1e-3,
    decay_steps=decay_steps,
    alpha=1e-2
)
model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=lr_schedule, clipvalue=0.5), loss="categorical_crossentropy", metrics = ["accuracy"])

model.summary()
hist = model.fit(x_train, y_train, batch_size=batch_size, epochs=epochs,
				validation_data=(x_val,y_val), callbacks=[es, mc])

#tf.keras.models.save_model(model,'model_{epoch:%d}.h5', overwrite=True, include_optimizer=True)
print(hist.history.keys())

loss_train = hist.history['loss']
loss_val = hist.history['val_loss']
plt.plot(loss_train, 'g', label='Training loss')
plt.plot(loss_val, 'b', label='Validation loss')
plt.savefig(train_outdir+"/train_val_loss.png")
plt.close()

accuracy_training = hist.history['accuracy']
accuracy_validation = hist.history['val_accuracy']

plt.plot(accuracy_training, 'g', label='Training loss')
plt.plot(accuracy_validation, 'b', label='Validation loss')
plt.savefig(train_outdir+'/acc_vs_epochs.png')
plt.close()

pred_train = np.argmax(model.predict(x_train), axis=1)
#pred_train = model.predict(x_train)
print("pred_train", pred_train)
print("orig train", y_train)
y_train = np.argmax(y_train, axis=1)
from sklearn.metrics import confusion_matrix
#train_result = pd.DataFrame(np.array([y_train.T[0], pred_train.T[1]]).T, columns=["True", "Pred"])
#pred_val = model.predict(x_val)
pred_val = np.argmax(model.predict(x_val), axis=1)
y_val = np.argmax(y_val, axis=1)
print("pred_val", pred_val)
print("orig val", y_val)
print("conf matrix on train set ")
print(confusion_matrix(y_train, pred_train))
print("conf matrix on val set ")
print(confusion_matrix(y_val, pred_val))
#val_result = pd.DataFrame(np.array([y_val.T[0], pred_val.T[1]]).T, columns=["True", "Pred"])
plot_confusion_matrix(y_val, pred_val, classes=class_names,
        title='Confusion matrix, without normalization', savename=train_outdir+"/confusion_matrix_val.pdf")
plot_confusion_matrix(y_val, pred_val, classes=class_names, normalize=True,
        title='Normalized confusion matrix', savename=train_outdir+"/norm_confusion_matrix_val.pdf")
plot_confusion_matrix(y_train, pred_train, classes=class_names,
        title='Confusion matrix, without normalization', savename=train_outdir+"/confusion_matrix_train.pdf")
plot_confusion_matrix(y_train, pred_train, classes=class_names, normalize=True,
        title='Normalized confusion matrix', savename=train_outdir+"/norm_confusion_matrix_train.pdf")
pred_val = model.predict(x_val)
pred_train = model.predict(x_train)
print(pred_val)
print(pred_val.T)
print(y_val)
print(y_val.T)
val_result = pd.DataFrame(np.array([y_val.T, pred_val.T[1]]).T, columns=["True", "Pred"])
train_result = pd.DataFrame(np.array([y_train.T, pred_train.T[1]]).T, columns=["True", "Pred"])
plot_output_dist(train_result, val_result, sig="tt",savedir=train_outdir)
val_result = pd.DataFrame(np.array([y_val.T, pred_val.T[2]]).T, columns=["True", "Pred"])
train_result = pd.DataFrame(np.array([y_train.T, pred_train.T[2]]).T, columns=["True", "Pred"])
plot_output_dist(train_result, val_result, sig="st",savedir=train_outdir)

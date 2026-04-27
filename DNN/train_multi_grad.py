import os
import sys

#os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"] = "1"

import uproot
import pandas as pd
import numpy as np
import tensorflow as tf
from tensorflow.python.eager import backprop
#tf.enable_eager_execution()
import matplotlib.pyplot as plt
import pickle
#from keras.callbacks import EarlyStopping, ModelCheckpoint
from keras.utils.np_utils import to_categorical
from utils.plots import *
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.model_selection import KFold
from keras.wrappers.scikit_learn import KerasClassifier
from sklearn.model_selection import cross_val_score, GridSearchCV


def min_max_scaling(series):
    return (series - series.min()) / (series.max() - series.min())


root_dir = os.getcwd().replace("DNN","") # Upper directory
# MODIFY !!!
processed = "June2023_GoingtoPrep_normalised_"
syst = "nom"
label = "top_lfv_multiClass"
class_names = ["bkg","sigTT", "sigST"]

print("Start multi LFV Training")
#epochs = 1000
epochs = 100
inputvars_st = [ # "Muon1_pt","Muon1_eta",
        "Tau1_pt","Tau1_mass","Tau1_eta",
        "Jet1_mass","Jet1_eta","Jet1_btagPNetB",
        "Jet2_pt","Jet2_mass","Jet2_eta","Jet2_btagPNetB",
        "Jet3_pt","Jet3_mass","Jet3_eta","Jet3_btagPNetB",
        "chi2","chi2_SMW_mass","chi2_SMTop_mass",
        "chi2_wqq_dEta","chi2_wqq_dPhi","chi2_wqq_dR",
        "leptau_mass","leptau_dEta","leptau_dPhi","leptau_dR",
        "PuppiMET_pt"
        ]
if ch == "muon": inputvars_st = ["Muon1_pt", "Muon1_eta"] + inputvars_st
else: inputvars_st = ["Electron1_pt", "Electron1_eta"] + inputvars_st

#"MET_pt" : helps to the expected limits, do not remove from the input vars.
sbratio = 1 # sig:bkg = 1:1

#kfold = KFold(n_splits=10, shuffle=True)

train_outdir = label+"_"+processed+"/"+syst
os.makedirs(train_outdir, exist_ok=True)

siglist_st = []
siglist_tt = []
if ch == "muon":
    siglist_st = ["TCMuTau-LFV-Scalar", "TCMuTau-LFV-Vector", "TCMuTau-LFV-Tensor", "TUMuTau-LFV-Scalar", "TUMuTau-LFV-Vector", "TUMuTau-LFV-Tensor"]
    siglist_tt = ["TTtoCMuTau-LFV-Scalar", "TTtoCMuTau-LFV-Vector", "TTtoCMuTau-LFV-Tensor", "TTtoUMuTau-LFV-Scalar", "TTtoUMuTau-LFV-Vector", "TTtoUMuTau-LFV-Tensor"]
else:
    siglist_st = ["TCETau-LFV-Scalar", "TCETau-LFV-Vector", "TCETau-LFV-Tensor", "TUETau-LFV-Scalar", "TUETau-LFV-Vector", "TUETau-LFV-Tensor"]
    siglist_tt = ["TTtoCETau-LFV-Scalar", "TTtoCETau-LFV-Vector", "TTtoCETau-LFV-Tensor", "TTtoUETau-LFV-Scalar", "TTtoUETau-LFV-Vector", "TTtoUETau-LFV-Tensor"]
years = ["v2022", "v2022EE", "v2023", "v2023_BPix", "v15_2024"]
years = ["v15_2024"]
project_dir = "/home/itseyes/github/NanoAODRun3/LFVAnalyzer/process_0323_btag_v3/"

df_sig_st_list = []
df_sig_tt_list = []
df_bkg_tt_list = []
#We use all the years together
for year in years:
   project_dir_y = project_dir+"/"+year+"/"

   #Concatinate ST signals
   for sig_tree in siglist_st:
      sig_tree_ = uproot.open(project_dir_y+"hist_"+sig_tree+".root")["Events"]
      df_sig_ = sig_tree_.arrays(inputvars_st,library="pd")
      df_sig_st_list.append(df_sig_)

   #Concatinate TT signals
   for sig_tree in siglist_tt:
      sig_tree_ = uproot.open(project_dir_y+"hist_"+sig_tree+".root")["Events"]
      df_sig_ = sig_tree_.arrays(inputvars_st,library="pd")
      df_sig_tt_list.append(df_sig_)

   bkg1_filedir_tt = project_dir_y+"hist_TTto2L2Nu.root"
   bkg2_filedir_tt = project_dir_y+"hist_TTtoLNu2Q.root"
   bkg1_tree_tt = uproot.open(bkg1_filedir_tt)["Events"]
   bkg2_tree_tt = uproot.open(bkg2_filedir_tt)["Events"]
   df_bkg1_tt = bkg1_tree_tt.arrays(inputvars_st,library="pd")
   df_bkg2_tt = bkg2_tree_tt.arrays(inputvars_st,library="pd")
   df_bkg_tt_list.append(df_bkg1_tt)
   df_bkg_tt_list.append(df_bkg2_tt)

df_sig_st = pd.concat(df_sig_st_list)
df_sig_tt = pd.concat(df_sig_tt_list)
df_bkg = pd.concat(df_bkg_tt_list)

ntotsig_tt = len(df_sig_st)
ntotsig_st = len(df_sig_tt)
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
print(nsig)

if nsig >= nbkg:
    if sbratio == 1:
        nsig = nbkg
    elif sbratio == 2:
        nsig = int(nbkg/2)
else:
    if sbratio == 1:
        nbkg = nsig
    elif nsig>=int(nbkg/2):
        if sbratio == 2:
            nsig = int(nbkg/2)
    else:
        "Error::Check the number of events!"
        sys.exit()

print("Take LFV : "+str(nsig)+" events")
print("Take TT  : "+str(nbkg)+" events")
df_sig_st = df_sig_st.sample(n=nsig)
df_sig_tt = df_sig_tt.sample(n=nsig)
df_bkg = df_bkg.sample(n=nbkg)
df_sig_st["category"] = 2
df_sig_tt["category"] = 1
df_bkg["category"] = 0

pd_data = pd.concat([df_sig_tt,df_sig_st,df_bkg])
pd_data = abs(pd_data)
colnames = pd_data.columns
print(pd_data.head())
print("Col names:",colnames)

for col in colnames:
    if col == "category": continue
    pd_data[col] = min_max_scaling(pd_data[col])

print(pd_data.head())


pd_sig_st = pd_data[pd_data['category'] == 2]
pd_sig_tt = pd_data[pd_data['category'] == 1]
pd_bkg = pd_data[pd_data['category'] == 0]
#for var in inputvars_st: 
#inputvars_st = [ "Muon1_eta",
#       ,"Tau1_mass","Tau1_eta",
#       "Jet1_eta","Jet1_btagDeepFlavB",
#       "Jet2_eta","Jet2_btagDeepFlavB",
#       "Jet3_eta","Jet3_btagDeepFlavB",
#       "chi2"
#       "chi2_wqq_dEta","chi2_wqq_dPhi","chi2_wqq_dR",
#       "mutau_dEta","mutau_dPhi","mutau_dR",
#for var in ["chi2_SMW_mass","chi2_SMTop_mass","Jet1_pt","Jet2_pt","Muon1_pt","mutau_mass","Tau1_pt"]: #3000 
#for var in ["Jet3_mass","Jet2_mass","Jet1_mass","Jet3_pt"]: #300 
#for var in ["Tau1_mass"]: #5 
#for var in ["Muon1_eta","Tau1_eta","Jet1_eta","Jet2_eta","Jet3_eta","chi2_wqq_dEta","chi2_wqq_dPhi","mutau_dEta","mutau_dPhi"]: 
#for var in ["chi2"]: 
## Plot signal and background overlayed plots
#plt.hist(pd_bkg[var], bins=30,range=[0, 10000] ,alpha=0.5, label='Background', color='red')
#plt.hist(pd_sig_st[var], bins=30,range=[0, 10000] , alpha=0.1, label='Signal', color='blue', ec="blue")
##plt.hist(abs(pd_bkg[var]), bins=30,alpha=0.5, label='Background', color='red')
##plt.hist(abs(pd_sig_st[var]), bins=30, alpha=0.1, label='Signal', color='blue', ec="blue")
#plt.xlabel(var)
#plt.ylabel('Count')
#plt.title('Signal vs. Background Overlayed Plots')
#plt.legend()
#plt.savefig('plots/signal_background_overlayed_plots_'+var+'.png')
#plt.close()



#print("Plotting corr_matrix total")
#plot_corrMatrix(pd_data,train_outdir,"total")
#print("Plotting corr_matrix ST")
#plot_corrMatrix(pd_sig_st,train_outdir,"sig_st")
#print("Plotting corr_matrix TT")
#plot_corrMatrix(pd_sig_tt,train_outdir,"sig_tt")
#print("Plotting corr_matrix bkg")
#plot_corrMatrix(pd_bkg,train_outdir,"bkg")

pd_data = pd_data.sample(frac=1).reset_index(drop=True)

print(pd_data.head())

#    # Min-Max Scale
#    scaler = MinMaxScaler()
#    pd_data = pd.DataFrame(scaler.fit_transform(pd_data), columns=pd_data.columns)

x_total = np.array(pd_data.filter(items = inputvars_st))
y_total = np.array(pd_data.filter(items = ['category']))

y_cat = to_categorical(y_total)
print("Y cat: ", y_cat)
# Splitting between training set and cross-validation set
numbertr = len(y_total)
trainlen = int(0.7*numbertr) # Fraction used for training

from sklearn.model_selection import train_test_split
x_train, x_val, y_train, y_val = train_test_split(x_total, y_cat, test_size=0.3)
print(len(x_train),len(x_val),len(y_train),len(y_val))
'''
x_train = x_total[:trainlen,0::]
y_train = y_total[:trainlen]
print("y_train:" , y_train)
x_val = x_total[trainlen:,0::]
y_val = y_total[trainlen:]
'''

#patience_epoch = 30
#patience_epoch = 5
# Early Stopping with Validation Loss for Best Model
#es = EarlyStopping(monitor='val_loss', mode='min', verbose=1, patience=patience_epoch)
#mc = ModelCheckpoint(train_outdir+'/best_model.h5', monitor='val_loss', mode='min', save_best_only=True)
print("xtrain shape:",x_train.shape)

model_dir = train_outdir+'/best_model.h5'
print("model dir", model_dir)
model = tf.keras.models.load_model(model_dir)
model.summary()
#pred_train = model.predict_classes(x_train)
#print("pred_train", pred_train)
#print("orig train", y_train)
#y_train = np.argmax(y_train, axis=1)
#from sklearn.metrics import confusion_matrix
#train_result = pd.DataFrame(np.array([y_train.T[0], pred_train.T[1]]).T, columns=["True", "Pred"])
#pred_val = model.predict_classes(x_val)
#y_val = np.argmax(y_val, axis=1)
#print("pred_val", pred_val)
#print("orig train", y_val)

###Feature importance###

print("feature importance starts")
input_data = tf.convert_to_tensor(x_val, dtype=tf.float32)
print(input_data)
name_inputvar = inputvars_st
outputDir = train_outdir
n_evts = len(x_val)
n_var = len(name_inputvar)
mean_grads = n_var*[0.0]
all_grads = []
#mean_jacobian = np.zeros(len(name_inputvar))
#jacobian_matrix = np.zeros((len(name_inputvar),len(name_inputvar)))
for i, event in enumerate(x_val):
  print(i,"/",n_evts)
  with tf.GradientTape() as tape:
    inputs = tf.Variable([event])
    tape.watch(inputs)
    prediction = model(inputs)[:, 1]
  first_order_gradients = tape.gradient(prediction, inputs)
  gradiants = tf.abs(first_order_gradients)
  numpy_array = gradiants.numpy()[0]
  all_grads.append(numpy_array)
  #print(i,numpy_array , len(numpy_array))
  for n in range(len(mean_grads)):
    mean_grads[n] += abs(numpy_array[n])/n_evts

print(mean_grads)
df = pd.DataFrame(all_grads)
print(df)
df.to_csv('data.csv', index=False)

#order_gradiants = tape2.gradient(first_order_gradients, input_data)
#print(first_order_gradients)
#print(second_order_gradiants)
#jacobian = tape2.jacobian(gradients_2, input_data)
#print(jacobian)
#for i in range(n_var):
#	mean_grads[i] += abs(gradients_2.numpy()[0][i])/n_evts
#	mean_jacobian[i] += abs(jacobian[i][0][i])/n_evts
#	for j in range(len(name_inputvar)): jacobian_matrix[i][j] += abs(jacobian[i][0][j])/n_evts

'''

gradiants = tf.abs(first_order_gradients)
numpy_array = gradiants.numpy()
df = pd.DataFrame(numpy_array)
print(df)
df.to_csv('data.csv', index=False)
feature_importance_first_order  = tf.reduce_mean(tf.abs(first_order_gradients), axis=0)
feature_importance_dict = dict(zip(name_inputvar, feature_importance_first_order.numpy()))
#feature_importance_Secondorder  = tf.reduce_mean(tf.abs(second_order_gradients), axis=0)
feature_importance_series = pd.Series(feature_importance_dict)

for i, importance_score in enumerate(feature_importance_first_order):
   print(f"Feature {i+1} , {name_inputvar[i]} Importance: {importance_score:.4f}")


print(feature_importance_series.index, feature_importance_series.values)

#plt.barh(feature_importance_series.index, feature_importance_series.values)
#plt.xlabel('Feature Importance')
#plt.ylabel('Feature Names')
#plt.savefig(train_outdir+'plots/first_order_gradient_importance.png')
#plt.title('First-Order Gradient Feature Importance')
#plt.show()
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
'''

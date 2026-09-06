# use: python eval_multi.py -O DNN_ -I top_lfv_date_time -C muon -P /home/itseyes/github/NanoAODRun3/LFVAnalyzer/process_0513_v7/ --xsecfile xsec.yaml --alpha 0.1
# use: python eval_multi.py -I top_lfv_date_time -C muon
# outdir: DNN_datetime

import os
import sys
os.environ["CUDA_VISIBLE_DEVICES"] = "2"
import uproot
import pandas as pd
import awkward as ak
import numpy as np
import argparse
import tensorflow as tf
import matplotlib.pyplot as plt
import tensorflow.keras.backend as K
from keras.callbacks import EarlyStopping, ModelCheckpoint
import multiprocessing
from datetime import datetime
import pickle
import yaml
from utils.feature_config import get_inputvars

gpus = tf.config.list_physical_devices('GPU')
if gpus:
    try:
        tf.config.set_logical_device_configuration(gpus[0], [tf.config.LogicalDeviceConfiguration(memory_limit=1024)])
    except RuntimeError as e:
        print(e)

default_proj_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "process_0715_v4")
parser = argparse.ArgumentParser(usage="%prog [options]")
parser.add_argument("-O", "--outdir", dest="outdir", type=str, default="DNN_", help="Evaluation folder in your working directory")
parser.add_argument("-I", "--indir", dest="indir", type=str, default="top_lfv_0715_v4", help="Training folder in your working directory")
parser.add_argument("-C", "--ch", dest="ch", type=str, default="muon")
parser.add_argument("-P", "--project-dir", dest="project_dir", type=str, default="/home/itseyes/github/anti_NanoAODRun3/LFVAnalyzer/process_0715_v4/", help="Processed NanoAOD histogram base directory")
parser.add_argument("--xsecfile", dest="xsecfile", type=str, default="", help="Optional YAML file. If set, evaluate only ROOT files listed in the YAML keys")
parser.add_argument("--alpha", dest="alpha", type=float, default=0.1, help="DNN score mixing: ((1-alpha)*p_ST + alpha*p_TT)/p_bkg")
options = parser.parse_args()
ch = options.ch
outdir = options.outdir+datetime.now().strftime("%m%d")

def min_max_scaling(series):
    return (series - series.min()) / (series.max() - series.min())


training_path = options.indir

def load_xsec_rootfiles(xsecfile):
    if not xsecfile:
        return None
    if yaml is None:
        raise ImportError("PyYAML is required to use --xsecfile")
    with open(xsecfile, "r") as handle:
        data = yaml.safe_load(handle)
    if not data:
        raise RuntimeError("Could not load any dataset from %s" % xsecfile)
    return set([name for name in data.keys() if name.endswith(".root")])

# Input variables loaded from central config (utils/feature_config.py)
inputvars = get_inputvars(ch)


discriminators = {"p_st" : 2, "p_tt" : 1 , "p_bkg" : 0 , "p_st_tt" : 999, "p_st_tt_ob" : 999 }

def run(inputs):
    year = inputs[0]
    input_file = inputs[1]
    discriminator_key = inputs[2]
    alpha = inputs[3]
    print(year, discriminator_key , alpha)

    binedges = [0,1,2,5,10,30,100]
    binedges2 = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    if "ob" in discriminator_key : binedges = [0,0.01,1,2,5,10,30,100]
    #binedges = [i for i in frange(0.0, 100.01, 0.01)]

    model_dir = os.path.join(training_path, ch+"_nom/best_model.h5")
    if not os.path.exists(model_dir):
        model_dir = os.path.join(training_path, ch+"_nom/best_model.keras")
    print ("model: ", model_dir)
    model = tf.keras.models.load_model(model_dir)

    scaler_dir = os.path.join(training_path, ch+"_nom/scaler.pkl")
    print ("scaler: ", scaler_dir)
    with open(scaler_dir, "rb") as f:
        scaler = pickle.load(f)

    hists_path = os.path.join(outdir+"/"+ch, year)
    if not os.path.isdir(hists_path):
        os.makedirs(hists_path, exist_ok=True)

    eventWeight = "eventWeight"
    weights = [eventWeight]

    outf_dir = os.path.join(hists_path, input_file.split("/")[-1:][0])

    #print("IN PATH : ", input_file) 
    print("OUT PATH : ", outf_dir)

    infile = uproot.open(input_file)
    tree = infile["Events"]
    branch_names = tree.keys()
    syst_list = [syst.split("__")[1] for syst in branch_names if "eventWeight__" in syst]
    hcounter = infile["hcounter"]
    nEvents = len(tree)
    print("1 - N EVENTS : ", nEvents)
    syst_hist_dnn_dict = {}
    syst_hist_dnn_entries_dict = {}
    hist_nevents_S4_dict = {}
    hist_nevents_S5_dict = {}
    syst_extend = []

    if (not "Muon" in input_file) and (not "egamma" in input_file.lower()):
        infile_forS = uproot.open(input_file.replace("_FF", ""))
        h_nevents_S4_nobtag = infile_forS["h_nevents_S4_nobtag"]  ### get it from removed FF 
        h_nevents_S4 = infile_forS["h_nevents_S4"]
        h_nevents_S2_nobtag = infile_forS["h_nevents_S2_nobtag"]
        h_nevents_S2 = infile_forS["h_nevents_S2"]
        #if any(string in input_file for string in ["Tau-LFV", "TTt", "_ST_t"]) and "__" not in input_file:
        #    ScaleWeightSum = infile['ScaleWeightSum']
        #    PSWeightSum = infile['PSWeightSum']
        #    LHEPdfWeightSum = infile['LHEPdfWeightSum']

    if nEvents == 0:
        print("No events : "+input_file)
        #Need to add empth histograms for technical reasons ....
        muon_pt = []
        tau_pt = []
        pred = []
        pred_st = []
        pred_tt = []
        pred_bkg = []
        pd_weight = []
        dnnhist_nom = np.histogram(pred, bins=binedges, weights=pd_weight, density=False)
        dnnhist_entries_nom = np.histogram(pred, bins=binedges, density=False)
        dnnhist_st = np.histogram(pred_st, bins=binedges2, weights=pd_weight, density=False)
        dnnhist_entries_st = np.histogram(pred_st, bins=binedges2, density=False)
        dnnhist_tt = np.histogram(pred_tt, bins=binedges2, weights=pd_weight, density=False)
        dnnhist_entries_tt = np.histogram(pred_tt, bins=binedges2, density=False)
        dnnhist_bkg = np.histogram(pred_bkg, bins=binedges2, weights=pd_weight, density=False)
        dnnhist_entries_bkg = np.histogram(pred_bkg, bins=binedges2, density=False)
        print ("syst_list: ", syst_list)
        for syst in syst_list:
            syst_hist_dnn_dict["h_dnn_pred_S5__"+syst] = dnnhist_nom
            syst_hist_dnn_entries_dict["h_dnn_entries_S5__"+syst] = dnnhist_entries_nom
            if 'btag' in syst:
                hist_nevents_S4_dict["h_nevents_S4__"+syst] = infile_forS["h_nevents_S4__"+syst]
                hist_nevents_S5_dict["h_nevents_S5__"+syst] = infile["h_nevents_S5__"+syst]
    else:
        muon_pt = None
        if ch == "muon": muon_pt = tree["Muon1_pt"].array()
        else: muon_pt = tree["Electron1_pt"].array()
        tau_pt = tree["Tau1_pt"].array()
        pd_data = tree.arrays(inputvars, library="pd")
        # NOTE: abs() removed — StandardScaler handles all ranges;
        # abs() would destroy eta/dEta/dPhi/phi sign information.
        pd_weight = tree.arrays(weights, library="np")
        pred_data = np.array(pd_data.filter(items = inputvars), dtype=np.float32)

        print(type(pred_data))

        try:
            print("len =", len(pred_data))
        except Exception as e:
            print("len failed:", e)

        if isinstance(pred_data, np.ndarray):
            print(pred_data.shape)
        if len(pred_data) == 0:
            print("No events : "+input_file)
            muon_pt = []
            tau_pt = []
            pred = []
            pred_st = []
            pred_tt = []
            pred_bkg = []
            pd_weight = []
            dnnhist_nom = np.histogram(pred, bins=binedges, weights=pd_weight, density=False)
            dnnhist_entries_nom = np.histogram(pred, bins=binedges, density=False)
            dnnhist_st = np.histogram(pred_st, bins=binedges2, weights=pd_weight, density=False)
            dnnhist_entries_st = np.histogram(pred_st, bins=binedges2, density=False)
            dnnhist_tt = np.histogram(pred_tt, bins=binedges2, weights=pd_weight, density=False)
            dnnhist_entries_tt = np.histogram(pred_tt, bins=binedges2, density=False)
            dnnhist_bkg = np.histogram(pred_bkg, bins=binedges2, weights=pd_weight, density=False)
            dnnhist_entries_bkg = np.histogram(pred_bkg, bins=binedges2, density=False)
            print ("syst_list: ", syst_list)
            for syst in syst_list:
                syst_hist_dnn_dict["h_dnn_pred_S5__"+syst] = dnnhist_nom
                syst_hist_dnn_entries_dict["h_dnn_entries_S5__"+syst] = dnnhist_entries_nom
                if 'btag' in syst:
                    hist_nevents_S4_dict["h_nevents_S4__"+syst] = infile_forS["h_nevents_S4__"+syst]
                    hist_nevents_S5_dict["h_nevents_S5__"+syst] = infile["h_nevents_S5__"+syst]

        else:
            pred_data = scaler.transform(pred_data)
            pred = model.predict(pred_data, batch_size=128, verbose=0)
            #pred_df = pd.DataFrame(pred, columns=['Prediction1', 'Prediction2', 'Prediction3'])
            #print("3 - PRED SHAPE : ", pred.shape)
            #result_df = pd.concat([pd_data, pred_df], axis=1)
            #result_df.to_csv(outf_dir.split(".root")[0]+'combined_data.csv', index=False)
            #print("PRED COMB SHAPE : ", result_df.shape)
            #print("DF saved here: ", outf_dir.split(".root")[0]+'combined_data.csv')

            if discriminator_key == "p_st_tt_ob"  :
                rmzeros = pred[:,0]
                rmzeros[rmzeros<=0.00] = 0.001
                pred[:,0] = rmzeros
                pred_st = pred[:,2].tolist()
                pred_tt = pred[:,1].tolist()
                pred_bkg = pred[:,0].tolist()
                if alpha <10: pred = ( ( (1-alpha) * pred[:,2] + alpha * pred[:,1] ) / pred[:,0] ).tolist()
                else: pred = ( ( pred[:,2] + pred[:,1] ) / pred[:,0] ).tolist()
                pred = np.array(pred)
                pred[pred>=100.0] = 99.999
                pred = pred.tolist()
            elif discriminator_key == "p_st_tt"  : pred = ( pred[:,2] + pred[:,1] ).tolist()
            elif discriminator_key == "p_max"  : pred = ( np.amax(pred, axis=1) ).tolist()
            elif discriminator_key == "p_mean"  : pred = ( np.mean(pred, axis=1) ).tolist()
            elif discriminator_key == "p_min"  : pred = ( np.min(pred, axis=1) ).tolist()
            else  : pred = ( pred[:, discriminators[discriminator_key]] ).tolist()

            #pred = ( ((1 - alpha) * pred[:,2] + alpha * pred[:,1]) / pred[:,0])
            #pred[pred >= 100.0] = 99.999
            #pred = pred.tolist()
            #print ("pred: ", pred)
            pd_weight = pd_weight[eventWeight].tolist()
            nom_weight = pd_weight.copy()
            dnnhist_nom = np.histogram(pred, bins=binedges, weights=pd_weight, density=False)
            dnnhist_entries_nom = np.histogram(pred, bins=binedges, density=False)
            dnnhist_st = np.histogram(pred_st, bins=binedges2, weights=pd_weight, density=False)
            dnnhist_entries_st = np.histogram(pred_st, bins=binedges2, density=False)
            dnnhist_tt = np.histogram(pred_tt, bins=binedges2, weights=pd_weight, density=False)
            dnnhist_entries_tt = np.histogram(pred_tt, bins=binedges2, density=False)
            dnnhist_bkg = np.histogram(pred_bkg, bins=binedges2, weights=pd_weight, density=False)
            dnnhist_entries_bkg = np.histogram(pred_bkg, bins=binedges2, density=False)

            print ("syst_list: ", syst_list)
            for syst in syst_list:
                pd_weight = tree.arrays(["eventWeight__"+syst], library="np")
                pd_weight = pd_weight["eventWeight__"+syst].tolist()

                # set SF for mu / tau 100 GeV above ro below
                if any(s_ in syst for s_ in ["mescale", "renscale", "facscale"]):
                    tmp_weight_mu1ta1 = []
                    tmp_weight_mu1ta2 = []
                    tmp_weight_mu2ta1 = []
                    tmp_weight_mu2ta2 = []
                    for ei in range(len(pd_weight)):
                        if muon_pt[ei] < 100 and tau_pt[ei] < 100:
                            tmp_weight_mu1ta1.append(pd_weight[ei])
                            tmp_weight_mu1ta2.append(nom_weight[ei])
                            tmp_weight_mu2ta1.append(nom_weight[ei])
                            tmp_weight_mu2ta2.append(nom_weight[ei])
                        elif muon_pt[ei] < 100 and tau_pt[ei] >= 100:
                            tmp_weight_mu1ta1.append(nom_weight[ei])
                            tmp_weight_mu1ta2.append(pd_weight[ei])
                            tmp_weight_mu2ta1.append(nom_weight[ei])
                            tmp_weight_mu2ta2.append(nom_weight[ei])
                        elif muon_pt[ei] >= 100 and tau_pt[ei] < 100:
                            tmp_weight_mu1ta1.append(nom_weight[ei])
                            tmp_weight_mu1ta2.append(nom_weight[ei])
                            tmp_weight_mu2ta1.append(pd_weight[ei])
                            tmp_weight_mu2ta2.append(nom_weight[ei])
                        elif muon_pt[ei] >= 100 and tau_pt[ei] >= 100:
                            tmp_weight_mu1ta1.append(nom_weight[ei])
                            tmp_weight_mu1ta2.append(nom_weight[ei])
                            tmp_weight_mu2ta1.append(nom_weight[ei])
                            tmp_weight_mu2ta2.append(pd_weight[ei])
                        else: print("Wrong event categorization for scale UNC!!")

                    if "up" in syst: typeS = "up"
                    else:            typeS = "down"

                    #do normal histo first
                    dnnhist_Wsyst = np.histogram(pred, bins=binedges, weights=pd_weight, density=False)
                    dnnhist_entries_Wsyst = np.histogram(pred, bins=binedges, density=False)
                    syst_hist_dnn_dict["h_dnn_pred_S5__"+syst] = dnnhist_Wsyst
                    syst_hist_dnn_entries_dict["h_dnn_entries_S5__"+syst] = dnnhist_entries_Wsyst

                    #pt binned unc
                    dnnhist_Wsyst11 = np.histogram(pred, bins=binedges, weights=tmp_weight_mu1ta1, density=False)
                    dnnhist_entries_Wsyst11 = np.histogram(pred, bins=binedges, density=False)
                    syst_hist_dnn_dict["h_dnn_pred_S5__"+syst.replace(typeS, "")+"mu1ta1"+typeS] = dnnhist_Wsyst11
                    syst_hist_dnn_entries_dict["h_dnn_entries_S5__"+syst.replace(typeS, "")+"mu1ta1"+typeS] = dnnhist_entries_Wsyst11
                    syst_extend.append(syst.replace(typeS, "")+"mu1ta1"+typeS)

                    dnnhist_Wsyst12 = np.histogram(pred, bins=binedges, weights=tmp_weight_mu1ta2, density=False)
                    dnnhist_entries_Wsyst12 = np.histogram(pred, bins=binedges, density=False)
                    syst_hist_dnn_dict["h_dnn_pred_S5__"+syst.replace(typeS, "")+"mu1ta2"+typeS] = dnnhist_Wsyst12
                    syst_hist_dnn_entries_dict["h_dnn_entries_S5__"+syst.replace(typeS, "")+"mu1ta2"+typeS] = dnnhist_entries_Wsyst12
                    syst_extend.append(syst.replace(typeS, "")+"mu1ta2"+typeS)

                    dnnhist_Wsyst21 = np.histogram(pred, bins=binedges, weights=tmp_weight_mu2ta1, density=False)
                    dnnhist_entries_Wsyst21 = np.histogram(pred, bins=binedges, density=False)
                    syst_hist_dnn_dict["h_dnn_pred_S5__"+syst.replace(typeS, "")+"mu2ta1"+typeS] = dnnhist_Wsyst21
                    syst_hist_dnn_entries_dict["h_dnn_entries_S5__"+syst.replace(typeS, "")+"mu2ta1"+typeS] = dnnhist_entries_Wsyst21
                    syst_extend.append(syst.replace(typeS, "")+"mu2ta1"+typeS)

                    dnnhist_Wsyst22 = np.histogram(pred, bins=binedges, weights=tmp_weight_mu2ta2, density=False)
                    dnnhist_entries_Wsyst22 = np.histogram(pred, bins=binedges, density=False)
                    syst_hist_dnn_dict["h_dnn_pred_S5__"+syst.replace(typeS, "")+"mu2ta2"+typeS] = dnnhist_Wsyst22
                    syst_hist_dnn_entries_dict["h_dnn_entries_S5__"+syst.replace(typeS, "")+"mu2ta2"+typeS] = dnnhist_entries_Wsyst22
                    syst_extend.append(syst.replace(typeS, "")+"mu2ta2"+typeS)

                else:
                    dnnhist_Wsyst = np.histogram(pred, bins=binedges, weights=pd_weight, density=False)
                    dnnhist_entries_Wsyst = np.histogram(pred, bins=binedges, density=False)
                    syst_hist_dnn_dict["h_dnn_pred_S5__"+syst] = dnnhist_Wsyst
                    syst_hist_dnn_entries_dict["h_dnn_entries_S5__"+syst] = dnnhist_entries_Wsyst
                if 'btag' in syst:
                    hist_nevents_S4_dict["h_nevents_S4__"+syst] = infile_forS["h_nevents_S4__"+syst]
                    hist_nevents_S5_dict["h_nevents_S5__"+syst] = infile["h_nevents_S5__"+syst]

    with uproot.recreate(outf_dir) as outf:
        if '__' not in outf_dir.split('/')[-1]:
            outt= outf.mktree("Events", {"Muon1_pt": np.float64, "Tau1_pt": np.float64, "dnn_pred": np.float64, "dnn_st": np.float64, "dnn_tt": np.float64, "dnn_bkg": np.float64})
            outt.extend({"Muon1_pt": muon_pt, "Tau1_pt": tau_pt, "dnn_pred": pred, "dnn_st": pred_st, "dnn_tt": pred_tt, "dnn_bkg": pred_bkg})
        outf["h_dnn_pred_S5"] = dnnhist_nom
        outf["h_dnn_entries_S5"] = dnnhist_entries_nom
        outf["h_dnn_pred_st_S5"] = dnnhist_st
        outf["h_dnn_entries_st_S5"] = dnnhist_entries_st
        outf["h_dnn_pred_tt_S5"] = dnnhist_tt
        outf["h_dnn_entries_tt_S5"] = dnnhist_entries_tt
        outf["h_dnn_pred_bkg_S5"] = dnnhist_bkg
        outf["h_dnn_entries_bkg_S5"] = dnnhist_entries_bkg
        outf["hcounter"] = hcounter

        if (not "Muon" in input_file) and (not "egamma" in input_file.lower()):
            outf["h_nevents_S4_nobtag"] = h_nevents_S4_nobtag
            outf["h_nevents_S4"] = h_nevents_S4
            outf["h_nevents_S2_nobtag"] = h_nevents_S2_nobtag
            outf["h_nevents_S2"] = h_nevents_S2

            #if any(string in input_file for string in ["Tau-LFV","TTt","_ST_t"]) and "__" not in input_file:
            #    outf["ScaleWeightSum"] = ScaleWeightSum
            #    outf["PSWeightSum"] = PSWeightSum
            #    outf["LHEPdfWeightSum"] = LHEPdfWeightSum

        if len(syst_extend) > 1:
            syst_list.extend(syst_extend)
            print ("syst_list.extend ", syst_extend)

        for syst in syst_list:
            if 'pdf' in syst and not 'alphas' in syst:
                systnum = int(syst.split("pdf")[1])
                if systnum > 100: continue
                systupdated="pdf"+str(systnum)+"up"
                outf["h_dnn_pred_S5__"+systupdated] = syst_hist_dnn_dict["h_dnn_pred_S5__"+syst]
                outf["h_dnn_entries_S5__"+systupdated] = syst_hist_dnn_entries_dict["h_dnn_entries_S5__"+syst]
                systupdated="pdf"+str(systnum)+"down"
                outf["h_dnn_pred_S5__"+systupdated] = dnnhist_nom
                outf["h_dnn_entries_S5__"+systupdated] = dnnhist_entries_nom
            else:
                outf["h_dnn_pred_S5__"+syst] = syst_hist_dnn_dict["h_dnn_pred_S5__"+syst]
                outf["h_dnn_entries_S5__"+syst] = syst_hist_dnn_entries_dict["h_dnn_entries_S5__"+syst]

            if 'btag' in syst:
                outf["h_nevents_S4__"+syst] = hist_nevents_S4_dict["h_nevents_S4__"+syst]
                outf["h_nevents_S5__"+syst] = hist_nevents_S5_dict["h_nevents_S5__"+syst]

    K.clear_session()

if __name__ == '__main__':

    print("Start Multi LFV Evaluation")
    discriminator = "p_st_tt_ob"
    #discriminator = "p_tt"
    alpha=options.alpha
    parameters = []
    xsec_rootfiles = load_xsec_rootfiles(options.xsecfile)
    for year in ["v12_2022", "v12_2022EE", "v12_2023", "v12_2023BPix", "v15_2024"]:
        print(year)
        project_dir = options.project_dir + "/" + ch + "/" + year + "/"
        flist = os.listdir(project_dir)
        flist = [i for i in flist if (".root" in i)]
        if xsec_rootfiles is not None:
            flist = [i for i in flist if i in xsec_rootfiles]
        for curfile in flist:
            parameters.append((year, project_dir + curfile, discriminator, alpha))

    parameters_sorted = [tup for tup in parameters if '__' not in tup[1]]
    parameters_sorted.extend([tup for tup in parameters if '__' in tup[1]])
    print ("params: ", parameters_sorted)
    pool = multiprocessing.get_context("spawn").Pool(8)
    pool.map(run, parameters_sorted)
    pool.close()
    pool.join()


def frange(start, stop, step, n=None):
    if step == 0:
        raise ValueError('step must not be 0')
    # how many decimal places are showing?
    if n is None:
        n = max([0 if '.' not in str(i) else len(str(i).split('.')[1])
                for i in (start, stop, step)])
    if step*(stop - start) > 0:  # a non-null incr/decr range
        if step < 0:
            for i in frange(-start, -stop, -step, n):
                yield -i
        else:
            steps = round((stop - start)/step)
            while round(step*steps + start, n) < stop:
                steps += 1
            for i in range(steps):
                yield round(start + i*step, n)

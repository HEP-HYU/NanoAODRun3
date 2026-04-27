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

gpus = tf.config.list_physical_devices('GPU')
if gpus:
    try:
        tf.config.set_logical_device_configuration(gpus[0], [tf.config.LogicalDeviceConfiguration(memory_limit=1024)])
    except RuntimeError as e:
        print(e)

parser = argparse.ArgumentParser(usage="%prog [options]")
parser.add_argument("-O", "--outdir", dest="outdir", type=str, default="test", help="Output folder in your working directory")
parser.add_argument("-I", "--indir", dest="indir", type=str, default="test", help="Iput folder in your working directory")
parser.add_argument("-C", "--ch", dest="ch", type=str, default="muon")
options = parser.parse_args()
ch = options.ch

def min_max_scaling(series):
    return (series - series.min()) / (series.max() - series.min())


training_path = options.indir
#training_path = ""

inputvars = [#"Muon1_pt", "Muon1_eta",
            "Tau1_pt","Tau1_mass","Tau1_eta",
            "Jet1_mass","Jet1_eta","Jet1_btagPNetB",
            "Jet2_pt","Jet2_mass","Jet2_eta","Jet2_btagPNetB",
            "Jet3_pt","Jet3_mass","Jet3_eta","Jet3_btagPNetB",
            "chi2","chi2_SMW_mass","chi2_SMTop_mass",
            "chi2_wqq_dEta","chi2_wqq_dPhi","chi2_wqq_dR",
            "leptau_mass","leptau_dEta","leptau_dPhi","leptau_dR",
            "PuppiMET_pt"
            ]
if ch == "muon": inputvars = ["Muon1_pt", "Muon1_eta"] + inputvars
else: inputvars = ["Electron1_pt", "Electron1_eta"] + inputvars


#discriminators = {"p_st" : 2, "p_tt" : 1 , "p_bkg" : 0 , "p_st_tt" : 999, "p_st_tt_ob" : 999 }

def run(inputs):
    year = inputs[0]
    input_file = inputs[1]
    #discriminator_key = inputs[2]
    #alpha = inputs[3]
    alpha = inputs[2]
    #print(year, discriminator_key , alpha)

    #binedges = [0,1,2,3,5,10,30,60]
    binedges = [i for i in frange(0.0, 100.01, 0.01)]

    model_dir = os.path.join(training_path, ch+"/nom/best_model.keras")
    model = tf.keras.models.load_model(model_dir)

    hists_path = os.path.join(options.outdir, year)
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
    syst_hist_dnn_dict = {}
    syst_hist_dnn_entries_dict = {}
    hist_nevents_S4_dict = {}
    hist_nevents_S5_dict = {}
    syst_extend = []

    if not "Muon" in input_file :
        infile_forS = uproot.open(input_file.replace("_FF", ""))
        h_nevents_S4_nobtag = infile_forS["h_nevents_S4_nobtag"]  ### get it from removed FF 
        h_nevents_S4 = infile_forS["h_nevents_S4"]
        h_nevents_S2_nobtag = infile_forS["h_nevents_S2_nobtag"]
        h_nevents_S2 = infile_forS["h_nevents_S2"]
        if any(string in input_file for string in ["Tau-LFV", "TTt", "_ST_t"]) and "__" not in input_file:
            ScaleWeightSum = infile['ScaleWeightSum']
            PSWeightSum = infile['PSWeightSum']
            LHEPdfWeightSum = infile['LHEPdfWeightSum']

    if nEvents == 0:
        #print("No events : "+input_file)
        #Need to add empth histograms for technical reasons ....
        muon_pt = []
        tau_pt = []
        pred = []
        pd_weight = []
        dnnhist_nom = np.histogram(pred, bins=binedges, weights=pd_weight, density=False)
        dnnhist_entries_nom = np.histogram(pred, bins=binedges, density=False)
        print ("syst_list: ", syst_list)
        for syst in syst_list:
            syst_hist_dnn_dict["h_dnn_pred_S5__"+syst] = dnnhist_nom
            syst_hist_dnn_entries_dict["h_dnn_entries_S5__"+syst] = dnnhist_entries_nom
            if 'btag' in syst:
                hist_nevents_S4_dict["h_nevents_S4__"+syst] = infile_forS["h_nevents_S4__"+syst]
                hist_nevents_S5_dict["h_nevents_S5__"+syst] = infile["h_nevents_S5__"+syst]
    else:
        muon_pt = 0
        if ch == "muon": muon_pt = tree["Muon1_pt"].array()
        else: muon_pt = tree["Electron1_pt"].array()
        tau_pt = tree["Tau1_pt"].array()
        pd_data = tree.arrays(inputvars, library="pd")
        pd_weight = tree.arrays(weights, library="np")
        pred_data = np.array(pd_data.filter(items = inputvars))

        pred = model.predict(pred_data, batch_size=128)
        #pred_df = pd.DataFrame(pred, columns=['Prediction1', 'Prediction2', 'Prediction3'])
        #print("PRED SHAPE : ", pred.shape)
        #result_df = pd.concat([pd_data, pred_df], axis=1)
        #result_df.to_csv(outf_dir.split(".root")[0]+'combined_data.csv', index=False)
        #print("PRED COMB SHAPE : ", result_df.shape)
        #print("DF saved here: ", outf_dir.split(".root")[0]+'combined_data.csv')

        rmzeros = pred[:,0]
        rmzeros[rmzeros <= 0.001] = 0.001
        pred[:,0] = rmzeros
        pred = ( ((1 - alpha) * pred[:,2] + alpha * pred[:,1]) / pred[:,0])
        pred[pred >= 100.0] = 99.999
        pred = pred.tolist()
        pd_weight = pd_weight[eventWeight].tolist()
        nom_weight = pd_weight.copy()
        dnnhist_nom = np.histogram(pred, bins=binedges, weights=pd_weight, density=False)
        dnnhist_entries_nom = np.histogram(pred, bins=binedges, density=False)

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
            outt= outf.mktree("Events", {"Muon1_pt": np.float64, "Tau1_pt": np.float64, "dnn_pred": np.float64})
            outt.extend({"Muon1_pt": muon_pt, "Tau1_pt": tau_pt, "dnn_pred": pred})
        outf["h_dnn_pred_S5"] = dnnhist_nom
        outf["h_dnn_entries_S5"] = dnnhist_entries_nom
        outf["hcounter"] = hcounter

        if not "SingleMuon" in input_file:
            outf["h_nevents_S4_nobtag"] = h_nevents_S4_nobtag
            outf["h_nevents_S4"] = h_nevents_S4

            if any(string in input_file for string in ["Tau-LFV","TTt","_ST_t"]) and "__" not in input_file:
                outf["ScaleWeightSum"] = ScaleWeightSum
                outf["PSWeightSum"] = PSWeightSum
                outf["LHEPdfWeightSum"] = LHEPdfWeightSum

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
    #discriminator = "p_st_tt_ob"
    alpha=0.1
    parameters = []
    #for year in ["v2022", "v2022EE", "v2023", "v2023_BPix", "v15_2024"]:
    for year in ["v15_2024"]:
        print(year)
        project_dir = "/home/itseyes/github/NanoAODRun3/LFVAnalyzer/process_0323_btag_v3/" + ch + "/" + year + "/"
        flist = os.listdir(project_dir)
        flist = [i for i in flist if (".root" in i)]
        for curfile in flist:
            #parameters.append((year, project_dir + curfile, discriminator, alpha))
            parameters.append((year, project_dir + curfile, alpha))

    parameters_sorted = [tup for tup in parameters if '__' not in tup[1]]
    parameters_sorted.extend([tup for tup in parameters if '__' in tup[1]])
    print ("params: ", parameters_sorted)
    pool = multiprocessing.get_context("spawn").Pool(12)
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

from __future__ import print_function
from ROOT import *
import ROOT
import sys, os
import ast, array
from collections import OrderedDict
gROOT.SetBatch(True)
hostname = os.environ["HOSTNAME"]

bkg_ds = "TTTo2L2Nu"
years = ["2016pre", "2016post", "2017", "2018"]
dss = ["TTTo2L2Nu", "TT_LFV_TCMuTau_Scalar", "TT_LFV_TCMuTau_Vector", "TT_LFV_TCMuTau_Tensor", "TT_LFV_TUMuTau_Scalar", "TT_LFV_TUMuTau_Vector", "TT_LFV_TUMuTau_Tensor"]
datasets = OrderedDict()
filelists = OrderedDict()
c = TCanvas('c','c',400,400)

for year in years:
    for ds in dss:
        htmp = TH1F(ds, ds, 40, 0, 2000)
        arr = array.array('d',[0,50,100,150,200,250,300,350,400,450,500,550,
                    600,800,1000,2000])
        #            600,650,700,800,900,1000,1100,1200,1400,1600,1800,2000])
        htmp = htmp.Rebin(len(arr)-1, htmp.GetName(), arr)
        htmp.GetXaxis().SetTitle('top quark p_{T} (GeV)')
        htmp.SetStats(0)
        datasets[year + '_' + ds] = [htmp, 831.76]

        #datasets['tthut'] = [hhut, 496.1] #From xsdb
        #datasets['tthct'] = [hhct, 496.1]


    for name, data in datasets.items():
      entries = 0
      with open(name + '.txt', 'r') as f:
        bins = f.readline(); #Should be first line
        bins = ast.literal_eval(bins)
        for i in range(len(bins)):
          h = data[0]
          h.SetBinContent(i+1, bins[i])
        entries = f.readline();
      if float(entries) > 0:
          h.Scale(data[1]/float(entries))
      h.Draw('hist e')
      c.Print(name+'.pdf')

    c.Clear()


    hratio = []
    hbkg = datasets[year + '_' + bkg_ds][0]

    for key, value in datasets.items():
        if 'TTTo' in key: continue

        hbkg_sig = hbkg.Clone('hbkg_sig')
        hbkg_sig.SetDirectory(0)

        hbkg_sig.Divide(value[0])
        hbkg_sig.Draw('hist e')
        c.Print(key+'.pdf')
        hratio.append(hbkg_sig)
        print('Integral of Powheg+P8 divided by ' + key, hbkg_sig.Integral()/hbkg_sig.GetNbinsX())
        tmp_string = ''
        for i in range(1,hbkg_sig.GetNbinsX()+1):
          tmp_string += '{'+str(hbkg_sig.GetBinCenter(i))+', '+str(round(hbkg_sig.GetBinContent(i),6))+'},'
        with open(key + '_ratio.txt', 'w') as f: f.write(tmp_string)

        c.Clear()

        l = TLegend(0.3,0.71,0.7,0.87)
        l.SetTextSize(0.03)
        l.SetLineColor(0)
        l.SetFillColor(0)
        l.AddEntry(hbkg_sig, 'Powheg+P8/' + key.replace('TT_LFV_',''), 'l')

        hbkg_sig.SetTitle('Ratio between MG5 signal and Powheg NLO tt+LL')
        hbkg_sig.SetLineColor(2)
        #hbkg_sig.Fit('pol4')
        #hbkg_sig.GetYaxis().SetRangeUser(0.9,1.9)
        hbkg_sig.GetYaxis().SetRangeUser(0.8,1.4)
        hbkg_sig.Draw('ep')
        l.Draw('same')
        c.Print(key + '_ratio.pdf')

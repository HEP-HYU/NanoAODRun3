from ROOT import *
import ROOT
import sys, os
import ast, array
from collections import OrderedDict
gROOT.SetBatch(True)
hostname = os.environ["HOSTNAME"]

bkg_ds = 'TTTo2L2Nu'
years = ["2016pre", "2016post", "2017", "2018"]
ranks = ["Scalar", "Vector", "Tensor"]
coupls = ["TCMuTau", "TUMuTau"]
dss = ["TTTo2L2Nu", "TT_LFV_COUP_RANK"]
filelists = OrderedDict()

h = TH1F('h', 'h', 40, 0, 2000)
arr = array.array('d',[0,50,100,150,200,250,300,350,400,450,500,550,
            600,800,1000,2000])
#            600,650,700,800,900,1000,1100,1200,1400,1600,1800,2000])
h = h.Rebin(len(arr)-1, h.GetName(), arr)
h.GetXaxis().SetTitle('top quark p_{T} (GeV)')
h.SetStats(0)

hbkg17 = h.Clone('hbkg17')
hbkg17.SetDirectory(0)

c = TCanvas('c','c',400,400)
sum_weight = [0] * (len(arr)-1)

for coup in coupls:
    for rank in ranks:
        hs = []
        hb = {}
        l = TLegend(0.3,0.71,0.7,0.87)
        l.SetTextSize(0.03)
        l.SetLineColor(0)
        l.SetFillColor(0)

        for year in years:
            entries = 0

            with open(year + '_' + bkg_ds +'.txt', 'r') as f:
                bins = f.readline(); #Should be first line
                bins = ast.literal_eval(bins)
                htmp = h.Clone('htmp')
                for i in range(len(bins)):
                    htmp.SetBinContent(i+1, bins[i])
                entries = f.readline();
                htmp.Scale(831.76/float(entries))
                hb[year + '_' + bkg_ds] = htmp

            with open(year + '_TT_LFV_' + coup + '_' + rank +'.txt', 'r') as f:
                bins = f.readline(); #Should be first line
                bins = ast.literal_eval(bins)
                htmp = h.Clone('htmp')
                for i in range(len(bins)):
                    htmp.SetBinContent(i+1, bins[i])
                entries = f.readline();
                if int(entries) > 0:
                    htmp.Scale(831.76/float(entries))
                    htmp2 = hb[year + '_' + bkg_ds].Clone('htmp2')
                    htmp2.Divide(htmp)
                    hs.append(htmp2)
                    l.AddEntry(htmp2, 'Powheg+P8/' + year + '_' + coup + '_' + rank, 'l')

        col = 1
        hs[0].SetTitle('Ratio between MG5 signal and Powheg NLO tt+LL')
        hs[0].SetLineColor(col)
        hs[0].SetLineWidth(2)
        hs[0].GetYaxis().SetRangeUser(0.5,1.4)
        hs[0].Draw('ep')

        for nbin in range(hs[0].GetNbinsX()):
            sum_weight[nbin] += hs[0].GetBinContent(nbin+1)

        for i in range(1, len(hs)):
            col += 1
            hs[i].SetLineColor(col)
            hs[i].SetLineWidth(2)
            hs[i].Draw('ep same')

            for nbin in range(hs[i].GetNbinsX()):
                sum_weight[nbin] += hs[i].GetBinContent(nbin+1)

        l.Draw('same')
        c.Print(coup + '_' + rank + '_ratios.pdf')


average_weight = []
for i in sum_weight:
    average_weight.append(round(i/24.0, 4))

print(average_weight)

#!/usr/bin/env python
import os, sys
from ROOT import *
import ROOT
from array import array
gROOT.SetBatch(True)
import numpy as np

draw_tau = False
#draw_tau = True
if draw_tau:
    pdgid = 15
    lepname = 'Tau'
    histname = 'h_tau1_pt_S5'
else:
    pdgid = 13
    lepname = 'Muon'
    histname = 'h_muon1_pt_S5'

bare_folder = 'v9_0714_1010'
rewgt_folder = 'v9_0714_1010_uforeweight'

year_list = {'2016pre':'19.5', '2016post':16.8, '2017':41.5, '2018':59.8, 'Run2':138}

file_list = ['hist_ST_LFV_TCMuTau_Scalar.root', 'hist_ST_LFV_TUMuTau_Scalar.root',
             'hist_ST_LFV_TCMuTau_Vector.root', 'hist_ST_LFV_TUMuTau_Vector.root',
             'hist_ST_LFV_TCMuTau_Tensor.root', 'hist_ST_LFV_TUMuTau_Tensor.root']

target_dir = 'ufo_reweight_comparison'
os.makedirs(target_dir, exist_ok=True)

for year, lumi in year_list.items():
    for file_name in file_list:

        f_smeftfr = TFile.Open(os.path.join(bare_folder, year, file_name), 'R')
        f_smeftsim = TFile.Open(os.path.join(rewgt_folder, year, file_name), 'R')

        hnowgt = f_smeftfr.Get(histname)
        hrewgt = f_smeftsim.Get(histname)

        hnowgt.AddBinContent(hnowgt.GetNbinsX(), hnowgt.GetBinContent(hnowgt.GetNbinsX()))
        hrewgt.AddBinContent(hrewgt.GetNbinsX(), hrewgt.GetBinContent(hrewgt.GetNbinsX()))

        #different # gen events, only take relative ratio
        sum_nowgt = hnowgt.Integral(0,-1)
        sum_rewgt = hrewgt.Integral(0,-1)

        print(sum_nowgt, sum_rewgt)

        c = TCanvas('c', 'c', 450, 450)
        pad1 = TPad('pad1', 'pad1', 0.0, 0.3, 1, 1.0)
        pad1.SetBottomMargin(0.02)
        pad1.Draw()
        c.cd()
        pad2 = TPad('pad2', 'pad2', 0.0, 0.0, 1, 0.28)
        pad2.SetBottomMargin(0.3)
        pad2.SetTopMargin(0.02)
        pad2.SetGridx()
        pad2.SetGridy()
        pad2.Draw()

        pad1.cd()
        #pad1.SetLogy()
        hrewgt.SetMarkerStyle(7)
        hrewgt.SetStats(0)
        hrewgt.SetTitle('')
        hrewgt.GetXaxis().SetLabelSize(0)
        hrewgt.GetXaxis().SetTitleOffset(5.0)
        #hrewgt.GetYaxis().SetRangeUser(5, hnowgt.GetMaximum()*20)
        hrewgt.GetYaxis().SetRangeUser(0, hnowgt.GetMaximum()*1.1)
        hrewgt.GetYaxis().SetTitle('A.U.')
        hrewgt.GetYaxis().SetTitleOffset(1.1)
        hrewgt.GetYaxis().SetTitleSize(0.045)
        hrewgt.DrawNormalized('ep')

        hnowgt.SetLineColor(ROOT.kRed+2)
        hnowgt.DrawNormalized('same hist')

        l = TLegend(0.4,0.6,0.83,0.87)
        l.SetNColumns(1);
        l.SetTextSize(0.04);
        l.SetLineColor(0);
        l.SetFillColor(0);
        l.AddEntry(hnowgt, 'SMEFTfr', 'f')
        l.AddEntry(hrewgt, 'SMEFTfr, reweighted', 'p')
        l.AddEntry('', 'SumW SMEFTfr rewgt / unwgt', '')
        l.AddEntry('', ' = {}'.format(str(round(sum_rewgt / sum_nowgt,2))), '')
        l.Draw('same')

        label = TPaveText()
        label.SetX1NDC(gStyle.GetPadLeftMargin())
        label.SetY1NDC(1.0-gStyle.GetPadTopMargin())
        label.SetX2NDC(1.0-gStyle.GetPadRightMargin()+0.03)
        label.SetY2NDC(1.0)
        label.SetTextFont(62)
        label.AddText("CMS simulation                         {}".format(lumi) + " fb^{-1} at #sqrt{s} = 13 TeV")
        label.SetFillStyle(0)
        label.SetBorderSize(0)
        label.SetTextSize(0.05)
        label.SetTextAlign(32)
        label.Draw("same")

        pad2.cd()

        hrewgt_r = hrewgt.Clone('hrewgt_r')
        hnowgt_r = hnowgt.Clone('hnowgt_r')
        hrewgt_r.Scale(1/sum_rewgt)
        hnowgt_r.Scale(1/sum_nowgt)

        hnowgt_r.Divide(hrewgt_r)
        hnowgt_r.SetLineColor(1)

        hnowgt_r.SetStats(0)
        hnowgt_r.SetTitle('')
        hnowgt_r.GetXaxis().SetTitle('S5 Reco {} pT (GeV)'.format(lepname))
        hnowgt_r.GetXaxis().SetTitleSize(0.11)
        hnowgt_r.GetXaxis().SetTitleOffset(1.0)
        hnowgt_r.GetXaxis().SetLabelSize(0.1)
        hnowgt_r.GetYaxis().SetTitle('Divided by rewgt. SMEFTfr')
        hnowgt_r.GetYaxis().SetRangeUser(0.4,1.6)
        hnowgt_r.GetYaxis().SetNdivisions(505)
        hnowgt_r.GetYaxis().SetTitleSize(0.08)
        hnowgt_r.GetYaxis().SetTitleOffset(0.45)
        hnowgt_r.GetYaxis().SetLabelSize(0.1)

        hnowgt_r.Draw('hist e')

        c.Print(target_dir + '/comp_ufo_{}.pdf'.format(file_name.replace('hist_','').replace('.root','')+'_'+lepname+'_'+year))

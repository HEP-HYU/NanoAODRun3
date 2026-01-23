import ROOT
from ROOT import *
import tdrstyle
import CMS_lumi

runs = ["16pre","16post","17","18"]
datafiles = ["PileupDATA_UL16pre","PileupDATA_UL16post","PileupDATA_UL17","PileupDATA_UL18"]
mcfiles = ["PileupMC_UL16","PileupMC_UL17","PileupMC_UL18"]
tdrstyle.setTDRStyle()
ROOT.gROOT.SetStyle("tdrStyle")

for run in runs:
    datarootfile = "PileupDATA_UL"+run+".root"
    mcrootfile = ""
    if "16" in run:
        mcrootfile = "PileupMC_UL16.root"
    else:
        mcrootfile = "PileupMC_UL"+run+".root"
    
    if run == "16pre":
        CMS_lumi.lumi_13TeV = "19.5 fb^{-1}"
    if run == "16post":
        CMS_lumi.lumi_13TeV = "16.8 fb^{-1}"
    elif run == "17":
        CMS_lumi.lumi_13TeV = "41.5 fb^{-1}"
    elif run == "18":
        CMS_lumi.lumi_13TeV = "59.8 fb^{-1}"

    fdata = TFile(datarootfile)
    h_pu = fdata.Get("pileup")
    h_pu_plus = fdata.Get("pileup_plus")
    h_pu_minus = fdata.Get("pileup_minus")
    fmc = TFile(mcrootfile)
    fmc.cd()
    h_pu_mc = fmc.Get("pu_mc")

    h_pu.SetLineStyle(1)
    h_pu_plus.SetLineStyle(3)
    h_pu_minus.SetLineStyle(7)
    h_pu_mc.SetLineStyle(1)
    h_pu.SetLineColor(1)
    h_pu_plus.SetLineColor(1)
    h_pu_minus.SetLineColor(1)
    h_pu_mc.SetLineColor(2)
    h_pu.SetLineWidth(2)
    h_pu_plus.SetLineWidth(2)
    h_pu_minus.SetLineWidth(2)
    h_pu_mc.SetLineWidth(2)

    h_pu.Scale(1/h_pu.Integral())
    h_pu_plus.Scale(1/h_pu_plus.Integral())
    h_pu_minus.Scale(1/h_pu_minus.Integral())
    h_pu_mc.Scale(1/h_pu_mc.Integral())
    
    m = max([h_pu.GetMaximum(),h_pu_plus.GetMaximum(),h_pu_minus.GetMaximum(),h_pu_mc.GetMaximum()])
    h_pu.SetMaximum(m*1.2)

    c1 = TCanvas('c','c',600,600)
    c1.SetTopMargin(0.1)
    c1.SetBottomMargin(0.13)
    c1.SetRightMargin(0.05)
    c1.SetLeftMargin(0.15)
    gStyle.SetOptStat(0)
    gStyle.SetOptTitle(0)

    h_pu.Draw("hist")
    h_pu_plus.Draw("histsame")
    h_pu_minus.Draw("histsame")
    h_pu_mc.Draw("histsame")

    h_pu.GetXaxis().SetTitle("True Number of Interactions")
    h_pu.GetXaxis().SetTitleSize(0.05)
    h_pu.GetYaxis().SetTitle("Probability")
    h_pu.GetYaxis().SetTitleSize(0.05)
    
    leg = TLegend(0.6,0.6,0.9,0.85)
    gStyle.SetLegendTextSize(0.03)
    leg.SetBorderSize(0)
    leg.AddEntry(h_pu,"Pileup Data")
    leg.AddEntry(h_pu_plus,"Pileup Data up")
    leg.AddEntry(h_pu_minus,"Pileup Data down")
    leg.AddEntry(h_pu_mc,"Pileup MC")
    leg.Draw()
    
    CMS_lumi.extraText = "Work in Progress"
    CMS_lumi.CMS_lumi(c1, 4, 0)
    c1.cd()
    c1.Update()
    c1.RedrawAxis()
     
    c1.Print("h_pileup"+run+".pdf")
    c1.Close()
    fdata.Close()
    fmc.Close()

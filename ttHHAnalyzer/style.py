#! /usr/bin/env python 

from ROOT import gROOT, TFile, TCanvas, TF1, TGraph, TPaveText, TLatex, TH1F, TLegend, TLine, TColor, TIter, gStyle
from ROOT import TGaxis
from array import array
import ROOT

gROOT.SetStyle("Plain")
gStyle.SetOptStat(1110)
gStyle.SetOptFit(1)
gStyle.SetOptStat(0)
gStyle.SetStatW(0.25)
gStyle.SetStatH(0.15)

gStyle.SetCanvasDefH(150)
gStyle.SetCanvasDefW(150)

gStyle.SetAxisColor(1, "XYZ")
gStyle.SetStripDecimals(1)
gStyle.SetTickLength(0.03, "XYZ")
gStyle.SetNdivisions(510, "XYZ")
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)

gStyle.SetPadTopMargin(0.1)
gStyle.SetPadBottomMargin(0.13)
gStyle.SetPadLeftMargin(0.13)
gStyle.SetPadRightMargin(0.1)

gStyle.SetTitleColor(1, "XYZ")
gStyle.SetTitleFont(42, "XYZ")
gStyle.SetTitleSize(0.035, "XYZ")
gStyle.SetTitleXOffset(1.1)
gStyle.SetTitleYOffset(1.55)

gStyle.SetLabelColor(1, "XYZ")
gStyle.SetLabelFont(42, "XYZ")
gStyle.SetLabelOffset(0.007, "XYZ")
gStyle.SetLabelSize(0.035, "XYZ")

gStyle.SetLegendTextSize(0.035)
gStyle.SetLegendBorderSize(0)

gStyle.SetErrorX(0)
gROOT.ForceStyle()

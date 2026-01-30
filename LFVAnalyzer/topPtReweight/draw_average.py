from __future__ import print_function
from ROOT import *
import ROOT
import sys, os
import ast, array
gROOT.SetBatch(True)

hbkg = TH1F('hbkg', 'hbkg', 40, 0, 2000)
arr = array.array('d',[0,50,100,150,200,250,300,350,400,450,500,550,
            600,800,1000,2000])
hbkg = hbkg.Rebin(len(arr)-1, hbkg.GetName(), arr)
hbkg.GetXaxis().SetTitle('top quark p_{T} (GeV)')
hbkg.SetStats(0)
hbkg2 = hbkg.Clone('hbkg2')
hsig  = hbkg.Clone('hsig')
hsig2 = hbkg.Clone('hsig2')
hsig3 = hbkg.Clone('hsig3')
hsig4 = hbkg.Clone('hsig4')

c = TCanvas('c','c',400,300)

with open('V9_7_ttbkg.txt', 'r') as f:
  bins = f.readline() #Should be first line
  bins = ast.literal_eval(bins)
  entries = f.readline()
  for i in xrange(len(bins)):
    hbkg.SetBinContent(i+1, bins[i])
  hbkg.Scale(831.76/float(entries))

with open('V10_3_ttbkg.txt', 'r') as f:
  bins = f.readline() #Should be first line
  bins = ast.literal_eval(bins)
  entries = f.readline()
  for i in xrange(len(bins)):
    hbkg2.SetBinContent(i+1, bins[i])
  hbkg2.Scale(831.76/float(entries))

hbkg.Add(hbkg2 ,1.0)
hbkg.Scale(0.5)

with open('V9_7_tthct.txt', 'r') as f:
  bins = f.readline() #Should be first line
  bins = ast.literal_eval(bins)
  entries = f.readline()
  for i in xrange(len(bins)):
    hsig.SetBinContent(i+1, bins[i])
  hsig.Scale(496.1/float(entries))

with open('V9_7_tthut.txt', 'r') as f:
  bins = f.readline() #Should be first line
  bins = ast.literal_eval(bins)
  entries = f.readline()
  for i in xrange(len(bins)):
    hsig2.SetBinContent(i+1, bins[i])
  hsig2.Scale(496.1/float(entries))

with open('V10_3_tthct.txt', 'r') as f:
  bins = f.readline() #Should be first line
  bins = ast.literal_eval(bins)
  entries = f.readline()
  for i in xrange(len(bins)):
    hsig3.SetBinContent(i+1, bins[i])
  hsig3.Scale(496.1/float(entries))

with open('V10_3_tthut.txt', 'r') as f:
  bins = f.readline() #Should be first line
  bins = ast.literal_eval(bins)
  entries = f.readline()
  for i in xrange(len(bins)):
    hsig4.SetBinContent(i+1, bins[i])
  hsig4.Scale(496.1/float(entries))

hsig.Add(hsig2, 1.0)
hsig.Add(hsig3, 1.0)
hsig.Add(hsig4, 1.0)
hsig.Scale(0.25)

hbkg_sig = hbkg.Clone('hbkg_sig')
hbkg_sig.SetDirectory(0)

hbkg_sig.Divide(hsig)

c.Clear()
l = TLegend(0.3,0.71,0.7,0.87)
l.SetTextSize(0.03)
l.SetLineColor(0)
l.SetFillColor(0)
l.AddEntry(hbkg_sig, 'Powheg+P8/LO signal', 'l')

hbkg_sig.SetTitle('Ratio between MG5 LO signal and Powheg NLO tt+LJ')
hbkg_sig.GetYaxis().SetRangeUser(0.9,1.9)
hbkg_sig.Draw('ep')
l.Draw('same')
c.Print('weight_average.pdf')

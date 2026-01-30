import os
from uncertainties import ufloat
import ctypes
from ROOT import *
import numpy as np
import array
from symmetrize import smoothing
gROOT.SetBatch(True)

input_ver = 'v9_0714_fake'
isQcd = False

yield_file = 'yields.tex'
if isQcd: yield_file = 'qcd/' + yield_file

year_list = ['2016pre', '2016post', '2017', '2018']
#no_clean = True
no_clean = False

#pt binning: 20, 0, 400
hist_name = 'h_tau1_pt_S3'
hist_gen_name = 'h_tau1_gen_pt_S3'
#bin_edges = [[1, -1]]
#bin_edges = [[1, 7], [8, -1]] #(0-140) (140-inf)
#bin_edges = [[1, 3], [4, 5], [6, 7], [8, -1]] #(0-60) (60-100) (100-140) (140-inf)


#los = loose tau, OS mu tau (C)
#lss = loose tau, SS mu tau (A)
#tos = tight tau, OS mu tau (D)
#tss = tight tau, SS mu tau (B)

sf_file = TFile.Open('ff_plots/test_sf.root', 'RECREATE')

for year in year_list:

    print("Year: ", year)

    init_dict = {}
    init_dict.fromkeys(['gentau_total', 'mc_total', 'ndata', 'nsub_mc', 'nsub_data'])

    nevents_region = {
        'los': {key: value[:] for key, value in init_dict.items()},
        'lss': {key: value[:] for key, value in init_dict.items()},
        'tos': {key: value[:] for key, value in init_dict.items()},
        'tss': {key: value[:] for key, value in init_dict.items()}
    }

    for region, nevents in nevents_region.items():

        plot_file_path = os.path.join(input_ver + '_' + region, 'figure_'+year, 'plots.root')
        if no_clean and (region == 'los' or region == 'lss'):
            plot_file_path = os.path.join(input_ver + '_' + region + '_org', 'figure_'+year, 'plots.root')
        plot_file = TFile.Open(plot_file_path)

        c_gen = plot_file.Get(hist_gen_name)
        p_gen = c_gen.GetPrimitive('pad_hi')
        s_gen = p_gen.GetPrimitive('mc_stack_0')
        h_gen = s_gen.GetStack().Last()

        c_all = plot_file.Get(hist_name)
        p_all = c_all.GetPrimitive('pad_hi')
        s_all = p_all.GetPrimitive('mc_stack_0')
        h_all = s_all.GetStack().Last()

        h_data = None
        prim_all = p_all.GetListOfPrimitives()
        for prim in prim_all:
            if prim.InheritsFrom("TH1"):
                h_data = prim

        nevents['gentau_total'] = h_gen
        nevents['mc_total'] = h_all
        nevents['ndata'] = h_data
        h_gen.SetDirectory(0)
        h_all.SetDirectory(0)
        h_data.SetDirectory(0)
        h_nsub_mc = h_all.Clone('h_nsub_mc')
        h_nsub_data = h_data.Clone('h_nsub_data')
        h_nsub_mc.Add(h_gen, -1.0)
        h_nsub_data.Add(h_gen, -1.0)
        h_nsub_mc.SetDirectory(0)
        h_nsub_data.SetDirectory(0)
        nevents['nsub_mc'] = h_nsub_mc
        nevents['nsub_data'] = h_nsub_data

    #print(nevents_region)

    c = TCanvas('c', 'c', 400, 400)

    legend = TLegend(0.1, 0.75, 0.3, 0.9)
    legend.SetFillStyle(0)

    ss_mc = nevents_region['tss']['nsub_mc'].Clone('ss_mc')
    ss_mc.Divide(nevents_region['lss']['nsub_mc'])
    ss_mc.SetLineWidth(2)
    ss_mc.SetLineColor(2)
    ss_mc.SetFillColor(0)
    ss_mc.GetYaxis().SetTitle('MisID probability')
    ss_mc.GetYaxis().SetRangeUser(0, 1)
    ss_mc.Draw('b e')

    os_mc = nevents_region['tos']['nsub_mc'].Clone('os_mc')
    os_mc.Divide(nevents_region['los']['nsub_mc'])
    os_mc.SetLineWidth(2)
    os_mc.SetLineColor(4)
    os_mc.SetFillColor(0)
    os_mc.Draw('b e same')

    legend.AddEntry(os_mc, 'OS MC', 'l')
    legend.AddEntry(ss_mc, 'SS MC', 'l')
    legend.Draw('same')

    latexLabel = TLatex()
    latexLabel.SetNDC()
    latexLabel.SetTextFont(62) # helvetica bold face
    latexLabel.SetTextSize(0.4 * c.GetTopMargin())
    latexLabel.DrawLatex(0.11, 0.92, "CMS private work")
    latexLabel.SetTextSize(0.4 * c.GetTopMargin())
    latexLabel.SetTextAlign(31)
    latexLabel.DrawLatex(0.9, 0.92, "%s"%(year))
    latexLabel.Clear()

    c.Print('ff_plots/misid_prob_os_vs_ss_' + year + '.png')
    c.Print('ff_plots/misid_prob_os_vs_ss_' + year + '.pdf')
    legend.Clear()
    c.Clear()


    #############
    # ss mc vs data
    #############
    ss_data = nevents_region['tss']['nsub_data'].Clone('ss_data')
    ss_data.Divide(nevents_region['lss']['nsub_data'])
    ss_data.SetMarkerStyle(1)
    ss_data.SetLineWidth(2)
    ss_data.SetLineColor(4)
    ss_data.SetFillColor(0)
    #ss_data.Rebin(2)
    #ss_data.Fit('pol1')
    #ss_data.GetFunction("pol1").SetLineColor(4)
    ss_data.GetYaxis().SetTitle('MisID probability')
    ss_data.GetYaxis().SetRangeUser(0, 1)
    ss_data.Draw('b e')

    ss_mc = nevents_region['tss']['nsub_mc'].Clone('ss_mc')
    ss_mc.Divide(nevents_region['lss']['nsub_mc'])
    ss_mc.SetLineWidth(2)
    ss_mc.SetLineColor(2)
    ss_mc.SetFillColor(0)
    #ss_mc.Rebin(2)
    ss_mc.Draw('b e same')

    legend.AddEntry(ss_mc, 'SS MC', 'l')
    legend.AddEntry(ss_data, 'SS Data', 'l')
    legend.Draw('same')

    latexLabel = TLatex()
    latexLabel.SetNDC()
    latexLabel.SetTextFont(62) # helvetica bold face
    latexLabel.SetTextSize(0.4 * c.GetTopMargin())
    latexLabel.DrawLatex(0.11, 0.92, "CMS private work")
    latexLabel.SetTextSize(0.4 * c.GetTopMargin())
    latexLabel.SetTextAlign(31)
    latexLabel.DrawLatex(0.9, 0.92, "%s"%(year))
    latexLabel.Clear()

    c.Print('ff_plots/misid_prob_ss_mc_vs_data_' + year + '.png')
    c.Print('ff_plots/misid_prob_ss_mc_vs_data_' + year + '.pdf')
    legend.Clear()
    c.Clear()


    ###############
    # store unbinned hist
    ###############
    os_data = nevents_region['los']['nsub_data'].Clone('os_data')
    ss_data = nevents_region['tss']['nsub_data'].Clone('ss_data')
    lss_data = nevents_region['lss']['nsub_data'].Clone('lss_data')
    ss_data.Divide(lss_data)

    tos_mc = nevents_region['tos']['nsub_mc'].Clone('tos_mc')
    os_data.Multiply(ss_data)
    os_data.Divide(tos_mc)

    sf_file.cd()
    os_data.SetName(year)
    os_data.Write()


    ###############
    # rebinned SF
    ###############
    arr = array.array('d',[0., 40, 120, 200, 280, 400])
    #arr = array.array('d',[0., 100, 200, 400])

    os_data = nevents_region['los']['nsub_data'].Clone('os_data')
    os_data = os_data.Rebin(len(arr)-1, os_data.GetName(), arr)

    ss_data = nevents_region['tss']['nsub_data'].Clone('ss_data')
    ss_data = ss_data.Rebin(len(arr)-1, ss_data.GetName(), arr)

    lss_data = nevents_region['lss']['nsub_data'].Clone('lss_data')
    lss_data = lss_data.Rebin(len(arr)-1, lss_data.GetName(), arr)
    ss_data.Divide(lss_data)

    tos_mc = nevents_region['tos']['nsub_mc'].Clone('tos_mc')
    tos_mc = tos_mc.Rebin(len(arr)-1, tos_mc.GetName(), arr)

    os_data.Multiply(ss_data)
    os_data.Divide(tos_mc)

    for i in range(os_data.GetNbinsX()):
        print("Bin", str(i+1), os_data.GetBinContent(i+1))

    os_data.SetMarkerStyle(1)
    os_data.SetLineWidth(2)
    os_data.SetLineColor(4)
    os_data.SetFillColor(0)
    os_data.Fit('pol1')
    os_data.GetFunction("pol1").SetLineColor(4)
    os_data.SetStats(1)
    gStyle.SetOptStat(000000000)
    gStyle.SetOptFit(11)
    os_data.GetYaxis().SetRangeUser(0, 2)
    os_data.GetYaxis().SetTitle('MisID SF')

    stats1 = os_data.GetListOfFunctions().FindObject("stats").Clone("stats1")
    c.Modified()
    c.Update();
    stats = c.GetPrimitive("stats")
    stats.SetName("Fit stats");
    stats.SetX1NDC(.6);
    stats.SetX2NDC(.9);
    stats.SetY1NDC(.8);
    stats.SetY2NDC(.9);

    os_data.Draw('b e')

    legend.AddEntry(os_data, 'misID SF', 'l')
    legend.Draw('same')

    latexLabel = TLatex()
    latexLabel.SetNDC()
    latexLabel.SetTextFont(62) # helvetica bold face
    latexLabel.SetTextSize(0.4 * c.GetTopMargin())
    latexLabel.DrawLatex(0.11, 0.92, "CMS private work")
    latexLabel.SetTextSize(0.4 * c.GetTopMargin())
    latexLabel.SetTextAlign(31)
    latexLabel.DrawLatex(0.9, 0.92, "%s"%(year))
    latexLabel.Clear()

    #c.Print('ff_plots/sf_' + year + '.png')
    #c.Print('ff_plots/sf_' + year + '.pdf')
    c.Print('ff_plots/test_sf_' + year + '.png')
    c.Print('ff_plots/test_sf_' + year + '.pdf')
    legend.Clear()
    c.Clear()


    ###############
    # store smooth hist
    ###############
    os_data = nevents_region['los']['nsub_data'].Clone('os_data')
    os_data = smoothing(os_data)

    ss_data = nevents_region['tss']['nsub_data'].Clone('ss_data')
    ss_data = smoothing(ss_data)

    lss_data = nevents_region['lss']['nsub_data'].Clone('lss_data')
    lss_data = smoothing(lss_data)

    tos_mc = nevents_region['tos']['nsub_mc'].Clone('tos_mc')
    #tos_mc = smoothing(tos_mc)

    os_data.Multiply(ss_data)
    os_data.Divide(tos_mc)

    sf_file.cd()
    os_data.SetName(year + '_smooth')
    os_data.Write()

    ###############
    # smooth rebinned SF
    ###############
    os_data = nevents_region['los']['nsub_data'].Clone('os_data')
    os_data = smoothing(os_data)
    arr = array.array('d',[0., 40, 120, 200, 280, 400])
    os_data = os_data.Rebin(len(arr)-1, os_data.GetName(), arr)

    ss_data = nevents_region['tss']['nsub_data'].Clone('ss_data')
    ss_data = smoothing(ss_data)
    ss_data = ss_data.Rebin(len(arr)-1, ss_data.GetName(), arr)

    lss_data = nevents_region['lss']['nsub_data'].Clone('lss_data')
    lss_data = smoothing(lss_data)
    lss_data = lss_data.Rebin(len(arr)-1, lss_data.GetName(), arr)
    ss_data.Divide(lss_data)

    tos_mc = nevents_region['tos']['nsub_mc'].Clone('tos_mc')
    tos_mc = tos_mc.Rebin(len(arr)-1, tos_mc.GetName(), arr)
    #tos_mc = smoothing(tos_mc)

    os_data.Multiply(ss_data)
    os_data.Divide(tos_mc)

    for i in range(os_data.GetNbinsX()):
        print("Bin", str(i+1), os_data.GetBinContent(i+1))

    os_data.SetMarkerStyle(1)
    os_data.SetLineWidth(2)
    os_data.SetLineColor(4)
    os_data.SetFillColor(0)
    os_data.Fit('pol1')
    os_data.GetFunction("pol1").SetLineColor(4)
    os_data.SetStats(1)
    gStyle.SetOptStat(000000000)
    gStyle.SetOptFit(11)
    os_data.GetYaxis().SetRangeUser(0, 2)
    os_data.GetYaxis().SetTitle('MisID SF')

    stats1 = os_data.GetListOfFunctions().FindObject("stats").Clone("stats1")
    c.Modified()
    c.Update();
    stats = c.GetPrimitive("stats")
    stats.SetName("Fit stats");
    stats.SetX1NDC(.6);
    stats.SetX2NDC(.9);
    stats.SetY1NDC(.8);
    stats.SetY2NDC(.9);

    os_data.Draw('b e')

    legend.AddEntry(os_data, 'misID SF', 'l')
    legend.Draw('same')

    latexLabel = TLatex()
    latexLabel.SetNDC()
    latexLabel.SetTextFont(62) # helvetica bold face
    latexLabel.SetTextSize(0.4 * c.GetTopMargin())
    latexLabel.DrawLatex(0.11, 0.92, "CMS private work")
    latexLabel.SetTextSize(0.4 * c.GetTopMargin())
    latexLabel.SetTextAlign(31)
    latexLabel.DrawLatex(0.9, 0.92, "%s"%(year))
    latexLabel.Clear()

    #c.Print('ff_plots/sf_' + year + '.png')
    #c.Print('ff_plots/sf_' + year + '.pdf')
    c.Print('ff_plots/test_sf_smooth_' + year + '.png')
    c.Print('ff_plots/test_sf_smooth_' + year + '.pdf')
    legend.Clear()
    c.Clear()

    print('/////////////////////////////////////////////////')

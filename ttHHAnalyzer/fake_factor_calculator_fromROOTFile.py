import os
from uncertainties import ufloat
import ctypes
from ROOT import *

input_ver = 'v9_0714_fake'
isQcd = False

yield_file = 'yields.tex'
if isQcd: yield_file = 'qcd/' + yield_file

year_list = ['2016pre', '2016post', '2017', '2018']
#no_clean = True
no_clean = False

#pt binning: 20, 0, 400
#hist_name = 'h_tau1_pt_S3'
#hist_gen_name = 'h_tau1_gen_pt_S3'
#bin_edges = [[1, -1]]
#bin_edges = [[1, 7], [8, -1]] #(0-140) (140-inf)
#bin_edges = [[1, 3], [4, 5], [6, 7], [8, -1]] #(0-60) (60-100) (100-140) (140-inf)

#dm binning: 1, 2, 11, 12 (dm = bin num - 1)
hist_name = 'h_tau1_decayMode_S3'
hist_gen_name = 'h_tau1_decayMode_gen_S3'
bin_edges = [[1, 1], [2, 2], [11, 11], [12, 12]]

#los = loose tau, OS mu tau (C)
#lss = loose tau, SS mu tau (A)
#tos = tight tau, OS mu tau (D)
#tss = tight tau, SS mu tau (B)


for year in year_list:

    for bin_edge in bin_edges:

        print("Year: ", year)

        nevents_region = {
            'los': {'gentau_total': ufloat(0.0, 0.0), 'mc_total': ufloat(0.0, 0.0), 'ndata': 0.0,
                    'nsub_mc': ufloat(0.0, 0.0), 'nsub_data': ufloat(0.0, 0.0)},
            'lss': {'gentau_total': ufloat(0.0, 0.0), 'mc_total': ufloat(0.0, 0.0), 'ndata': 0.0,
                    'nsub_mc': ufloat(0.0, 0.0), 'nsub_data': ufloat(0.0, 0.0)},
            'tos': {'gentau_total': ufloat(0.0, 0.0), 'mc_total': ufloat(0.0, 0.0), 'ndata': 0.0,
                    'nsub_mc': ufloat(0.0, 0.0), 'nsub_data': ufloat(0.0, 0.0)},
            'tss': {'gentau_total': ufloat(0.0, 0.0), 'mc_total': ufloat(0.0, 0.0), 'ndata': 0.0,
                    'nsub_mc': ufloat(0.0, 0.0), 'nsub_data': ufloat(0.0, 0.0)}
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
            gen_err = ctypes.c_double(0.)
            gen_int = h_gen.IntegralAndError(bin_edge[0], bin_edge[1], gen_err)

            c_all = plot_file.Get(hist_name)
            p_all = c_all.GetPrimitive('pad_hi')
            s_all = p_all.GetPrimitive('mc_stack_0')
            h_all = s_all.GetStack().Last()
            all_err = ctypes.c_double(0.)
            all_int = h_all.IntegralAndError(bin_edge[0], bin_edge[1], all_err)

            prim_all = p_all.GetListOfPrimitives()
            for prim in prim_all:
                if prim.InheritsFrom("TH1"):
                    h_data = prim
                    data_int = h_data.Integral(bin_edge[0], bin_edge[1])

            nevents['gentau_total'] = ufloat(gen_int, gen_err.value)
            nevents['mc_total'] = ufloat(all_int, all_err.value)
            nevents['nsub_mc'] = nevents['mc_total'] - nevents['gentau_total']
            nevents['ndata'] = int(data_int)
            nevents['nsub_data'] = nevents['ndata'] - nevents['gentau_total']

        #print(nevents_region)

        FF = nevents_region['tss']['nsub_data'] / nevents_region['lss']['nsub_data']
        print('FF: {:.5f}'.format(FF.n), 'pm {:.6f}'.format(FF.std_dev), '-> ' + str(round(100*FF.std_dev/FF.n, 3)) + '%')

        fake_nevent = nevents_region['los']['nsub_data'] * FF
        print('Fake event in tos (SR) = los * (tss / lss) = {:.2f}'.format(fake_nevent))

        SF = fake_nevent/nevents_region['tos']['nsub_mc']
        print('SF = new fake / mc fake = {:.4f}'.format(SF.n), 'pm {:5f}'.format(SF.std_dev), '-> ' + str(round(100*SF.std_dev/SF.n, 3)) + '%')
        print('UNC : {0:.4g}'.format(SF.std_dev/SF.n), ', statup {0:.4g},'.format(SF.n+SF.std_dev/SF.n), 'statdown {0:.4g}'.format(SF.n-SF.std_dev/SF.n))

        new_mc = nevents_region['tos']['gentau_total'] + fake_nevent
        print('Recalculated MC = {:.2f}'.format(new_mc))
        print('data/mc = {0:.2g}'.format(nevents_region['tos']['ndata']/nevents_region['tos']['mc_total'].n), '-> {0:.2g}'.format(nevents_region['tos']['ndata']/new_mc.n))

        fake_closure = nevents_region['los']['nsub_mc'] * (nevents_region['tss']['nsub_mc'] / nevents_region['lss']['nsub_mc'])
        print('Fake closure: {:.2f}'.format(fake_closure))

        SF_closure = (nevents_region['tos']['nsub_mc'] - fake_closure)/nevents_region['tos']['nsub_mc']
        print('Difference: {0:.4g} %'.format(100 * (nevents_region['tos']['nsub_mc'].n - fake_closure.n)/nevents_region['tos']['nsub_mc'].n))
        print('UNC : {0:.4g}'.format(SF_closure.n), ', systup {0:.4g},'.format(SF.n+SF_closure.n), 'systdown {0:.4g}'.format(SF.n-SF_closure.n))

        print('        map_ff["{0}"]["nom"]  = {1:.4g};'.format(year, SF.n))
        print('        map_ff["{0}"]["stat"] = {1:.4g};'.format(year, SF.std_dev/SF.n))
        print('        map_ff["{0}"]["syst"] = {1:.4g};'.format(year, SF_closure.n))

        print('/////////////////////////////////////////////////')

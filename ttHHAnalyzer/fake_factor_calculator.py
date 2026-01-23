import os
from uncertainties import ufloat

input_ver = 'v9_0714_fake'
isQcd = False

yield_file = 'yields.tex'
if isQcd: yield_file = 'qcd/' + yield_file

year_list = ['2016pre', '2016post', '2017', '2018']

#los = loose tau, OS mu tau (C)
#lss = loose tau, SS mu tau (A)
#tos = tight tau, OS mu tau (D)
#tss = tight tau, SS mu tau (B)


for year in year_list:

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

        with open(os.path.join(input_ver + '_' + region, 'figure_'+year, yield_file)) as y:
            lines = y.readlines()
            for line in lines:
                if 'Total' in line:
                    list_parse = line.split('&')
                    nevents['gentau_total'] = ufloat(float(list_parse[1][list_parse[1].find('$')+1:list_parse[1].find('{')]), float(list_parse[1][list_parse[1].find('\pm')+4:list_parse[1].find('}$')]))
                    nevents['mc_total'] = ufloat(float(list_parse[2][list_parse[2].find('$')+1:list_parse[2].find('{')]), float(list_parse[2][list_parse[2].find('\pm')+4:list_parse[2].find('}$')]))
                    nevents['nsub_mc'] = nevents['mc_total'] - nevents['gentau_total']
                if 'Data' in line and 'rediction' not in line:
                    list_parse = line.split('&')
                    nevents['ndata'] = int(list_parse[1])
                    nevents['nsub_data'] = nevents['ndata'] - nevents['gentau_total']

    #print(nevents_region)

    FF = nevents_region['tss']['nsub_data'] / nevents_region['lss']['nsub_data']
    print('FF: {:.5f}'.format(FF.n), 'pm {:.6f}'.format(FF.std_dev), '-> ' + str(round(100*FF.std_dev/FF.n, 3)) + '%')

    fake_nevent = nevents_region['los']['nsub_data'] * FF
    print('Fake event in tos (SR) = los * (tss / lss) = {:.2f}'.format(fake_nevent))

    SF = fake_nevent/nevents_region['tos']['nsub_mc']
    print('SF = new fake / mc fake = {:.4f}'.format(SF.n), 'pm {:5f}'.format(SF.std_dev), '-> ' + str(round(100*SF.std_dev/SF.n, 3)) + '%')
    print('UNC : {0:.4g} '.format(SF.std_dev/SF.n), 'statup {0:.4g},'.format(1+SF.std_dev/SF.n), 'statdown {0:.4g}'.format(1-SF.std_dev/SF.n))

    new_mc = nevents_region['tos']['gentau_total'] + fake_nevent
    print('Recalculated MC = {:.2f}'.format(new_mc), ' from {:.2f}'.format(nevents_region['tos']['mc_total']), ' ratio: {:.2f}'.format(new_mc/nevents_region['tos']['mc_total']))
    print('data/mc = {0:.2g}'.format(nevents_region['tos']['ndata']/nevents_region['tos']['mc_total'].n), '-> {0:.2g}'.format(nevents_region['tos']['ndata']/new_mc.n))

    fake_closure = nevents_region['los']['nsub_mc'] * (nevents_region['tss']['nsub_mc'] / nevents_region['lss']['nsub_mc'])
    print('Fake closure: {:.2f}'.format(fake_closure))

    SF_closure = (nevents_region['tos']['nsub_mc'] - fake_closure)/nevents_region['tos']['nsub_mc']
    print('Difference: {0:.4g} %'.format(100 * (nevents_region['tos']['nsub_mc'] - fake_closure)/nevents_region['tos']['nsub_mc']))
    print('UNC : {0:.4g} '.format(SF_closure.n), 'systup {0:.4g},'.format(1+SF_closure.n), 'systdown {0:.4g}'.format(1-SF_closure.n))

    print('/////////////////////////////////////////////////')

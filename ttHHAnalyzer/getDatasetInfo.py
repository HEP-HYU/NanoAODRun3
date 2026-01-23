import sys, os, socket, csv
from urllib.request import urlopen
import pandas as pd

if len(sys.argv) < 2:
    print("Provide run era: 2016pre, 2016post, 2017, 2018, 2022, 2023_BPix, 2024")
    sys.exit()
else: era = sys.argv[1]

gdocbase = "https://docs.google.com/spreadsheets/d/1A2cKv080xvfVEJBhkDo6NzQdvXbeGu95/pub?gid=%s&single=true&output=csv"

def getCampaignSummary():
    url = gdocbase % "545547353"
    print("Retrieving campaign info...")
    print("Source URL = ", url)
    df = pd.read_csv(urlopen(url)).set_index('Campaign')
    return df

url, era = '', ''
if len(sys.argv) < 2:
    print("Missing year tag")
    sys.exit()
else:
    gidMap = getCampaignSummary()
    era = sys.argv[1]
    eraData = gidMap.loc[era]
    print("="*40)
    print(eraData)
    print("="*40)
    url = gdocbase % int(eraData["GID"])
    do_all = False
    if len(sys.argv) > 2: do_all = sys.argv[2] == '-f'

if url != '':
    print("Retrieving dataset info...")
    print("Source URL =", url)
    content = list(csv.reader(urlopen(url).read().decode('utf8').splitlines()))

def str2float(x):
    x = x.strip()
    if x == '': return 0.0
    return float(x)
def str2int(x):
    x = x.strip()
    if x == '': return 0
    return int(x)

ds = []
titleBabel = {
    "Name":"name","DataSetName":"DataSetName",
    "Cross section (pb)":"xsec", "Raw events":"nevt",
    "Path":"path", 
}
titleRule = {
    "xsec":str2float, "nevt":str2int,
}

## analyze CSV file
titleIdx = {}
for titleCSV in titleBabel:
    titleJSON = titleBabel[titleCSV]
    if titleCSV not in content[0]: continue
    titleIdx[titleJSON] = content[0].index(titleCSV)

for l in content[1:]:
    if len(l) == 0 or len(l[0]) == 0: continue

    item = {}
    for title in titleIdx:
        val = l[titleIdx[title]]
        if title not in titleRule: item[title] = val
        else: item[title] = titleRule[title](val)
    ds.append(item)

## Save results
outDir = os.path.join(os.getcwd(), "data/dataset", era)
if not os.path.exists(outDir):
    os.makedirs(outDir)

import json
f = open("%s/dataset.json" % outDir, "w")
print("Saving dataset info to %s/dataset.json" % outDir)
print(json.dumps(ds, indent=4, sort_keys=True), file=f)
f.close()


exit()
era_path = {'2016pre': ()}

from ROOT import *


hostname = os.environ["HOSTNAME"]
if any(i in hostname for i in ['htop', 'compute-0', 'gpu-0']):
    storageprefix = "/data1/common/NanoAOD"
elif 'sdfarm.kr' in hostname:
    storageprefix = "root://cms-xrdr.sdfarm.kr:1094//store"
elif "knu.ac.kr" in hostname:
    storageprefix = "root://cluster142.knu.ac.kr:1094///store/"
    if os.system("voms-proxy-info -exists --valid 8:00") != 0:
        os.system("voms-proxy-init -voms cms --valid 144:00")
else:
    print("Unknown server, cannot create file list!")
    sys.exit()

emptydirectory = []
print("Listing dataset location")
print("-"*40)
for d in ds:
    num_lines = 0
    if os.path.isfile("%s/dataset_%s.txt" % (outDir, d['name'])):
        num_lines = sum(1 for line in open("%s/dataset_%s.txt" % (outDir, d['name']), "r"))
    if num_lines > 5 and do_all == False: continue

    dataSetName = d['DataSetName']
    dataSetName = dataSetName.split('/')
    if len(dataSetName) < 4: continue

    pdName, sdName = dataSetName[1], dataSetName[2]
    if sdName.startswith("Run20"):
        runEra = sdName.split("-")[0]
        path = "%s/data/%s/%s/NANOAOD/%s" % (storageprefix, runEra, pdName, sdName.replace(runEra + '-', ''))
    else:
        mcCampaign = sdName.split("-")[0]
        path = "%s/mc/%s/%s/NANOAODSIM/%s" % (storageprefix, mcCampaign, pdName, sdName.replace(mcCampaign + '-', ''))

    d['files'] = []
    if os.path.exists(path):
        for d_ in os.listdir(path):
            for f_ in os.listdir(os.path.join(path, d_)):
                d['files'].append(os.path.join(path, d_, f_))
    d['files'].sort()

    f = open("%s/dataset_%s.txt" % (outDir, d['name']), "w")
    #for key in d:
    #    if key == "files": continue
    #    print("# " + key + " = " + str(d[key]), end="\n", file=f)
    for ff in d['files']:
        print(os.path.join(d['path'], ff), end="\n", file=f)
    f.close()
    if len(d['files']) == 0:
        emptydirectory.append(d['name'])
print("-"*40)


if len(emptydirectory) > 0:
    print('List of datasets with no data')
    for i in emptydirectory:
        print(i)
    print("-"*40)

print("Saved dataset location info to %s/dataset_SampleName.txt" % (outDir))

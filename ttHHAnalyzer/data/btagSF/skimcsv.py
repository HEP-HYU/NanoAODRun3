import pandas as pd
import json

pd.set_option("mode.chained_assignment",  None)
infiles = {"2016post": "reshaping_deepJet_2016post.csv",
           "2016pre" : "reshaping_deepJet_2016pre.csv",
           "2017"    : "reshaping_deepJet_2017.csv",
           "2018"    : "reshaping_deepJet_2018.csv"}
#files from btv gitlab, v2 SFs

btag_var = ["hf","lf","hfstats1","hfstats2","lfstats1","lfstats2","cferr1","cferr2"]
#jes_var = ["Absolute", "Absolute_year", "BBEC1", "BBEC1_year", "EC2", "EC2_year", "FlavorQCD", "HF", "HF_year", "RelativeBal", "RelativeSample_year"]
#EC1: 1.3<|eta|<2.5; EC2: 2.5<|eta|<3.0; HF: |eta|>3
jes_var = ["Absolute", "Absolute_year", "BBEC1", "BBEC1_year", "FlavorQCD", "RelativeBal", "RelativeSample_year"]

for year, infile in infiles.items():
    print(infile)
    df = pd.read_csv(infile)

    btag_systs = ["central"]
    for var in btag_var:
        btag_systs.append("up_" + var)
        btag_systs.append("down_" + var)
    newdf = df.loc[df["sysType"].isin(btag_systs)]
    newdf.loc[df["OperatingPoint"] == "shape", "OperatingPoint"] = 3

    out_btag_systs = ["central"]
    for syst in btag_systs:
        if syst == "central": continue
        if "up" in syst:
            newdf["sysType"] = newdf["sysType"].replace(syst, syst.lstrip("up_") + "up")
            out_btag_systs.append(syst.lstrip("up_") + "up")
        elif "down" in syst:
            newdf["sysType"] = newdf["sysType"].replace(syst, syst.lstrip("down_") + "down")
            out_btag_systs.append(syst.lstrip("down_") + "down")
    newdf.to_csv("skimmed_btag_" + year + ".csv", index=False)
    print(json.dumps(out_btag_systs))


    jes_systs   = ["central"]
    for var in jes_var:
        if "year" in var:
            jes_systs.append("up_jes" + var.replace("year", year[:4]))
            jes_systs.append("down_jes" + var.replace("year", year[:4]))
        else:
            jes_systs.append("up_jes" + var)
            jes_systs.append("down_jes" + var)
    #print(json.dumps(jes_systs))
    newjesdf = df.loc[df["sysType"].isin(jes_systs)]
    newjesdf.loc[df["OperatingPoint"] == "shape","OperatingPoint"] = 3

    out_jes_systs = ["central"]
    for syst in jes_systs:
        if "up" in syst:
            newjesdf["sysType"] = newjesdf["sysType"].replace(syst, (syst.lstrip("up_") + "up"))
            out_jes_systs.append((syst.lstrip("up_") + "up"))
        elif "down" in syst:
            newjesdf["sysType"] = newjesdf["sysType"].replace(syst, (syst.lstrip("down_") + "down"))
            out_jes_systs.append((syst.lstrip("down_") + "down"))
    newjesdf.to_csv("skimmed_jes_" + year + ".csv", index=False)
    print(json.dumps(out_jes_systs))

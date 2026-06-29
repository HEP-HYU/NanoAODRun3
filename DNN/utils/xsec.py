# Cross Sections (pb)
# Canonical cross-section reference for the LFV analysis.
# Values must be kept consistent with plotIt YAML configurations.
# Units: pb
#
# NOTE: Previously this dict had duplicate keys for LFV_TT_TToCMuTau_*
# (identical values, last definition was silently winning).
# Duplicates removed 2026-06: values are unchanged.
xsec = {
    "DYJetsToLL_M-10to50"             : 18610.00,
    "DYJetsToLL_M-50-amcatnloFXFX"    : 6077.22,
    "DYJetsToLL_M-50-madgraphMLM"     : 6077.22,
    "WJetsToLNu_madgraphMLM"          : 61526.7,
    "WJetsToLNu_inclHT100"            : 59373.3,
    "WJetsToLNu_HT-100To200"          : 1345   * 0.993,
    "WJetsToLNu_HT-200To400"          : 359.7  * 1.002,
    "WJetsToLNu_HT-400To600"          : 48.91  * 1.009,
    "WJetsToLNu_HT-600To800"          : 12.05  * 1.120,
    "WJetsToLNu_HT-800To1200"         : 5.501  * 1.202,
    "WJetsToLNu_HT-1200To2500"        : 1.329  * 1.332,
    "WJetsToLNu_HT-2500ToInf"         : 0.03216 * 4.200,
    "ST_t-channel_top"                : 136.02,
    "ST_t-channel_antitop"            : 80.95,
    "ST_tW_top"                       : 35.85,
    "ST_tW_antitop"                   : 35.85,
    "TTTo2L2Nu"                       : 88.29,
    "TTToHadronic"                    : 377.96,
    "TTToSemiLeptonic"                : 365.34,
    # LFV TT signals (tc coupling)
    "LFV_TT_TToCMuTau_Scalar"         : 2.688 * 1E-03,
    "LFV_TT_TToCMuTau_Vector"         : 21.50 * 1E-03,
    "LFV_TT_TToCMuTau_Tensor"         : 129.0 * 1E-03,
    # LFV TT signals (tu coupling)
    "LFV_TT_TToUMuTau_Scalar"         : 2.688 * 1E-03,
    "LFV_TT_TToUMuTau_Vector"         : 21.50 * 1E-03,
    "LFV_TT_TToUMuTau_Tensor"         : 129.0 * 1E-03,
    # LFV ST signals (tc coupling)
    "LFV_ST_TCMuTau_Scalar"           : 7.40  * 1E-03,
    "LFV_ST_TCMuTau_Vector"           : 36.8  * 1E-03,
    "LFV_ST_TCMuTau_Tensor"           : 178.4 * 1E-03,
    # LFV ST signals (tu coupling)
    "LFV_ST_TUMuTau_Scalar"           : 83.8  * 1E-03,
    "LFV_ST_TUMuTau_Vector"           : 393   * 1E-03,
    "LFV_ST_TUMuTau_Tensor"           : 1796  * 1E-03,
    "LQ_signal"                       : 0.003,
    "WWTo2L2Nu"                       : 12.178,
    "WWToLNuQQ"                       : 49.997,
    "WZTo2L2Q"                        : 5.595,
    "WZTo3LNu"                        : 4.42965,
    "ZZTo2L2Q"                        : 3.220,
    "TTWJetsToLNu"                    : 0.2043,
    "TTWJetsToQQ"                     : 0.4062,
    "TTZToLLNuNu_M-10"               : 0.2529,
    "TTZToQQ"                         : 0.5297,
    # QCD (muon-enriched)
    "QCD_Pt-15to20_MuEnrichedPt5"     : 1273190000 * 0.003,
    "QCD_Pt-20to30_MuEnrichedPt5"     : 558528000  * 0.0053,
    "QCD_Pt-30to50_MuEnrichedPt5"     : 139803000  * 0.01182,
    "QCD_Pt-50to80_MuEnrichedPt5"     : 19222500   * 0.02276,
    "QCD_Pt-80to120_MuEnrichedPt5"    : 2758420    * 0.03844,
    "QCD_Pt-120to170_MuEnrichedPt5"   : 469797     * 0.05362,
    "QCD_Pt-170to300_MuEnrichedPt5"   : 117989     * 0.07335,
    "QCD_Pt-300to470_MuEnrichedPt5"   : 7820.25    * 0.10196,
    "QCD_Pt-470to600_MuEnrichedPt5"   : 645.528    * 0.12242,
    "QCD_Pt-600to800_MuEnrichedPt5"   : 187.109    * 0.13412,
    "QCD_Pt-800to1000_MuEnrichedPt5"  : 32.3486    * 0.14552,
    "QCD_Pt-1000toInf_MuEnrichedPt5"  : 10.4305    * 0.15544,
}  # pb

import os, sys
import numpy as np
import statsmodels.api as sm
from statsmodels.nonparametric.kde import KDEUnivariate
import ROOT
from ROOT import *


def smoothing(hin):

    htmp = hin.Clone()

    x_vals = np.arange(htmp.GetNbinsX())
    y_vals = np.zeros(htmp.GetNbinsX())

    for i in range(htmp.GetNbinsX()):
        y_vals[i] = htmp.GetBinContent(i+1)

    #local linear regression (locally weighted polynomial regression)
    lowess = sm.nonparametric.lowess
    smoothed_vals = np.zeros(y_vals.shape)
    smoothed_vals = lowess(y_vals, x_vals, frac=1./3, return_sorted=False)


    for x_position in range(htmp.GetNbinsX()):
        hin.SetBinContent(x_position+1, max(0, smoothed_vals[x_position]))

    return hin




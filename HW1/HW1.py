# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 16:57:15 2019

@author: WHITE
"""
#%% IMPORT MODULES
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import stats
#%% MAKE RANDOM NUMBER
n = 1000
np.random.seed(n)
x = np.random.normal(0, 1, n)
#%% MAKE FUNCTIONS
def STAT(x):
    mean = np.nanmean(x)
    var  = np.nanvar(x)
    skew = stats.moment(x, 3) / (np.nanstd(x) ** 3)
    kurt = stats.moment(x, 4) / (np.nanstd(x) ** 4)
    print("Mean: %.4f\nVariance: %.4f\nSkewness: %.4f\nKurtosis: %.4f" % (mean, var, skew, kurt))
#%% SOLVE
STAT(x)

hist, bins = np.histogram(x, bins = np.linspace(-3., 3., 16), density = False)
hist_norm  = hist / n
mpl.rc("font", weight = "bold")
plt.bar([(bins[i] + bins[i + 1]) / 2 for i in range(len(bins) - 1)], hist_norm, color = "grey", edgecolor = "k", width = bins[-1] - bins[-2])
plt.xlim(-3, 3)
plt.ylim(0)
plt.xlabel("X", fontdict = {"weight": "bold"})
plt.ylabel("Probability", fontdict = {"weight": "bold"})
plt.title("n = " + str(n), fontdict = {"weight": "bold"})
#%% HW1-7
n = 101
x = np.arange(n)
y = np.sin(2 * np.pi * x / 10)

STAT(y)

hist, bins = np.histogram(y, bins = np.linspace(-1., 1., 16), density = False)
hist_norm  = hist / n
mpl.rc("font", weight = "bold")
plt.bar([(bins[i] + bins[i + 1]) / 2 for i in range(len(bins) - 1)], hist_norm, color = "grey", edgecolor = "k", width = bins[-1] - bins[-2])
plt.xlim(-1, 1)
plt.ylim(0)
plt.xlabel("Y", fontdict = {"weight": "bold"})
plt.ylabel("Probability", fontdict = {"weight": "bold"})
plt.title("n = " + str(n), fontdict = {"weight": "bold"})
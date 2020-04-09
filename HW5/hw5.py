# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 12:49:52 2019

@author: WHITE
"""
#%% IMPORT MODULES
import numpy as np
import matplotlib as mpl
mpl.rc("font", weight = "bold")
import matplotlib.pyplot as plt
from statsmodels.tsa.stattools import acovf
from scipy.signal import periodogram
#%% FUNCTIONS
def FUNC1(i, t):
    return np.sin(2 * np.pi * i * t / 10)
def FUNC2(i, t):
    return np.sin(2 * np.pi * i * t / 10) + np.sin(2 * np.pi * i * t / 5) / 2
def WHITE(n):
    return np.random.normal(0, 1, n)
def RED(n, a):
    result = np.zeros(n)
    error  = np.random.normal(0, 1, n)
    for ni in range(n - 1):
        result[ni + 1] = a * result[ni] + ((1 - a ** 2) ** 0.5) * error[ni + 1]
    return result
#%% DRAW FUNCTION
def DRAWFUNC(time, x, tlim):
    fig, sub = plt.subplots()
    sub.set_title("Original function", fontdict = {"weight": "bold"})
    sub.plot(time, x, ls = "-", c = "k")
    sub.set_xlabel("Time", fontdict = {"weight": "bold"})
    sub.set_xlim(min(time), tlim)
def DRAWAUTOCOV(x, llim):
    ac  = acovf(x, unbiased = True, demean = True)
    lag = np.arange(len(ac))
    fig, sub = plt.subplots()
    sub.set_title("Autocovariance", fontdict = {"weight": "bold"})
    sub.plot(lag, ac, ls = "-", c = "k")
    sub.set_xlabel("Time lag", fontdict = {"weight": "bold"})
    sub.set_xlim(min(lag), llim)
def DRAWPS(x, sampling_frequency, window, separate, scaling, xscale, yscale, xunit, yunit):
    for ss in range(separate):
        xlen = int(len(x) / separate)
        if ss == 0:
            freq, ps = periodogram(x[xlen * ss : xlen * (ss + 1)], fs = sampling_frequency, window = window, detrend = "linear", scaling = scaling)
        else:
            _, ps_ss = periodogram(x[xlen * ss : xlen * (ss + 1)], fs = sampling_frequency, window = window, detrend = "linear", scaling = scaling)
            ps = ps + ps_ss
    ps = ps / separate
    fig, sub = plt.subplots()
    sub.set_title("Power spectrum (" + window + ")", fontdict = {"weight": "bold"})
    freq = freq[1 :]
    ps   = ps[1 :]
    xlabel = "Frequency"
    if xscale == "log":
        freq   = np.log10(freq)
        xlabel = "log(" + xlabel + ")"
    xlabel = xlabel + " [" + str(xunit) + "]"
    ylabel = "Phi"
    if yscale == "log":
        ps     = np.log10(ps)
        ylabel = "log(" + ylabel + ")"
    if xscale == "log" and yscale == "log":
        nm = separate * 2
        if window == "hann":
            dof = 8 * nm / 3
        elif window == "parzen":
            dof = 3.7086 * nm
        print("D.O.F:", dof)
        chi1 = input("Chi square 0.025:")
        chi2 = input("Chi square 0.975:")
        conf1 = np.log10(dof / float(chi1))
        conf2 = np.log10(dof / float(chi2))
        sub.fill_between(freq, ps + conf1, ps + conf2, facecolor = "r", alpha = 0.8, interpolate = True)
        sub.text(min(freq), max(ps), "D.O.F.=%.1f" % dof, ha = "left", va = "top")
    if scaling == "density":
        yunit = yunit + "/" + xunit
    ylabel = ylabel + " [" + str(yunit) + "]"
    sub.plot(freq, ps, ls = "-", lw = 0.5, c = "k")
    sub.set_xlabel(xlabel, fontdict = {"weight": "bold"})
    sub.set_ylabel(ylabel, fontdict = {"weight": "bold"})
    sub.set_xlim(min(freq), max(freq))
    first_idx = np.argmax(ps)
    sub.text(freq[first_idx], ps[first_idx], "%.1f" %(ps[first_idx]), ha = "center", va = "bottom")
#%% PROBLEM#1
i = np.arange(1000 + 1)
t = 1
x = FUNC1(i = i, t = t)
DRAWFUNC(time = i, x = x, tlim = 50)
DRAWAUTOCOV(x = x, llim = 50)
DRAWPS(x = x, flim = 0.5, sampling_frequency = 1, window = "bartlett", nfft = None, scaling = "density")
#%% PROBLEM#2
i = np.arange(1000 + 1)
t = 1
x = FUNC2(i = i, t = t)
DRAWFUNC(time = i, x = x, tlim = 50)
DRAWAUTOCOV(x = x, llim = 50)
DRAWPS(x = x, flim = 0.5, sampling_frequency = 1, window = "bartlett", nfft = None, scaling = "density")
#%% PROBLEM#3
i = np.arange(1000)
x = WHITE(1000)
DRAWFUNC(time = i, x = x, tlim = 50)
DRAWAUTOCOV(x = x, llim = 50)
DRAWPS(x = x, flim = 0.5, sampling_frequency = 1, window = "bartlett", nfft = None, scaling = "density")
#%% PROBLEM#4
i = np.arange(1000)
x = RED(n = 1000, a = 0.5)
DRAWFUNC(time = i, x = x, tlim = 50)
DRAWAUTOCOV(x = x, llim = 50)
DRAWPS(x = x, flim = 0.5, sampling_frequency = 1, window = "bartlett", nfft = None, scaling = "density")
#%% PROBLEM#5
data = np.loadtxt("C:/Users/WHITE/Google 드라이브/SNU/2019.1/객관적자료분석/HW5/pus9802.dat")
juli = data[:, 0]
tide = data[:, 1]
DRAWFUNC(time = juli, x = tide, tlim = max(juli))
#%% PROBLEM#5-1
separate = 5
DRAWPS(x = tide, sampling_frequency = 1, window = "hann",\
       separate = separate, scaling = "density", xscale = "log", yscale = "log", xunit = "cph", yunit = r"cm$^2$")
DRAWPS(x = tide, sampling_frequency = 1, window = "hann",\
       separate = separate, scaling = "spectrum", xscale = "log", yscale = None, xunit = "cph", yunit = r"cm$^2$")
DRAWPS(x = tide, sampling_frequency = 1, window = "parzen",\
       separate = separate, scaling = "density", xscale = "log", yscale = "log", xunit = "cph", yunit = r"cm$^2$")
DRAWPS(x = tide, sampling_frequency = 1, window = "parzen",\
       separate = separate, scaling = "spectrum", xscale = "log", yscale = None, xunit = "cph", yunit = r"cm$^2$")
#%% PROBLEM#5-2
separate = 1
DRAWPS(x = tide, sampling_frequency = 1, window = "hann",\
       separate = separate, scaling = "density", xscale = "log", yscale = "log", xunit = "cph", yunit = r"cm$^2$")
DRAWPS(x = tide, sampling_frequency = 1, window = "hann",\
       separate = separate, scaling = "spectrum", xscale = "log", yscale = None, xunit = "cph", yunit = r"cm$^2$")
DRAWPS(x = tide, sampling_frequency = 1, window = "parzen",\
       separate = separate, scaling = "density", xscale = "log", yscale = "log", xunit = "cph", yunit = r"cm$^2$")
DRAWPS(x = tide, sampling_frequency = 1, window = "parzen",\
       separate = separate, scaling = "spectrum", xscale = "log", yscale = None, xunit = "cph", yunit = r"cm$^2$")
#%% PROBLEM#5-3
separate = 5
DRAWPS(x = tide[::12], sampling_frequency = 1, window = "hann",\
       separate = separate, scaling = "density", xscale = "log", yscale = "log", xunit = "cph", yunit = r"cm$^2$")
DRAWPS(x = tide[::12], sampling_frequency = 1, window = "hann",\
       separate = separate, scaling = "spectrum", xscale = "log", yscale = None, xunit = "cph", yunit = r"cm$^2$")
DRAWPS(x = tide[::12], sampling_frequency = 1, window = "parzen",\
       separate = separate, scaling = "density", xscale = "log", yscale = "log", xunit = "cph", yunit = r"cm$^2$")
DRAWPS(x = tide[::12], sampling_frequency = 1, window = "parzen",\
       separate = separate, scaling = "spectrum", xscale = "log", yscale = None, xunit = "cph", yunit = r"cm$^2$")
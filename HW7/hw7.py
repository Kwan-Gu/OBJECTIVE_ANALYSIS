# -*- coding: utf-8 -*-
"""
Created on Fri May 10 12:51:26 2019

@author: WHITE
"""
#%% IMPORT MODULES
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import periodogram
#%% RAW DATA
raw_file = np.loadtxt("C:/Users/WHITE/Google 드라이브/SNU/2019.1/객관적자료분석/HW7/pus9802.dat")
day = raw_file[:, 0]
raw = raw_file[:, 1]

fig, sub = plt.subplots(figsize = (18, 4))
sub.plot(day, raw, marker = None, ls = "-", c = "k", lw = 0.5)
sub.set_xlim(min(day), max(day))
sub.set_ylim(-50, 200)
sub.set_xlabel("Day", fontdict = {"weight": "bold"})
sub.set_ylabel("Sea level [cm]", fontdict = {"weight": "bold"})
sub.set_title("Original", fontdict = {"weight": "bold"})
plt.savefig("./original.png", dpi = 500, bbox_inches = "tight")
#%% LOW PASS FILTER
low_file = np.loadtxt("C:/Users/WHITE/Google 드라이브/SNU/2019.1/객관적자료분석/HW7/pus_low.dat")
low_day = low_file[:, 0]
low     = low_file[:, 1]

fig, sub = plt.subplots(figsize = (18, 4))
sub.plot(low_day, low, marker = None, ls = "-", c = "r", lw = 0.5)
sub.set_xlim(min(day), max(day))
sub.set_ylim(-50, 200)
sub.set_xlabel("Day", fontdict = {"weight": "bold"})
sub.set_ylabel("Sea level [cm]", fontdict = {"weight": "bold"})
sub.set_title("Low pass filter", fontdict = {"weight": "bold"})
plt.savefig("./low.png", dpi = 500, bbox_inches = "tight")
#%% HIGH PASS FILTER
high_file = np.loadtxt("C:/Users/WHITE/Google 드라이브/SNU/2019.1/객관적자료분석/HW7/pus_high.dat")
high_day  = high_file[:, 0]
high      = high_file[:, 1]

fig, sub = plt.subplots(figsize = (18, 4))
sub.plot(high_day, high, marker = None, ls = "-", c = "b", lw = 0.5)
sub.set_xlim(min(day), max(day))
sub.set_ylim(-50, 200)
sub.set_xlabel("Day", fontdict = {"weight": "bold"})
sub.set_ylabel("Sea level [cm]", fontdict = {"weight": "bold"})
sub.set_title("High pass filter", fontdict = {"weight": "bold"})
plt.savefig("./high.png", dpi = 500, bbox_inches = "tight")
#%% COMPARE
com1 = raw[54 : -54] - (low - np.mean(low))
com2 = high
com  = com1 - com2

fig, sub = plt.subplots(figsize = (18, 4))
sub.plot(high_day, com, marker = None, ls = "-", c = "k", lw = 1)
sub.set_xlim(min(day), max(day))
sub.set_ylim(-50, 200)
sub.set_xlabel("Day", fontdict = {"weight": "bold"})
sub.set_ylabel("Sea level [cm]", fontdict = {"weight": "bold"})
sub.set_title("Compare", fontdict = {"weight": "bold"})
plt.savefig("./compare.png", dpi = 500, bbox_inches = "tight")
#%% BAND PASS FILTER
band_file = np.loadtxt("C:/Users/WHITE/Google 드라이브/SNU/2019.1/객관적자료분석/HW7/pus_band.dat")
band_day  = band_file[:, 0]
band      = band_file[:, 1]

fig, sub = plt.subplots(figsize = (18, 4))
sub.plot(band_day, band, marker = None, ls = "-", c = "g", lw = 0.5)
sub.set_xlim(min(day), max(day))
sub.set_ylim(-50, 200)
sub.set_xlabel("Day", fontdict = {"weight": "bold"})
sub.set_ylabel("Sea level [cm]", fontdict = {"weight": "bold"})
sub.set_title("Band pass filter", fontdict = {"weight": "bold"})
plt.savefig("./band.png", dpi = 500, bbox_inches = "tight")
#%% POWER SPECTRUM
fr_raw, ps_raw = periodogram(raw, fs = 1, window = "parzen", scaling = "density")
fr_low, ps_low = periodogram(low, fs = 1, window = "parzen", scaling = "density")
fr_high, ps_high = periodogram(high, fs = 1, window = "parzen", scaling = "density")
fr_band, ps_band = periodogram(band, fs = 1, window = "parzen", scaling = "density")

fig, sub = plt.subplots(figsize = (7, 5))
sub.plot(fr_raw[1 :], ps_raw[1 :], c = "k", ls = "-", lw = 0.5, marker = None, label = "Original")
sub.plot(fr_low[1 :], ps_low[1 :], c = "r", ls = "-", lw=  0.5, marker = None, label = "Low filter")
sub.plot(fr_high[1 :], ps_high[1 :], c = "b", ls = "-", lw=  0.5, marker = None, label = "High filter")
sub.plot(fr_band[1 :], ps_band[1 :], c = "g", ls = "-", lw = 0.5, marker = None, label = "Band filter")
sub.set_xscale("log")
sub.set_yscale("log")
sub.set_xlim(left = 0.0001, right = 0.5)
sub.set_ylim(10**(-14), 10**9)
sub.set_xlabel("Frequency [cph]", fontdict = {"weight": "bold"})
sub.set_ylabel(r"Phi [cm$^2$/cph]", fontdict = {"weight": "bold"})
sub.legend(prop = {"weight": "bold"})
plt.savefig("./ps.png", dpi = 500, bbox_inches = "tight")
# -*- coding: utf-8 -*-
"""
Created on Fri May  3 13:52:55 2019

@author: WHITE
"""
#%% IMPORT MODULES
import numpy as np
import matplotlib.pyplot as plt
#%% READ DATA
sok = np.loadtxt("./sok_result.txt", skiprows = 1)
fr_sok = sok[:, 1]
pr_sok = sok[:, 2]
ch_sok = sok[:, 3]
ph_sok = sok[:, 4]
for i in range(len(ph_sok)):
    if ph_sok[i] < 0:
        ph_sok[i] = ph_sok[i] + 360
gn_sok = sok[:, 5]

nish = np.loadtxt("./nish_result.txt", skiprows = 1)
fr_nish = nish[:, 1]
pr_nish = nish[:, 2]
ch_nish = nish[:, 3]
ph_nish = nish[:, 4]
for i in range(len(ph_nish)):
    if ph_nish[i] < 0:
        ph_nish[i] = ph_nish[i] + 360
gn_nish = nish[:, 5]

if np.sum(fr_sok == fr_nish) == len(fr_sok):
    freq = fr_sok * 24
    xlim = 0.5 # cpd
    fr_idx = np.where(freq < xlim)[0][-1]
else:
    raise IndexError("len(fr_sok) and len(fr_nish) are different!!!")
#%% DEGREE OF FREEDOM
dof  = 37.09149
clev = 1 - np.power(0.05, 1 / (dof - 1))
#%% PLOT COHERENCE
fig, sub = plt.subplots()
sub.plot(freq[: fr_idx], ch_sok[: fr_idx], c = "k", ls = "-", lw = 1, label = "Sokcho")
sub.plot(freq[: fr_idx], ch_nish[: fr_idx], c = "grey", ls = "-", lw = 1, label = "Nishinoomote")
sub.axhline(clev, c = "k", ls = "--", label = "95% confidence level", lw = 1)
sub.set_xscale("log")
sub.set_xlim(right = xlim)
sub.set_ylim(0, 1)
sub.set_xlabel("Frequency (cpd)", fontdict = {"weight": "bold"})
sub.set_ylabel("Coherence", fontdict = {"weight": "bold"})
sub.legend(loc = "best", prop = {"weight": "bold"})
#%% PLOT GAIN
fig, sub = plt.subplots()
sub.plot(freq[: fr_idx], gn_sok[: fr_idx], c = "k", ls = "-", lw = 1, label = "Sokcho")
sub.plot(freq[: fr_idx], gn_nish[: fr_idx], c = "grey", ls = "-", lw = 1, label = "Nishinoomote")
sub.set_xscale("log")
sub.set_xlim(right = xlim)
sub.set_xlabel("Frequency (cpd)", fontdict = {"weight": "bold"})
sub.set_ylabel("Gain [cm/mbar]", fontdict = {"weight": "bold"})
sub.legend(loc = "best", prop = {"weight": "bold"})
#%% PLOT PHASE RELATION
fig, sub = plt.subplots()
sub.plot(freq[: fr_idx], ph_sok[: fr_idx], c = "k", ls = "-", lw = 1, label = "Sokcho")
sub.plot(freq[: fr_idx], ph_nish[: fr_idx], c = "grey", ls = "-", lw = 1, label = "Nishinoomote")
sub.set_xscale("log")
sub.set_xlim(right = xlim)
sub.set_xlabel("Frequency (cpd)", fontdict = {"weight": "bold"})
sub.set_ylabel("Phase [degree]", fontdict = {"weight": "bold"})
sub.legend(loc = "best", prop = {"weight": "bold"})
# -*- coding: utf-8 -*-
"""
Created on Thu May 16 14:24:20 2019

@author: WHITE
"""
#%% IMPORT MODULES
import numpy as np
import matplotlib.pyplot as plt
import pycwt as wavelet
#%% HW8
dt  = 1
i   = np.arange(0, 200 + dt, dt)
eta = np.zeros(len(i), dtype = float)

# HW8-1

for ii in i:
    if ii >= 0 and ii <= 100:
        eta[ii] = np.sin(2 * np.pi * ii / 20)
    else:
        eta[ii] = np.sin(2 * np.pi * ii / 10) / 2

# HW8-2
"""
for ii in i:
    eta[ii] = np.sin(2 * np.pi * ii / 20) + np.sin(2 * np.pi * ii / 10) / 2
"""
fig, sub = plt.subplots(figsize = (10, 4))
sub.plot(i, eta, ls = "-", c = "k", lw = 1)
sub.set_xlim(min(i), max(i))
sub.set_xticks(np.arange(0, 200 + 20, 20))
sub.set_xlabel("i", fontdict = {"weight": "bold"})
sub.set_ylabel("Eta", fontdict = {"weight": "bold"})
sub.set_title("Data", fontdict = {"weight": "bold"})

mother = wavelet.Morlet(f0 = 6)
alpha, _, _ = wavelet.ar1(eta)
dj = 0.25
s0 = 2 * dt
J  = 7 / dj

wave, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(signal = eta, dt = dt,\
                                                      dj = dj, s0 = s0, J = J,\
                                                      wavelet = mother)

power     = np.abs(wave) ** 2
fft_power = np.abs(fft) ** 2
period    = 1 / freqs

signif, fft_theor = wavelet.significance(signal = 1.0, dt = dt,\
                                         scales = scales, alpha = alpha,\
                                         significance_level = 0.95,\
                                         wavelet = mother)
sig95 = np.ones([1, len(eta)]) * signif[:, None]
sig95 = power / sig95

phase = np.zeros(np.shape(wave), dtype = float)
for jj in range(len(period)):
    for ii in range(len(i)):
        real = wave[jj, ii].real
        imag = wave[jj, ii].imag
        phase[jj, ii] = np.arctan(imag / real)

fig, sub = plt.subplots(figsize = (10, 4))
sub1 = sub.twinx()
contf = sub.contourf(i, period, power, cmap = "Greys", zorder = 7)
cont  = sub.contour(i, period, sig95,\
                    [-99, 1], colors = "k", linewidths = 2, zorder = 9)
sub.set_yscale("log")
sub1.set_yscale("log")
sub.fill(np.concatenate([i, i[-1 :] + dt, i[-1 :] + dt,\
                         i[: 1] - dt, i[: 1] - dt]),\
         np.concatenate([coi, [1e-9], period[-1 :],
                         period[-1 :], [1e-9]]),
         "None", hatch = "x", edgecolor = "k", zorder = 11)
sub.set_xlim(min(i), max(i))
sub.set_ylim(min(period), max(period))
sub1.set_ylim(min(scales), max(scales))
sub.invert_yaxis()
sub1.invert_yaxis()
sub.set_xticks(np.arange(0, 200 + 20, 20))
sub.set_xlabel("i", fontdict = {"weight": "bold"})
sub.set_ylabel("Period", fontdict = {"weight": "bold"})
sub1.set_ylabel("Scale", fontdict = {"weight": "bold"})
sub.set_title("Wavelet Power Spectrum", fontdict = {"weight": "bold"})
fig.colorbar(contf, pad = 0.1, label = "Wavelet power")

fig, sub = plt.subplots(figsize = (10, 4))
sub1 = sub.twinx()
contf = sub.contourf(i, period, phase, cmap = "Greys", zorder = 7)
sub.set_yscale("log")
sub1.set_yscale("log")
sub.fill(np.concatenate([i, i[-1 :] + dt, i[-1 :] + dt,\
                         i[: 1] - dt, i[: 1] - dt]),\
         np.concatenate([coi, [1e-9], period[-1 :],
                         period[-1 :], [1e-9]]),
         "None", hatch = "x", edgecolor = "k", zorder = 11)
sub.set_xlim(min(i), max(i))
sub.set_ylim(min(period), max(period))
sub1.set_ylim(min(scales), max(scales))
sub.invert_yaxis()
sub1.invert_yaxis()
sub.set_xticks(np.arange(0, 200 + 20, 20))
sub.set_xlabel("i", fontdict = {"weight": "bold"})
sub.set_ylabel("Period", fontdict = {"weight": "bold"})
sub1.set_ylabel("Scale", fontdict = {"weight": "bold"})
sub.set_title("Phase", fontdict = {"weight": "bold"})
fig.colorbar(contf, pad = 0.1, label = "Phase")

fig, sub = plt.subplots(figsize = (10, 4))
sub.plot(fftfreqs, fft_power, c = "k", ls = "-")
sub.set_xlim(min(fftfreqs), max(fftfreqs))
sub.set_xlabel("Frequency", fontdict = {"weight": "bold"})
sub.set_ylabel(r"Phi []$^2$", fontdict = {"weight": "bold"})
sub.set_title("Power spectrum", fontdict = {"weight": "bold"})
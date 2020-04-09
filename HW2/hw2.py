# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 10:58:17 2019

@author: WHITE
"""
#%% IMPORT MODULES
import numpy as np
import matplotlib.pyplot as plt
#%% SET MATRIX
A = [[1, 2, 3], [3, 4, 5], [5, 6, 7]]
#%% FIND RANK OF MATRIX
rank = np.linalg.matrix_rank(A)
#%% FIND SUBSPACES
U, S, V = np.linalg.svd(A)
#%% HW
t  = np.arange(101)
xi = np.sin(2 * np.pi * t / 10)
yi = np.cos(2 * np.pi * t / 10)
a1, a0 = np.polyfit(xi, yi, deg = 1)

fig, sub = plt.subplots()
sub.plot(xi, yi, "ok", markersize = 5)
sub.plot(xi, a1 * xi + a0, "-r", lw = 1)
sub.set_xlim(min(xi), max(xi))
sub.set_xlabel("X")
sub.set_ylabel("Y")

r    = np.corrcoef(yi, a1 * xi + a0)
r    = r[0, 1]
n    = len(xi)
t    = (r * (n - 2) ** 0.5) / ((1 - r ** 2) ** 0.5)
z    = 0.5 * np.log((1 + r) / (1 - r))
s    = 1 / (n - 3) ** 0.5

mmin = z - 1.96 * s
mmax = z + 1.96 * s

rhomin = (np.exp(2 * mmin) - 1) / (np.exp(2 * mmin) + 1)
rhomax = (np.exp(2 * mmax) - 1) / (np.exp(2 * mmax) + 1)
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 09:09:39 2019

@author: WHITE
"""
#%% IMPORT MODULES
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import griddata
#%% MAKE FUNCTIONS
def FUNC1(x, t):
    y = np.zeros((len(x), len(t)), dtype = float)
    for xi, xx in enumerate(x):
        for ti, tt in enumerate(t):
            y[xi, ti] = (xx * np.cos(2 * np.pi * tt / 100) / 50)\
                        + (xx * np.cos(2 * np.pi * tt / 10) / 500)
    return y

def FUNC2(x, t):
    y = np.zeros((len(x), len(t)), dtype = float)
    for xi, xx in enumerate(x):
        for ti, tt in enumerate(t):
            y[xi, ti] = np.cos((2 * np.pi * xx / 10) - (2 * np.pi * tt /10))
    return y
#%% PLOT DATA FUNCTION
def PLOTDATA(x, Y, t, ts):
    fig, sub = plt.subplots()
    cs = np.array(np.linspace(0, 0.8, len(ts)), dtype = str)
    plots, labels = [], []
    for tt, cc in zip(ts, cs):
        ti = np.where(t == tt)[0][0]
        plot = sub.plot(x, Y[:, ti],\
                        c = cc, linestyle = "-", marker = "o", markersize = 3)
        plots.append(plot[0])
        labels.append("t = %i" % tt)
    sub.legend(plots, labels, loc = 1, prop = {"weight": "bold"})
    sub.set_xlabel("X", fontdict = {"weight": "bold"})
    sub.set_ylabel("Y", fontdict = {"weight": "bold"})
    sub.set_xlim(min(x), max(x))
    sub.set_ylim(np.min(Y) * 1.1, np.max(Y) * 1.1)
#%% EOF FUNCTION
def EOF(x, t, Y, mode_max, mode_num):
    n = np.shape(Y)[1]
    
    U, s, Vt = np.linalg.svd(Y)
    E        = U * -1
    S        = np.zeros(np.shape(Y))
    for si in range(len(s)):
        S[si, si] = s[si]
    Z   = np.matmul(S, Vt) * -1
    R   = np.matmul(Y, Y.T) / n
    
    L   = np.matmul(E.T, np.matmul(R, E))
    L2  = np.zeros(np.shape(L), dtype = float)
    Lm2 = np.zeros(np.shape(L), dtype = float)
    
    for ll in range(len(L)):
        L2[ll, ll]  = L[ll, ll] ** 0.5
        Lm2[ll, ll] = L[ll, ll] ** -0.5
    
    E = np.matmul(E, L2)
    Z = np.matmul(Lm2, Z)
    
    eigenvalues    = np.array([L[l, l] for l in range(len(L))]) 
    total_variance = np.sum(eigenvalues)
    percent        = 100 * eigenvalues / total_variance
    cumpercent     = [np.sum(percent[: i + 1]) for i in range(len(percent))]
    
    modes = np.arange(len(eigenvalues)) + 1
    
    fig, sub = plt.subplots()
    plot1 = sub.plot(modes[: mode_max], eigenvalues[: mode_max], "ko-")
    for i in range(mode_max):
        plt.text(modes[i] - 0.05, eigenvalues[i] + np.max(eigenvalues) * 0.02 ,\
                 "%.1f" % abs(eigenvalues[i]), ha = "right", va = "bottom")
    sub.set_xlim(0, mode_max + 1)
    sub.set_ylim(0, max(eigenvalues) * 1.2)
    sub.set_xticks(np.arange(mode_max) + 1)
    sub.set_xlabel("Mode", fontdict = {"weight": "bold"})
    sub.set_ylabel("Eigenvalue", fontdict = {"weight": "bold"})
    
    subt = plt.twinx(sub)
    plot2 = subt.plot(modes[: mode_max], percent[: mode_max],\
                      linestyle = "--", c = "grey", marker = "^")
    plot3 = subt.plot(modes[: mode_max], cumpercent[: mode_max],\
                      linestyle = "-", c = "grey", marker = "^")
    subt.set_ylim(0, 105)
    subt.set_ylabel("Percent (%)", fontdict = {"weight": "bold", "color": "grey"})
    
    sub.legend([plot1[0], plot2[0], plot3[0]],\
               ["Eigenvalue", "Percentage", "Cumulative"],\
               loc = "best", prop = {"weight": "bold"})
    
    fig, subs = plt.subplots(nrows = mode_num, ncols = 2)
    for mm in range(mode_num):
        sub1 = subs[mm, 0]
        sub1.plot(x, E[:, mm], "k-", lw = 1)
        sub1.set_ylim(np.min(E) * 1.1, np.max(E) * 1.1)
        sub1.set_title("EOF %i mode" % (mm + 1), fontdict = {"weight": "bold"})
        if x[0].dtype == float or x[0].dtype == int:
            sub1.set_xlim(min(x), max(x))
        else:
            sub1.set_xlim(x[0], x[-1])
        
        sub2 = subs[mm, 1]
        sub2.plot(t, Z[mm, :], "k-", lw = 1)
        sub2.set_ylim(np.min(Z) * 1.1, np.max(Z) * 1.1)
        sub2.set_title("PC %i mode" % (mm + 1), fontdict = {"weight": "bold"})
        if t[0].dtype == float or t[0].dtype == int:
            sub2.set_xlim(min(t), max(t))
        else:
            sub1.set_xlim(t[0], t[-1])
        if mm == mode_num - 1:
            sub1.set_xlabel("X", fontdict = {"weight": "bold"})
            sub2.set_xlabel("T", fontdict = {"weight": "bold"})
    fig.tight_layout()
    return E, Z
    
def RANDOM_EOF(x_num, n):
    x = np.arange(x_num) + 1
    t = np.arange(n) + 1
    
    Y = []
    for xn in range(x_num):
        np.random.seed((xn + 1) * n)
        Y = np.concatenate((Y, np.random.normal(0, 1, n)))
    Y = np.reshape(Y, (x_num, n))
    EOF(x = x, t = t, Y = Y, mode_max = x_num, mode_num = 3)
#%% HW3-1
x = np.arange(-50, 50 + 1, 1)
t = np.arange(1, 100 + 1, 1)
Y = FUNC1(x, t)
PLOTDATA(x = x, Y = Y, t = t, ts = [1, 10, 30, 50])
E, Z = EOF(x = x, t = t, Y = Y, mode_max = 5, mode_num = 3)
#%% HW3-2
x = np.arange(-50, 50 + 1, 1)
t = np.arange(1, 100 + 1, 1)
Y = FUNC2(x, t)
PLOTDATA(x = x, Y = Y, t = t, ts = [1, 5])
E, Z = EOF(x = x, t = t, Y = Y, mode_max = 5, mode_num = 3)
#%% HW3-3-1
RANDOM_EOF(x_num = 3, n = 10)
#%% HW3-3-2
RANDOM_EOF(x_num = 3, n = 100)
#%% HW3-3-3
RANDOM_EOF(x_num = 3, n = 1000)
#%% HW3-4
main_path = "C:/Users/WHITE/Google 드라이브/SNU/2019.1/객관적자료분석/HW3/DATA/"
station_info = np.loadtxt(main_path + "station.dat", dtype = str)
lons  = np.array(station_info[:, 0], dtype = float)
lats  = np.array(station_info[:, 1], dtype = float)
names = np.array(station_info[:, 2], dtype = str)

file = np.loadtxt(main_path + "Eta_bp.dat", dtype = float)
time = file[:, 0  ]
data = file[:, 1 :].T

E, Z = EOF(x = names, t = time, Y = data, mode_max = 10, mode_num = 3)

fig, subs = plt.subplots(nrows = 1, ncols = 3, figsize = (15, 5))
for mm, sub in enumerate(subs):
    m = Basemap(llcrnrlon = min(lons) - 3,\
                urcrnrlon = max(lons) + 3,\
                llcrnrlat = min(lats) - 3,\
                urcrnrlat = max(lats) + 3,\
                resolution = "l", ax = sub)
    m.drawcoastlines(linewidth = 0.5, color = "k")
    m.fillcontinents(color = "0.8", zorder = 5)
    m.drawmeridians(np.arange(0, 360, 5), labels = [0, 0, 0, 1],\
                    fontsize = 10, linewidth = 0, color = "None")
    m.drawparallels(np.arange(-90, 90, 5), labels = [1, 0, 0, 0],\
                    fontsize = 10, linewidth = 0, color = "None")
    
    xx, yy = m(lons, lats)
    m.scatter(xx, yy, color = "k", s = 5, zorder = 7)
    
    grid_lons = np.linspace(min(lons), max(lons), 100)
    grid_lats = np.linspace(min(lats), max(lats), 100)
    grid_lons, grid_lats = np.meshgrid(grid_lons, grid_lats)
    
    grid_data = griddata((lons, lats), E[:, mm], (grid_lons, grid_lats),\
                         method = "linear")
    xxx, yyy = m(grid_lons, grid_lats)
    conf = m.contourf(grid_lons, grid_lats, grid_data,\
                      cmap = "RdYlBu_r",\
                      vmin = -np.max(np.abs(E)), vmax = np.max(np.abs(E)))
    con  = m.contour(grid_lons, grid_lats, grid_data,\
                     colors = "k", linewidths = 0.5)
    plt.clabel(con, fontsize = 6, fmt = "%.2f")
    m.colorbar(conf)
    sub.set_title("EOF %i mode" % (mm + 1), fontdict = {"weight": "bold"})
fig.tight_layout()
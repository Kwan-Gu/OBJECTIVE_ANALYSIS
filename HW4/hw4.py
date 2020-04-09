# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 10:55:56 2019

@author: WHITE
"""
#%% IMPORT MODULES
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from mpl_toolkits.basemap import Basemap
#%% FUNCTIONS
def FUNC1(x, t):
    p = np.zeros((len(x), len(t)), dtype = float)
    q = np.zeros((len(x), len(t)), dtype = float)
    for i, xx in enumerate(x):
        for j, tt in enumerate(t):
            p[i, j] = (1 / 50 ) * xx * np.cos(2 * tt * np.pi / 100)
            q[i, j] = (1 / 500) * xx * np.cos(2 * tt * np.pi / 10 )
    return p, q # [x, t]

def FUNC2(x, t):
    p = np.zeros((len(x), len(t)), dtype = float)
    q = np.zeros((len(x), len(t)), dtype = float)
    for i, xx in enumerate(x):
        for j, tt in enumerate(t):
            p[i, j] = np.sin(2 * np.pi * xx / 10 - 2 * np.pi *tt / 10)
            q[i, j] = np.cos(2 * np.pi * xx / 10 - 2 * np.pi *tt / 10)
    return p, q # [x, t]
#%% PLOT DATA FUNCTION
def PLOTDATA(x, t, Y, ts):
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
#%% SVD FUNCTION
def SVD(x1, x2, Y1, Y2, t, mode_max, mode_num, svd_fig):
    if np.shape(Y1)[1] != np.shape(Y2)[1]:
        raise ValueError("N of Y1 and Y2 are different!")
    else:
        n = np.shape(Y1)[1]
    
    U, s, Vt = np.linalg.svd(np.matmul(Y1, Y2.T) / n)
    V = Vt.T
    S        = np.zeros((np.shape(Y1)[0], np.shape(Y2)[0]), dtype = float)
    for si in range(len(s)):
        S[si, si] = s[si]
    Z1   = np.matmul(U.T, Y1)
    Z2   = np.matmul(Vt, Y2)

    squared_s = np.power(s, 2)
    total_variance = np.sum(squared_s)
    percent        = 100 * squared_s / total_variance
    cumpercent     = [np.sum(percent[: i + 1]) for i in range(len(percent))]
    
    modes = np.arange(len(squared_s)) + 1
    
    fig, sub = plt.subplots()
    plot1 = sub.plot(modes[: mode_max], squared_s[: mode_max],\
                     "ko-", zorder = 3)
    for i in range(mode_max):
        plt.text(modes[i] - 0.05, squared_s[i] + np.max(squared_s) * 0.02 ,\
                 "%.1e" % abs(squared_s[i]),\
                 ha = "right", va = "bottom", zorder = 3)
    sub.set_xlim(0, mode_max + 1)
    sub.set_ylim(0, max(squared_s) * 1.2)
    sub.set_xticks(np.arange(mode_max) + 1)
    sub.set_xlabel("Mode", fontdict = {"weight": "bold"})
    sub.set_ylabel("Squared sigma", fontdict = {"weight": "bold"})
    
    subt = plt.twinx(sub)
    plot2 = subt.plot(modes[: mode_max], percent[: mode_max],\
                      linestyle = "--", c = "grey", marker = "^", zorder = 3)
    plot3 = subt.plot(modes[: mode_max], cumpercent[: mode_max],\
                      linestyle = "-", c = "grey", marker = "^", zorder = 3)
    subt.set_ylim(0, 105)
    subt.set_ylabel("Percent (%)",\
                    fontdict = {"weight": "bold", "color": "grey"})
    
    sub.legend([plot1[0], plot2[0], plot3[0]],\
               ["Squared sigma", "Percentage", "Cumulative"],\
               loc = "center right", prop = {"weight": "bold"})
    
    fig, subs = plt.subplots(figsize = (10, mode_num * 3),\
                             nrows = mode_num, ncols = 2)
    for mm in range(mode_num):
        sub1 = subs[mm, 0]
        if svd_fig == True:
            sub1.plot(x1, U[:, mm], "k-", lw = 1, label = "U")
            sub1.plot(x2, V[:, mm], "b-", lw = 0.5, label = "V")
            sub1.set_ylim(np.nanmin([U, V]) * 1.1, np.max([U, V]) * 1.1)
            sub1.set_title("SVD %i mode" % (mm + 1),\
                           fontdict = {"weight": "bold"})
            if x1[0].dtype == float or x1[0].dtype == int:
                sub1.set_xlim(np.nanmin([x1, x2]), np.nanmax([x1, x2]))
            else:
                sub1.set_xlim(x1[0], x1[-1])
            sub1.legend(prop = {"weight": "bold"})
        elif svd_fig == False:
            pass
        
        sub2 = subs[mm, 1]
        sub2.plot(t, Z1[mm, :], "k-", lw = 1, label = r"Eta")
        sub2.plot(t, Z2[mm, :], "b-", lw = 0.5, label = r"Pair")
        sub2.set_ylim(np.min(np.concatenate((Z1, Z2))) * 1.1,\
                      np.max(np.concatenate((Z1, Z2))) * 1.1)
        sub2.set_title("PC %i mode" % (mm + 1), fontdict = {"weight": "bold"})
        if t[0].dtype == float or t[0].dtype == int:
            sub2.set_xlim(min(t), max(t))
        else:
            sub2.set_xlim(t[0], t[-1])
        sub2.legend(prop = {"weight": "bold"})
        if mm == mode_num - 1:
            sub1.set_xlabel("X", fontdict = {"weight": "bold"})
            sub2.set_xlabel("T", fontdict = {"weight": "bold"})
    fig.tight_layout()
    return U, V, Z1, Z2
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
    subt.set_ylabel("Percent (%)",\
                    fontdict = {"weight": "bold", "color": "grey"})
    
    sub.legend([plot1[0], plot2[0], plot3[0]],\
               ["Eigenvalue", "Percentage", "Cumulative"],\
               loc = "best", prop = {"weight": "bold"})
    
    fig, subs = plt.subplots(nrows = mode_num, ncols = 2)
    for mm in range(mode_num):
        sub1 = subs[mm, 0]
        sub1.plot(x, E[:, mm], "k-", lw = 1)
        sub1.set_ylim(np.min(E) * 1.1, np.max(E) * 1.1)
        sub1.set_title("EOF %i mode" % (mm + 1),\
                       fontdict = {"weight": "bold"})
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
#%% HW4-1
x = np.arange(-50,  50 + 1, 1)
t = np.arange(  1, 100 + 1, 1)
p, q = FUNC1(x, t)

PLOTDATA(x = x, Y = p, t = t, ts = [1, 25, 50])
PLOTDATA(x = x, Y = q, t = t, ts = [1, 3, 5])
SVD(x1 = x, x2 = x, t = t, Y1 = p, Y2 = q,\
    mode_max = 5, mode_num = 3, svd_fig = True)
#%% HW4-2
x = np.arange(-50,  50 + 1, 1)
t = np.arange(  1, 100 + 1, 1)
p, q = FUNC2(x, t)

PLOTDATA(x = x, Y = p, t = t, ts = [1, 5])
PLOTDATA(x = x, Y = q, t = t, ts = [1, 5])
SVD(x1 = x, x2 = x, t = t, Y1 = p, Y2 = q,\
    mode_max = 5, mode_num = 3, svd_fig = True)
#%% HW4-3
main_path = "C:/Users/WHITE/Google 드라이브/SNU/2019.1/객관적자료분석/HW4/DATA/"

eta_loc   = np.loadtxt(main_path + "lonlat_Eta.dat", dtype = str)
eta_lons  = np.array(eta_loc[:, 0], dtype = float)
eta_lats  = np.array(eta_loc[:, 1], dtype = float)
eta_names = np.array(eta_loc[:, 2], dtype = str)

air_loc  = np.loadtxt(main_path + "lonlat_Pair.dat", dtype = float)
air_lons = air_loc[:, 0]
air_lats = air_loc[:, 1]

eta_file = np.loadtxt(main_path + "Eta.dat", dtype = float)
eta_time = eta_file[:, 0  ]
eta_data = eta_file[:, 1 :].T

air_file = np.loadtxt(main_path + "Pair.dat", dtype = float)
air_time = air_file[:, 0  ]
air_data = air_file[:, 1 :].T

E, P, EZ, PZ = SVD(x1 = eta_names, x2 = np.arange(len(air_lons)),\
                   Y1 = eta_data, Y2 = air_data, t = eta_time,\
                   mode_max = 5, mode_num = 3, svd_fig = False)

fig, subs = plt.subplots(nrows = 1, ncols = 3, figsize = (15, 5))
for mm, sub in enumerate(subs):
    m = Basemap(llcrnrlon = min(air_lons) - 3,\
                urcrnrlon = max(air_lons) + 3,\
                llcrnrlat = min(air_lats) - 3,\
                urcrnrlat = max(air_lats) + 3,\
                resolution = "l", ax = sub)
    m.drawcoastlines(linewidth = 0.5, color = "k")
    m.fillcontinents(color = "0.8", zorder = 5)
    m.drawmeridians(np.arange(0, 360, 5), labels = [0, 0, 0, 1],\
                    fontsize = 10, linewidth = 0, color = "None")
    m.drawparallels(np.arange(-90, 90, 5), labels = [1, 0, 0, 0],\
                    fontsize = 10, linewidth = 0, color = "None")
    
    eta_xx, eta_yy = m(eta_lons, eta_lats)
    m.scatter(eta_xx, eta_yy, color = "k", marker = "o", s = 5, zorder = 7)
    
    eta_grid_lons = np.linspace(min(eta_lons), max(eta_lons), 1000)
    eta_grid_lats = np.linspace(min(eta_lats), max(eta_lats), 1000)
    eta_grid_lons, eta_grid_lats = np.meshgrid(eta_grid_lons, eta_grid_lats)
    
    air_grid_lons = np.linspace(min(air_lons), max(air_lons), 100)
    air_grid_lats = np.linspace(min(air_lats), max(air_lats), 100)
    air_grid_lons, air_grid_lats = np.meshgrid(air_grid_lons, air_grid_lats)
        
    method = "linear"
    grid_eta = griddata((eta_lons, eta_lats), E[:, mm],\
                        (eta_grid_lons, eta_grid_lats),\
                        method = method)
    grid_air = griddata((air_lons, air_lats), P[:, mm],\
                        (air_grid_lons, air_grid_lats),\
                        method = method)
        
    eta_xxx, eta_yyy = m(eta_grid_lons, eta_grid_lats)
    air_xxx, air_yyy = m(air_grid_lons, air_grid_lats)
    eta_conf = m.contourf(eta_xxx, eta_yyy, grid_eta,\
                          cmap = "RdYlBu_r",\
                          vmin = -np.nanmax(np.abs(grid_eta)),\
                          vmax = np.nanmax(np.abs(grid_eta)))
    air_con  = m.contour(air_xxx, air_yyy, grid_air,\
                         colors = "k", linewidths = 0.5)
    plt.clabel(air_con, fontsize = 6, fmt = "%.2f")
    m.colorbar(eta_conf)
    sub.set_title("SVD %i mode" % (mm + 1), fontdict = {"weight": "bold"})
fig.tight_layout()
#%% HW4-3 EOF ETA
main_path = "C:/Users/WHITE/Google 드라이브/SNU/2019.1/객관적자료분석/HW4/DATA/"

eta_loc   = np.loadtxt(main_path + "lonlat_Eta.dat", dtype = str)
eta_lons  = np.array(eta_loc[:, 0], dtype = float)
eta_lats  = np.array(eta_loc[:, 1], dtype = float)
eta_names = np.array(eta_loc[:, 2], dtype = str)

eta_file = np.loadtxt(main_path + "Eta.dat", dtype = float)
eta_time = eta_file[:, 0  ]
eta_data = eta_file[:, 1 :].T

E, Z = EOF(x = eta_names, Y = eta_data, t = eta_time,\
           mode_max = 5, mode_num = 3)

fig, subs = plt.subplots(nrows = 1, ncols = 3, figsize = (15, 5))
for mm, sub in enumerate(subs):
    m = Basemap(llcrnrlon = min(eta_lons) - 3,\
                urcrnrlon = max(eta_lons) + 3,\
                llcrnrlat = min(eta_lats) - 3,\
                urcrnrlat = max(eta_lats) + 3,\
                resolution = "l", ax = sub)
    m.drawcoastlines(linewidth = 0.5, color = "k")
    m.fillcontinents(color = "0.8", zorder = 5)
    m.drawmeridians(np.arange(0, 360, 5), labels = [0, 0, 0, 1],\
                    fontsize = 10, linewidth = 0, color = "None")
    m.drawparallels(np.arange(-90, 90, 5), labels = [1, 0, 0, 0],\
                    fontsize = 10, linewidth = 0, color = "None")
    
    xx, yy = m(eta_lons, eta_lats)
    m.scatter(xx, yy, color = "k", s = 5, zorder = 7)
    
    grid_lons = np.linspace(min(eta_lons), max(eta_lons), 100)
    grid_lats = np.linspace(min(eta_lats), max(eta_lats), 100)
    grid_lons, grid_lats = np.meshgrid(grid_lons, grid_lats)
    
    grid_data = griddata((eta_lons, eta_lats), E[:, mm],\
                         (grid_lons, grid_lats),\
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
#%% HW4-3 EOF PAIR
main_path = "C:/Users/WHITE/Google 드라이브/SNU/2019.1/객관적자료분석/HW4/DATA/"

air_loc  = np.loadtxt(main_path + "lonlat_Pair.dat", dtype = float)
air_lons = air_loc[:, 0]
air_lats = air_loc[:, 1]

air_file = np.loadtxt(main_path + "Pair.dat", dtype = float)
air_time = air_file[:, 0  ]
air_data = air_file[:, 1 :].T

E, Z = EOF(x = np.arange(len(air_loc)), Y = air_data, t = air_time,\
           mode_max = 5, mode_num = 3)

fig, subs = plt.subplots(nrows = 1, ncols = 3, figsize = (15, 5))
for mm, sub in enumerate(subs):
    m = Basemap(llcrnrlon = min(air_lons) - 3,\
                urcrnrlon = max(air_lons) + 3,\
                llcrnrlat = min(air_lats) - 3,\
                urcrnrlat = max(air_lats) + 3,\
                resolution = "l", ax = sub)
    m.drawcoastlines(linewidth = 0.5, color = "k")
    m.fillcontinents(color = "0.8", zorder = 5)
    m.drawmeridians(np.arange(0, 360, 5), labels = [0, 0, 0, 1],\
                    fontsize = 10, linewidth = 0, color = "None")
    m.drawparallels(np.arange(-90, 90, 5), labels = [1, 0, 0, 0],\
                    fontsize = 10, linewidth = 0, color = "None")
    
    xx, yy = m(air_lons, air_lats)
    m.scatter(xx, yy, color = "k", s = 5, zorder = 7)
    
    grid_lons = np.linspace(min(air_lons), max(air_lons), 100)
    grid_lats = np.linspace(min(air_lats), max(air_lats), 100)
    grid_lons, grid_lats = np.meshgrid(grid_lons, grid_lats)
    
    grid_data = griddata((air_lons, air_lats), E[:, mm],\
                         (grid_lons, grid_lats),\
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
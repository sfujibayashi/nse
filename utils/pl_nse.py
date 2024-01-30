import numpy as np
#import h5py
import sys
#import subprocess
import os
from scipy.interpolate import interp2d

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.ticker as ticker
#plt.switch_backend('agg')
#os.environ["PATH"] += os.pathsep + '/data/home/sfujibayashi/libs/texlive/2020/bin/x86_64-linux'

params={
    'font.size'           : 24.0     ,
    'font.family'         : 'DeJaVu Sans'  ,
#    'font.family'         : 'Times New Roman'  ,
    'xtick.major.size'    : 2        ,
    'xtick.major.width'   : 1.5      ,
    'xtick.labelsize'     : 24.0     ,
    'xtick.direction'     : 'out'     ,
    'ytick.major.size'    : 2        ,
    'ytick.major.width'   : 1.5      ,
    'ytick.labelsize'     : 24.0     ,
    'ytick.direction'     : 'out'     ,
    'xtick.major.pad'     : 2        ,
    'xtick.minor.pad'     : 2        ,
    'ytick.major.pad'     : 2        ,
    'ytick.minor.pad'     : 2        ,
    'axes.linewidth'      : 1.5      ,
    'text.usetex'         : True    }

plt.rcParams.update(params)

clight=2.99792458e10

fn = "nse"

rho =np.unique(np.loadtxt(fn,comments="#",usecols=0))
temp=np.unique(np.loadtxt(fn,comments="#",usecols=1))
ye  =np.unique(np.loadtxt(fn,comments="#",usecols=2))
nrho = len(rho)
ntemp= len(temp)
nye   = len(ye)
print(rho)
print(temp)
print(ye)
print(nrho,ntemp,nye)

mexc_av = np.loadtxt(fn,comments="#",usecols=3).reshape(nye,ntemp,nrho)
abar    = np.loadtxt(fn,comments="#",usecols=5).reshape(nye,ntemp,nrho)

size_fig=10.0
fig = plt.figure(figsize=(size_fig, size_fig*0.625))
ax = fig.add_subplot(111)

#fig2 = plt.figure(figsize=(size_fig, size_fig*0.625))
#ax2 = fig2.add_subplot(111)

cmap = cm.get_cmap("jet",256)
norm = colors.Normalize(vmin=0.0,vmax=150.0)

iye = 0
print(ye[iye])
abar_slice = abar[iye,:,:]
ax.text()
ax.pcolormesh(rho,temp,abar_slice,cmap=cmap  ,norm=norm)

# ax2.pcolormesh(rho,temp,abar_slice,cmap=cmap  ,norm=norm)
# for itemp in range(ntemp):
#     for irho in range(nrho):
#         col = cmap(norm(abar[iye,itemp,irho]))
#         ax2.scatter(rho[irho],temp[itemp],color=col)

#ax.set_xlim()
#ax.set_ylim(-r_max/runi_plot,r_max/runi_plot)
s_map = cm.ScalarMappable(norm=norm, cmap=cmap)
s_map.set_array([0.0])

cbar=fig.colorbar(s_map,pad=0.01)
cbar.set_label("$\\langle A\\rangle_\\mathrm{h}$",rotation = -90, labelpad= 30)


ax.set_xlim(1e3,1e12)
ax.set_ylim(1e8,1e11)
ax.set_xscale("log")
ax.set_yscale("log")

#ax2.set_xscale("log")
#ax2.set_yscale("log")

#ax.set_xlabel(xlabel)
#ax.xaxis.set_major_locator(ticker.MultipleLocator(xmajorl))
#ax.xaxis.set_minor_locator(ticker.MultipleLocator(xminorl))

plt.show()

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
    'xtick.direction'     : 'in'     ,
    'ytick.major.size'    : 2        ,
    'ytick.major.width'   : 1.5      ,
    'ytick.labelsize'     : 24.0     ,
    'ytick.direction'     : 'in'     ,
    'xtick.major.pad'     : 2        ,
    'xtick.minor.pad'     : 2        ,
    'ytick.major.pad'     : 2        ,
    'ytick.minor.pad'     : 2        ,
    'axes.linewidth'      : 1.5      ,
    'text.usetex'         : True    }

plt.rcParams.update(params)

clight=2.99792458e10

with open("fort.99",mode="r") as f:
    array = f.readline().split()
    print(array)
    rho = float(array[1])
    temp= float(array[2])
    ye  = float(array[3])
    array = f.readline().split()
    x_n = float(array[1])
    x_p = float(array[2])
    x_y = float(array[3])
    

data=np.loadtxt("fort.98",comments="#")
xn_history=data[:,1]
xp_history=data[:,2]

data=np.loadtxt("fort.99",comments="#")

xn=np.unique(data[:,0])
xp=np.unique(data[:,1])

nxn = len(xn)
nxp = len(xp)

logxsum = data[:,2].reshape(nxp,nxn)
logysum = data[:,3].reshape(nxp,nxn)
dxn = data[:,4].reshape(nxp,nxn)
dxp = data[:,5].reshape(nxp,nxn)

# print(nxn,nxp)
# #sys.exit()
# for ixn in range(nxn):
#     for ixp in range(nxp):
#         print(("%12.4e"*4) % (xn[ixn],xp[ixp],logxsum[ixp,ixn],logysum[ixp,ixn]) )
# sys.exit()

dxn_unit= dxn/np.sqrt(dxn**2+dxp**2+1e-20)
dxp_unit= dxp/np.sqrt(dxn**2+dxp**2+1e-20)

fig = plt.figure(figsize=(13.0, 13.0*0.5))
##ax1 = fig.add_subplot(121)
##ax2 = fig.add_subplot(122)
buf=0.05
bar=0.03
ax1 = fig.add_axes([buf ,buf ,0.4,0.9])
ax2 = fig.add_axes([0.50,buf ,0.4,0.9])
axc = fig.add_axes([0.9+0.01 , buf ,bar ,0.9])

pcm = ax1.pcolormesh(xn,xp,logxsum,cmap="seismic",shading='auto',norm=colors.Normalize(vmin=-3.0, vmax=3.0) )
pcm = ax2.pcolormesh(xn,xp,logysum,cmap="seismic",shading='auto',norm=colors.Normalize(vmin=-3.0, vmax=3.0) )
ax1.contour(xn,xp,logxsum,levels=[0.0],colors="g")
ax1.contour(xn,xp,logysum,levels=[0.0],colors="orange")
ax2.contour(xn,xp,logxsum,levels=[0.0],colors="g")
ax2.contour(xn,xp,logysum,levels=[0.0],colors="orange")
ax1.plot(xn_history,xp_history,lw=1.0,marker="o",markersize=2,color="y")
ax2.plot(xn_history,xp_history,lw=1.0,marker="o",markersize=2,color="y")
ax1.plot(xn_history[-1],xp_history[-1],lw=2.0,marker="*",color="y")
ax2.plot(xn_history[-1],xp_history[-1],lw=2.0,marker="*",color="y")

skip=10
st=int(skip/2)
quiv1= ax1.quiver(xn[st::skip], xp[st::skip], dxn_unit[st::skip,st::skip], dxp_unit[st::skip,st::skip],scale=20,scale_units='width',minlength=1,width=0.002,headwidth=4,headlength=6)
quiv2= ax2.quiver(xn[st::skip], xp[st::skip], dxn_unit[st::skip,st::skip], dxp_unit[st::skip,st::skip],scale=20,scale_units='width',minlength=1,width=0.002,headwidth=4,headlength=6)

props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax1.text(0.99,0.99,"$\\rho=$%9.1e, $T=$%9.1e, $Y_e=$%5.3f" % (rho,temp,ye), ha="right", va="top", transform=ax1.transAxes, fontsize=14, bbox=props)
ax2.text(0.99,0.99,"$X_n=$%9.2e, $X_p=$%9.2e, $X_\\alpha=$%9.2e" % (x_n,x_p,x_y), ha="right", va="top", transform=ax2.transAxes, fontsize=14, bbox=props)

fig.subplots_adjust(left=0.20,bottom=0.20,right=0.95, top=0.95)
#ax1.set_aspect('equal')
#ax2.set_aspect('equal')
#ax.tick_params(axis='both',width=1.5,length=6,which='major',top=True,right=True,pad=10.0)
#ax.tick_params(axis='both',width=1.5,length=3,which='minor',top=True,right=True)
ax1.set_xlabel("$x_n = (\\mu_n-m_\\mathrm{u}c^2)/kT/\\ln 10$")
ax1.set_ylabel("$x_p = (\\mu_p-m_\\mathrm{u}c^2)/kT/\\ln 10$")
ax2.set_xlabel("$x_n = (\\mu_n-m_\\mathrm{u}c^2)/kT/\\ln 10$")
#ax2.set_ylabel("$x_p = (\\mu_p-m_\\mathrm{u}c^2)/kT/\\ln 10$")

cbar=fig.colorbar(pcm,cax=axc)

#ax.set_xlim(-r_max/runi_plot,r_max/runi_plot)
#ax.set_ylim(-r_max/runi_plot,r_max/runi_plot)

#ax.xaxis.set_major_locator(ticker.MultipleLocator(xmajorl))
#ax.xaxis.set_minor_locator(ticker.MultipleLocator(xminorl))

fig.savefig("conv_rho%7.1e_temp%7.1e_ye%5.3f.png" % (rho,temp,ye))
plt.show()

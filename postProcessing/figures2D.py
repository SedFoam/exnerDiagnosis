#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Filename: figures2D.py
Date: 2025-03-27
Version: 1.0
License: GPL-3.0-or-later
Copyright 2025 Alban Gilletta, Cyrille Bonamy, Julien Chauchat
This file is part of ExnerAnalysis.

ExnerAnalysis is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

ExnerAnalysis is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with ExnerAnalysis. If not, see <https://www.gnu.org/licenses/>.

Contact: alban.gilletta-de-saint-joseph@univ-grenoble-alpes.fr
"""

#
# Import section
#

import subprocess
import sys
import numpy as np
import os
import fluidfoam
from pylab import figure, subplot, axis, xlabel, ylabel, show, savefig, plot
from pylab import title, matplotlib
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, LogNorm
import matplotlib.cm
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from massBalance2D import Twodimsimu

case = "2D"
basepath = "../Configuration/"

meshCreated = True          # True if the file read_mesh.h5 has been created
resultCreated = True        # True if the file read_massBalance'+str(self.asint).replace('.','p')+'.h5 has been created

nya = 64
nyb = 200 
read_point_nc = True        # True if the file pointerpostproc.nc has been created
save_point_nc = False       # False if the file pointerpostproc.nc has been saved

simu1 = Twodimsimu(basepath, case, nya, nyb, read_point_nc, save_point_nc, meshCreated, resultCreated, 0.57)
simu2 = Twodimsimu(basepath, case, nya, nyb, True         , False        , True       , resultCreated, 0.58)
simu3 = Twodimsimu(basepath, case, nya, nyb, True         , False        , True       , resultCreated, 0.59)
simu4 = Twodimsimu(basepath, case, nya, nyb, True         , False        , True       , resultCreated, 0.60)
simu5 = Twodimsimu(basepath, case, nya, nyb, True         , False        , True       , resultCreated, 0.61)
simu6 = Twodimsimu(basepath, case, nya, nyb, True         , False        , True       , resultCreated, 0.62)

# Output result at time 50 is required for figure 5

#tprofile = '50'
#alpha = fluidfoam.readscalar(basepath+case, tprofile, "alpha.a", structured=True, verbose=False, precision=13)
#aUa = fluidfoam.readvector(basepath+case, tprofile, "alphaUa", structured=True, verbose=False, precision=13)

matplotlib.rc('text', usetex = True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath} \usepackage{amssymb} \usepackage{physics} \usepackage[dvipsnames]{xcolor}')

# plot creation and parameters
matplotlib.rcParams.update({'font.size': 12})
mpl.rcParams["lines.linewidth"] = 1.5
mpl.rcParams["lines.markersize"] = 10

plt.figure(figsize=[12, 15])
gs = plt.GridSpec(1, 1)
gs.update(left=0.08, right=0.975, top=0.95, bottom=0.15, wspace=0.1,
          hspace=0.1)
font_legend = matplotlib.font_manager.FontProperties(family='Comic Sans MS',
                                                     weight='bold',
                                                     style='normal', size=14)
fignum1 = plt.figure(num=1,figsize=(12,15),dpi=300)
ax1 = fignum1.add_subplot(111)
plt.xlim(-4e-8,4e-8) 
plt.ylim(-4e-8,4e-8)
plt.xlabel(r'\textcircled{1} Bed evolution (m s$^{-3}$)',fontsize = '12') 
plt.ylabel(r'-(\textcircled{2} Storage evolution+\textcircled{3} Lateral flux+\textcircled{4} Bed flux+\textcircled{5} Top flux) (m$^{3}$ s$^{-1}$)', fontsize = '12') 
plt.grid()

fignum2 = plt.figure(num=2,figsize=(12,15),dpi=300)
ax2 = fignum2.add_subplot(111)
plt.xlim(0,150)
plt.ylim(-3.5e-8,3.5e-8)
plt.xlabel(r'Time (s)',fontsize = '12')
plt.ylabel(r'Volume evolution (m$^{3}$ s$^{-1}$)', fontsize = '12')
plt.grid()

fignum3 = plt.figure(num=3,figsize=(12,15),dpi=300)
ax3 = fignum3.add_subplot(111)
plt.xlim(0,150)
plt.ylim(-2e-8,0.5e-8)
plt.xlabel(r'Time (s)',fontsize = '12')
plt.ylabel(r'Residual (m s$^{-3}$)', fontsize = '12')
plt.grid()
ax3_0 = ax3.inset_axes([0.25,0.05,0.4,0.25])
ax3_0.axis([0, 150, -125, 125])
ax3_0.grid(linestyle=':')
ax3_0.set_ylabel(r'Residual ($\%$)', fontsize = '12')

fignum4 = plt.figure(num=4,figsize=(12,15),dpi=300)
ax4 = fignum4.add_subplot(111)
plt.xlim(0,150)
plt.ylim(-2e-8,0.5e-8)
plt.xlabel(r'Time (s)',fontsize = '12')
plt.ylabel(r'Volume evolution (m$^{3}$ s$^{-1}$)', fontsize = '12')
plt.grid()

#fignum5 = plt.figure(num=5,figsize=(12,15),dpi=300)
#ax5 = fignum5.add_subplot(111)
#plt.xlim(0.05,0.15)
#plt.ylim(-1.5e-10,0.75e-10)
#plt.xlabel(r'x (m)',fontsize = '12')
#plt.ylabel(r'$-\phi u_z^s \Delta x \Delta y$ (m$^{3}$ s$^{-1}$)', fontsize = '12')
#plt.grid()

#
# FIGURE 1
#

color = simu4.time[0:-1]
sc1 = ax1.scatter(simu4.phi_dzbed_sur_dt, -(simu4.intdVolPhi_sur_dt + simu4.intSlat_PindS + simu4.intS_PindS + simu4.intS2_PindS + simu4.intStop_PindS) ,marker ='o', s=50, c=color, cmap='seismic', edgecolors='none', vmin= simu4.time[0], vmax = simu4.time[-2])
ax1.plot(np.linspace(-10,10,100),np.linspace(-10,10,100),'--k',label ='y=x')
ax1.plot(np.linspace(-20,20,10),np.linspace(0,0,10),'-k', lw = 1)
ax1.plot(np.linspace(0,0,10),np.linspace(-20,20,10),'-k', lw = 1)
cbar1 = plt.colorbar(sc1,ax=ax1)
cbar1.set_label(r'Time (s)',fontsize= '12')
ax1.legend(prop={"size": 12})

#
# FIGURE 2
#
# Forward [:-1] Backward [1:]

ax2.plot(simu4.time[:-1], simu4.phi_dzbed_sur_dt , color = 'orange', label = r'\textcircled{1} Bed evolution')
ax2.plot(simu4.time[:-1], simu4.intdVolPhi_sur_dt , color = 'b', label = r'\textcircled{2} Storage evolution')
ax2.plot(simu4.time[:-1], simu4.intSlat_PindS, color = 'g', label = r' \textcircled{3} Lateral flux')
ax2.plot(simu4.time[:-1], simu4.intS_PindS + simu4.intS2_PindS, color = 'm', label = r' \textcircled{4} Bed flux')
ax2.plot(simu4.time[:-1], simu4.intStop_PindS, color = 'c', label = r' \textcircled{5} Top flux')
ax2.plot(simu4.time[:-1], simu4.phi_dzbed_sur_dt + simu4.intdVolPhi_sur_dt + simu4.intSlat_PindS + simu4.intS_PindS + simu4.intS2_PindS + simu4.intStop_PindS, color = 'r', label = r'\textcircled{1}+\textcircled{2}+\textcircled{3}+\textcircled{4}+\textcircled{4}')
ax2.legend(prop={"size": 12})

#
# FIGURE 3
#

ax3.plot(simu1.time[:-1], simu1.phi_dzbed_sur_dt + simu1.intdVolPhi_sur_dt + simu1.intSlat_PindS + simu1.intS_PindS + simu1.intS2_PindS + simu1.intStop_PindS, color = 'r', label = r'$\phi_b=0.57$' )
ax3.plot(simu2.time[:-1], simu2.phi_dzbed_sur_dt + simu2.intdVolPhi_sur_dt + simu2.intSlat_PindS + simu2.intS_PindS + simu2.intS2_PindS + simu2.intStop_PindS, color = 'b', label = r'$\phi_b=0.58$' )
ax3.plot(simu3.time[:-1], simu3.phi_dzbed_sur_dt + simu3.intdVolPhi_sur_dt + simu3.intSlat_PindS + simu3.intS_PindS + simu3.intS2_PindS + simu3.intStop_PindS, color = 'g', label = r'$\phi_b=0.59$' )
ax3.plot(simu4.time[:-1], simu4.phi_dzbed_sur_dt + simu4.intdVolPhi_sur_dt + simu4.intSlat_PindS + simu4.intS_PindS + simu4.intS2_PindS + simu4.intStop_PindS, color = 'm', label = r'$\phi_b=0.60$' )
ax3.plot(simu5.time[:-1], simu5.phi_dzbed_sur_dt + simu5.intdVolPhi_sur_dt + simu5.intSlat_PindS + simu5.intS_PindS + simu5.intS2_PindS + simu5.intStop_PindS, color = 'c', label = r'$\phi_b=0.61$' )
ax3.plot(simu6.time[:-1], simu6.phi_dzbed_sur_dt + simu6.intdVolPhi_sur_dt + simu6.intSlat_PindS + simu6.intS_PindS + simu6.intS2_PindS + simu6.intStop_PindS, color = 'y', label = r'$\phi_b=0.62$' )
ax3.legend(prop={"size": 12})

ax3_0.plot(simu1.time[:-1], (simu1.phi_dzbed_sur_dt + simu1.intdVolPhi_sur_dt + simu1.intSlat_PindS + simu1.intS_PindS + simu1.intS2_PindS + simu1.intStop_PindS)/simu1.phi_dzbed_sur_dt*100, color = 'r' )
ax3_0.plot(simu2.time[:-1], (simu2.phi_dzbed_sur_dt + simu2.intdVolPhi_sur_dt + simu2.intSlat_PindS + simu2.intS_PindS + simu2.intS2_PindS + simu2.intStop_PindS)/simu2.phi_dzbed_sur_dt*100, color = 'b' )
ax3_0.plot(simu3.time[:-1], (simu3.phi_dzbed_sur_dt + simu3.intdVolPhi_sur_dt + simu3.intSlat_PindS + simu3.intS_PindS + simu3.intS2_PindS + simu3.intStop_PindS)/simu3.phi_dzbed_sur_dt*100, color = 'g' )
ax3_0.plot(simu4.time[:-1], (simu4.phi_dzbed_sur_dt + simu4.intdVolPhi_sur_dt + simu4.intSlat_PindS + simu4.intS_PindS + simu4.intS2_PindS + simu4.intStop_PindS)/simu4.phi_dzbed_sur_dt*100, color = 'm' )
ax3_0.plot(simu5.time[:-1], (simu5.phi_dzbed_sur_dt + simu5.intdVolPhi_sur_dt + simu5.intSlat_PindS + simu5.intS_PindS + simu5.intS2_PindS + simu5.intStop_PindS)/simu5.phi_dzbed_sur_dt*100, color = 'c' )
ax3_0.plot(simu6.time[:-1], (simu6.phi_dzbed_sur_dt + simu6.intdVolPhi_sur_dt + simu6.intSlat_PindS + simu6.intS_PindS + simu6.intS2_PindS + simu6.intStop_PindS)/simu6.phi_dzbed_sur_dt*100, color = 'y' )

#
# FIGURE 4
#

ax4.plot(simu1.time[:-1], simu1.intS_PindS, ls = '-', color = 'r', label = r'$\phi_b=0.57$')
ax4.plot(simu2.time[:-1], simu2.intS_PindS, ls = '-', color = 'b', label = r'$\phi_b=0.58$')
ax4.plot(simu3.time[:-1], simu3.intS_PindS, ls = '-', color = 'g', label = r'$\phi_b=0.59$')
ax4.plot(simu4.time[:-1], simu4.intS_PindS, ls = '-', color = 'm', label = r'$\phi_b=0.60$')
ax4.plot(simu5.time[:-1], simu5.intS_PindS, ls = '-', color = 'c', label = r'$\phi_b=0.61$')
ax4.plot(simu6.time[:-1], simu6.intS_PindS, ls = '-', color = 'y', label = r'$\phi_b=0.62$')

ax4.plot(simu1.time[:-1], simu1.intS2_PindS, ls = '--', color = 'r')
ax4.plot(simu2.time[:-1], simu2.intS2_PindS, ls = '--', color = 'b')
ax4.plot(simu3.time[:-1], simu3.intS2_PindS, ls = '--', color = 'g')
ax4.plot(simu4.time[:-1], simu4.intS2_PindS, ls = '--', color = 'm')
ax4.plot(simu5.time[:-1], simu5.intS2_PindS, ls = '--', color = 'c')
ax4.plot(simu6.time[:-1], simu6.intS2_PindS, ls = '--', color = 'y')
ax4.legend(prop={"size": 12})

#
# FIGURE 5
#

#postProcess
#aWa57, aWa58, aWa59, aWa60, aWa61, aWa62 = np.zeros((6, simu1.Nx, simu1.Nt))
#for i in range(simu1.imin,simu1.imax+1):
#    for t in range(simu1.Nt-1):
#        bedcondi57 = np.where(alpha[i,:,0] <= simu1.asint)
#        bedcondi58 = np.where(alpha[i,:,0] <= simu2.asint)
#        bedcondi59 = np.where(alpha[i,:,0] <= simu3.asint)
#        bedcondi60 = np.where(alpha[i,:,0] <= simu4.asint)
#        bedcondi61 = np.where(alpha[i,:,0] <= simu5.asint)
#        bedcondi62 = np.where(alpha[i,:,0] <= simu6.asint)
#        aWa57[i,t] = - (aUa[1, i, bedcondi57[0][0], 0]*simu1.deltaY[bedcondi57[0][0]-1]+aUa[1, i, bedcondi57[0][0]-1, 0]*simu1.deltaY[bedcondi57[0][0]])/(simu1.deltaY[bedcondi57[0][0]]+simu1.deltaY[bedcondi57[0][0]-1])*simu1.deltaX*simu1.deltaZ
#        aWa58[i,t] = - (aUa[1, i, bedcondi58[0][0], 0]*simu1.deltaY[bedcondi58[0][0]-1]+aUa[1, i, bedcondi58[0][0]-1, 0]*simu1.deltaY[bedcondi58[0][0]])/(simu1.deltaY[bedcondi58[0][0]]+simu1.deltaY[bedcondi58[0][0]-1])*simu1.deltaX*simu1.deltaZ
#        aWa59[i,t] = - (aUa[1, i, bedcondi59[0][0], 0]*simu1.deltaY[bedcondi59[0][0]-1]+aUa[1, i, bedcondi59[0][0]-1, 0]*simu1.deltaY[bedcondi59[0][0]])/(simu1.deltaY[bedcondi59[0][0]]+simu1.deltaY[bedcondi59[0][0]-1])*simu1.deltaX*simu1.deltaZ
#        aWa60[i,t] = - (aUa[1, i, bedcondi60[0][0], 0]*simu1.deltaY[bedcondi60[0][0]-1]+aUa[1, i, bedcondi60[0][0]-1, 0]*simu1.deltaY[bedcondi60[0][0]])/(simu1.deltaY[bedcondi60[0][0]]+simu1.deltaY[bedcondi60[0][0]-1])*simu1.deltaX*simu1.deltaZ
#        aWa61[i,t] = - (aUa[1, i, bedcondi61[0][0], 0]*simu1.deltaY[bedcondi61[0][0]-1]+aUa[1, i, bedcondi61[0][0]-1, 0]*simu1.deltaY[bedcondi61[0][0]])/(simu1.deltaY[bedcondi61[0][0]]+simu1.deltaY[bedcondi61[0][0]-1])*simu1.deltaX*simu1.deltaZ

#ax5.plot(np.arange(simu1.imin,simu1.imax+1)*simu1.deltaX+simu1.deltaX/2, aWa57[simu1.imin:simu1.imax+1,50], color = 'r', label = r'$\phi_b=0.57$')
#ax5.plot(np.arange(simu1.imin,simu1.imax+1)*simu1.deltaX+simu1.deltaX/2, aWa58[simu1.imin:simu1.imax+1,50], color = 'b', label = r'$\phi_b=0.58$')
#ax5.plot(np.arange(simu1.imin,simu1.imax+1)*simu1.deltaX+simu1.deltaX/2, aWa59[simu1.imin:simu1.imax+1,50], color = 'g', label = r'$\phi_b=0.59$')
#ax5.plot(np.arange(simu1.imin,simu1.imax+1)*simu1.deltaX+simu1.deltaX/2, aWa60[simu1.imin:simu1.imax+1,50], color = 'm', label = r'$\phi_b=0.60$')
#ax5.plot(np.arange(simu1.imin,simu1.imax+1)*simu1.deltaX+simu1.deltaX/2, aWa61[simu1.imin:simu1.imax+1,50], color = 'c', label = r'$\phi_b=0.61$')
#ax5.plot(np.arange(simu1.imin,simu1.imax+1)*simu1.deltaX+simu1.deltaX/2, aWa62[simu1.imin:simu1.imax+1,50], color = 'y', label = r'$\phi_b=0.62$')

#ax5.legend(prop={"size": 12})

for i in plt.get_fignums():
    if i == 1:
        plt.figure(i).set_size_inches(10,8)
        plt.figure(i).savefig('Figures/2D/'+'massBalance_2D60'+'.eps', 
                facecolor='w', edgecolor='w', dpi=300, format='eps')
    if i == 2:
        plt.figure(i).set_size_inches(10,8)
        plt.figure(i).savefig('Figures/2D/'+'contrib_2D60'+'.eps', 
                facecolor='w', edgecolor='w', dpi=300, format='eps')
    if i == 3:
        plt.figure(i).set_size_inches(10,8)
        plt.figure(i).savefig('Figures/2D/'+'residual_2D'+'.eps', 
                facecolor='w', edgecolor='w', dpi=300, format='eps')
    if i == 4:
        plt.figure(i).set_size_inches(10,8)
        plt.figure(i).savefig('Figures/2D/'+'contribBedFluxComponent_2D'+'.eps', 
                facecolor='w', edgecolor='w', dpi=300, format='eps')
#    if i == 5:
#        plt.figure(i).set_size_inches(10,8)
#        plt.figure(i).savefig('Figure/2D/'+'profilebedFluxVert_2D'+'.eps', 
#                facecolor='w', edgecolor='w', dpi=300, format='eps')

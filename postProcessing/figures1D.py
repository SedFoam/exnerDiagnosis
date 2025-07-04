#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Filename: figures1D.py
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
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from massBalance1D import Onedimsimu

case = "1D"
basepath = "../Configuration/"

meshCreated = True           # True if the file read_mesh.h5 has been created
resultCreated = True        # True if the file read_massBalance'+str(self.asint).replace('.','p')+'.h5 has been created

simu57 = Onedimsimu(basepath, case, meshCreated, resultCreated, 0.57)
simu58 = Onedimsimu(basepath, case, meshCreated, resultCreated, 0.58)
simu59 = Onedimsimu(basepath, case, meshCreated, resultCreated, 0.59)
simu60 = Onedimsimu(basepath, case, meshCreated, resultCreated, 0.60)

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
plt.xlim(-1e-11,1e-11) 
plt.ylim(-1e-11,1e-11)
plt.xlabel(r'\textcircled{1} Bed evolution (m$^{3}$ s$^{-1}$)',fontsize = '12') 
plt.ylabel(r'-(\textcircled{2} Storage evolution+\textcircled{4} Bed flux+\textcircled{5} Top flux) (m$^{3}$ s$^{-1}$)', fontsize = '12') 
plt.grid()

fignum2 = plt.figure(num=2,figsize=(12,15),dpi=300)
ax2 = fignum2.add_subplot(111)
plt.xlim(0,1800)
plt.ylim(-1.1e-11,1.1e-11)
plt.xlabel(r'Time (s)',fontsize = '12')
plt.ylabel(r'Volume evolution (m$^{3}$ s$^{-1}$)', fontsize = '12')
plt.grid()

fignum3 = plt.figure(num=3,figsize=(12,15),dpi=300)
ax3 = fignum3.add_subplot(111)
plt.xlim(0,1800)
plt.ylim(-7.5e-13,1e-12)
plt.xlabel(r'Time (s)',fontsize = '12')
plt.ylabel(r'Residual (m$^{3}$ s$^{-1}$)', fontsize = '12')
plt.grid()
ax3_0 = ax3.inset_axes([0.1,0.72,0.4,0.25])
ax3_0.axis([0, 1800, -1.5e1, 1.5e1])
ax3_0.grid(linestyle=':')
ax3_0.set_ylabel(r'Residual ($\%$)', fontsize = '12')
        
#
# FIGURE 1
#

color = simu59.time[0:-1]
sc1 = ax1.scatter(simu59.phi_dzbed_sur_dt, -(simu59.intdVolPhi_sur_dt + simu59.intS_PindS + simu59.intStop_PindS), marker ='o', s=50, c=color, cmap='seismic', edgecolors='none', vmin= simu59.time[0], vmax = simu59.time[-2])
ax1.plot(np.linspace(-10,10,100),np.linspace(-10,10,100),'--k',label ='y=x')
ax1.plot(np.linspace(-20,20,10),np.linspace(0,0,10),'-k', lw = 1)
ax1.plot(np.linspace(0,0,10),np.linspace(-20,20,10),'-k', lw = 1)
cbar1 = plt.colorbar(sc1,ax=ax1)
cbar1.set_label(r'Time (s)',fontsize= '12')
ax1.legend(prop={"size": 12})

#
# FIGURE 2
#
# If you use Euler forward temporal scheme -> Forward [:-1]; Euler backward temporal scheme -> Backward [1:]

ax2.plot(simu59.time[:-1], simu59.phi_dzbed_sur_dt, color = 'orange', label = r'\textcircled{1} Bed evolution')
ax2.plot(simu59.time[:-1], simu59.intdVolPhi_sur_dt, color = 'b', label = r'\textcircled{2} Storage evolution')
ax2.plot(simu59.time[:-1], simu59.intS_PindS, color = 'm', label = r' \textcircled{4} Bed flux')
ax2.plot(simu59.time[:-1], simu59.intStop_PindS, color = 'c', label = r' \textcircled{5} Top flux')
ax2.plot(simu59.time[:-1],simu59.phi_dzbed_sur_dt + simu59.intdVolPhi_sur_dt + simu59.intS_PindS + simu59.intStop_PindS, color = 'r', label = r'\textcircled{1}+\textcircled{2}+\textcircled{4}+\textcircled{5}')
ax2.legend(prop={"size": 12})

#
# FIGURE 3
#

ax3.plot(simu57.time[:-1], simu57.phi_dzbed_sur_dt + simu57.intdVolPhi_sur_dt + simu57.intS_PindS + simu57.intStop_PindS, color = 'r', label = r'$\phi_b=0.57$' )
ax3.plot(simu58.time[:-1], simu58.phi_dzbed_sur_dt + simu58.intdVolPhi_sur_dt + simu58.intS_PindS + simu58.intStop_PindS, color = 'b', label = r'$\phi_b=0.58$' )
ax3.plot(simu59.time[:-1], simu59.phi_dzbed_sur_dt + simu59.intdVolPhi_sur_dt + simu59.intS_PindS + simu59.intStop_PindS, color = 'g', label = r'$\phi_b=0.59$' )
ax3.plot(simu60.time[:-1], simu60.phi_dzbed_sur_dt + simu60.intdVolPhi_sur_dt + simu60.intS_PindS + simu60.intStop_PindS, color = 'm', label = r'$\phi_b=0.60$' )
ax3.legend(prop={"size": 12})

ax3_0.plot(simu57.time[:-1], (simu57.phi_dzbed_sur_dt + simu57.intdVolPhi_sur_dt + simu57.intS_PindS + simu57.intStop_PindS)/np.mean(simu57.phi_dzbed_sur_dt)*100, color = 'r' )
ax3_0.plot(simu58.time[:-1], (simu58.phi_dzbed_sur_dt + simu58.intdVolPhi_sur_dt + simu58.intS_PindS + simu58.intStop_PindS)/np.mean(simu58.phi_dzbed_sur_dt)*100, color = 'b' )
ax3_0.plot(simu59.time[:-1], (simu59.phi_dzbed_sur_dt + simu59.intdVolPhi_sur_dt + simu59.intS_PindS + simu59.intStop_PindS)/np.mean(simu59.phi_dzbed_sur_dt)*100, color = 'g' )
ax3_0.plot(simu60.time[:-1], (simu60.phi_dzbed_sur_dt + simu60.intdVolPhi_sur_dt + simu60.intS_PindS + simu60.intStop_PindS)/np.mean(simu60.phi_dzbed_sur_dt)*100, color = 'm' )

for i in plt.get_fignums():
    if i == 1:
        plt.figure(i).set_size_inches(10,8)
        plt.figure(i).savefig('Figures/1D/'+'massBalance_1D59'+'.eps', 
                facecolor='w', edgecolor='w', dpi=300, format='eps')
    if i == 2:
        plt.figure(i).set_size_inches(10,8)
        plt.figure(i).savefig('Figures/1D/'+'contrib_1D59'+'.eps', 
                facecolor='w', edgecolor='w', dpi=300, format='eps')
    if i == 3:
        plt.figure(i).set_size_inches(10,8)
        plt.figure(i).savefig('Figures/1D/'+'residual_1D'+'.eps', 
                facecolor='w', edgecolor='w', dpi=300, format='eps')

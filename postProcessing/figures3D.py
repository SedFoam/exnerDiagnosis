#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Filename: figures3D.py
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
from netCDF4 import Dataset
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from massBalance3D import ThreedimSqrsimu
from func1 import create_point_2Dcyl, save_point, read_point
from func1 import uns2cartvec, uns2cart

caseM1 = "M1"
caseM2 = "M2"
caseM3 = "M3"
caseM4 = "M4"
caseM5 = "M5"

caseList = [caseM1, caseM2, caseM3, caseM4, caseM5]

gZ = 41.432951063305104
nzb = 52
gZ1 = 0.5830419509076277 
nz1 = 176 

read_point_nc = True        # True if the file pointerpostproc.nc has been created
save_point_nc = False       # False if the file pointerpostproc.nc has been saved
meshCreated = True          # True if the file read_mesh.h5 has been created
resultCreated = True        # True if the file read_massBalance'+str(self.asint).replace('.','p')+'.h5 has been created

basepath = "../Configuration/3D/"

simuList = []
k=-1
for case in caseList:
    k = k + 1
    if case == caseM1 :
        g0X = 1/2.75889735683314
        N0X = 12
        N1X = 6
        g2X = 4.936196796988305
        N2X = 18
        g0Y = 1/3.7513348259705044
        N0Y = 15
        N1Y = 9
        g2Y = 3.7513348259705044
        N2Y = 15
    elif case == caseM2 :
        g0X = 1/2.7869985823022327
        N0X = 24
        N1X = 12
        g2X = 4.99980614393634
        N2X = 36
        g0Y = 1/2.7869985823022327
        N0Y = 24
        N1Y = 12
        g2Y = 2.7869985823022327
        N2Y = 24
    elif case == caseM3 :
        g0X = 1/2.7960591876441585
        N0X = 48
        N1X = 24
        g2X = 5.02351891521504
        N2X = 72
        g0Y = 1/2.7960591876441585
        N0Y = 48
        N1Y = 24
        g2Y = 2.7960591876441585
        N2Y = 48
    elif case == caseM4 :
        g0X = 1/2.8005003010678777
        N0X = 96
        N1X = 48
        g2X = 14.151419136811377
        N2X = 72
        g0Y = 1/8.388849212076593
        N0Y = 48
        N1Y = 48
        g2Y = 8.388849212076593
        N2Y = 48
    elif case == caseM5 :
        g0X = 1/2.8019677508239584
        N0X = 144
        N1X = 72
        g2X = 24.845683120716316
        N2X = 72
        g0Y = 1/15.080304387910225
        N0Y = 48
        N1Y = 72
        g2Y = 15.080304387910225
        N2Y = 48
    simuList.append(ThreedimSqrsimu(basepath, case, g0X, g2X, N0X, N1X, N2X, g0Y, g2Y, N0Y, N1Y, N2Y, gZ, gZ1, nzb, nz1, read_point_nc, save_point_nc, meshCreated, resultCreated, 0.60))

simu1 = simuList[0]
simu2 = simuList[1]
simu3 = simuList[2]
simu4 = simuList[3]
simu5 = simuList[4]

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
plt.xlim(-3e-5,0e-5) 
plt.ylim(-3e-5,0e-5)
plt.xlabel(r'\textcircled{1} (m$^{3}$ s$^{-1}$)',fontsize = '20')
plt.ylabel(r'-(\textcircled{2} + \textcircled{3} + \textcircled{4} + \textcircled{5}) (m$^{3}$ s$^{-1}$)', fontsize = '20') 
plt.tick_params(axis='x', labelsize=20)
plt.tick_params(axis='y', labelsize=20)
plt.grid()

fignum2 = plt.figure(num=2,figsize=(12,15),dpi=300)
ax2 = fignum2.add_subplot(111)
plt.xlim(0,20)
plt.ylim(-2.5e-5,2.5e-5)
plt.xlabel(r'Time (s)',fontsize = '20')
plt.ylabel(r'Volume evolution (m$^{3}$ s$^{-1}$)', fontsize = '20')
plt.tick_params(axis='x', labelsize=20)
plt.tick_params(axis='y', labelsize=20)
plt.grid()

fignum3 = plt.figure(num=3,figsize=(12,15),dpi=300)
ax3 = fignum3.add_subplot(111)
plt.xlim(0,20)
plt.ylim(-1e-5,1e-5)
plt.xlabel(r'Time (s)',fontsize = '20')
plt.ylabel(r'Residual (m$^{3}$ s$^{-1}$)', fontsize = '20')
plt.tick_params(axis='x', labelsize=20)
plt.tick_params(axis='y', labelsize=20)
plt.grid()
ax3_0 = ax3.inset_axes([0.55,0.05,0.4,0.35])
ax3_0.axis([0, 20, -3.5e2, 6e1])
ax3_0.grid(linestyle=':')
ax3_0.set_ylabel(r'Residual ($\%$)', fontsize = '12')

fignum4 = plt.figure(num=4,figsize=(12,15),dpi=300)
ax4 = fignum4.add_subplot(111)
plt.xlim(0,20)
plt.ylim(-2.5e-5,2.5e-5)
plt.xlabel(r'Time (s)',fontsize = '20')
plt.ylabel(r'Volume evolution (m$^{3}$ s$^{-1}$)', fontsize = '20')
plt.tick_params(axis='x', labelsize=20)
plt.tick_params(axis='y', labelsize=20)
plt.grid()

#
# FIGURE 1
#

color = simu1.time[:-1]
sc1 = ax1.scatter(simu1.phi_dzbed_sur_dt, -(simu1.intdVolPhi_sur_dt + simu1.intSlat_PindS + simu1.intS_PindS + simu1.intS2_PindS + simu1.intS3_PindS + simu1.intStop_PindS), marker ='x', s=30, c=color, cmap='seismic', edgecolors='none', vmin= simu3.time[0], vmax = simu3.time[-1],label = 'M1')
ax1.scatter(simu2.phi_dzbed_sur_dt, -(simu2.intdVolPhi_sur_dt + simu2.intSlat_PindS + simu2.intS_PindS + simu2.intS2_PindS + simu2.intS3_PindS + simu2.intStop_PindS), marker ='+', s=30, c=color, cmap='seismic', edgecolors='none', vmin= simu3.time[0], vmax = simu3.time[-1],label = 'M2')
ax1.scatter(simu3.phi_dzbed_sur_dt, -(simu3.intdVolPhi_sur_dt + simu3.intSlat_PindS + simu3.intS_PindS + simu3.intS2_PindS + simu3.intS3_PindS + simu3.intStop_PindS), marker ='^', s=30, c=color, cmap='seismic', edgecolors='none', vmin= simu3.time[0], vmax = simu3.time[-1],label = 'M3')
ax1.scatter(simu4.phi_dzbed_sur_dt, -(simu4.intdVolPhi_sur_dt + simu4.intSlat_PindS + simu4.intS_PindS + simu4.intS2_PindS + simu4.intS3_PindS + simu4.intStop_PindS), marker ='*', s=30, c=color, cmap='seismic', edgecolors='none', vmin= simu3.time[0], vmax = simu3.time[-1],label = 'M4')
ax1.scatter(simu5.phi_dzbed_sur_dt, -(simu5.intdVolPhi_sur_dt + simu5.intSlat_PindS + simu5.intS_PindS + simu5.intS2_PindS + simu5.intS3_PindS + simu5.intStop_PindS), marker ='o', s=30, c=color, cmap='seismic', edgecolors='none', vmin= simu3.time[0], vmax = simu3.time[-1],label = 'M5')
ax1.plot(np.linspace(-10,10,100),np.linspace(-10,10,100),'--k',label ='y=x')
ax1.plot(np.linspace(-20,20,10),np.linspace(0,0,10),'-k', lw = 1)
ax1.plot(np.linspace(0,0,10),np.linspace(-20,20,10),'-k', lw = 1)
cbar1 = plt.colorbar(sc1,ax=ax1)
cbar1.ax.tick_params(labelsize=20)
cbar1.set_label(r'Time (s)',fontsize= '20')
ax1.legend(prop={"size": 17})

#
# FIGURE 2
#

ax2.plot(np.linspace(0,20,10),np.linspace(0,0,10),'-k', lw = 1)
ax2.plot(simu4.time[:-1], simu4.phi_dzbed_sur_dt , color = 'orange', label = r'\textcircled{1} Bed evolution')
ax2.plot(simu4.time[:-1], simu4.intdVolPhi_sur_dt , color = 'b', label = r'\textcircled{2} Storage evolution')
ax2.plot(simu4.time[:-1], simu4.intSlat_PindS, color = 'g', label = r' \textcircled{3} Lateral flux')
ax2.plot(simu4.time[:-1], simu4.intS_PindS + simu4.intS2_PindS + simu4.intS3_PindS, color = 'm', label = r' \textcircled{4} Bed flux')
ax2.plot(simu4.time[:-1], simu4.intStop_PindS, color = 'c', label = r' \textcircled{5} Top flux')
ax2.plot(simu4.time[:-1], simu4.phi_dzbed_sur_dt + simu4.intdVolPhi_sur_dt + simu4.intSlat_PindS + simu4.intS_PindS + simu4.intS2_PindS + simu4.intS3_PindS + simu4.intStop_PindS, color = 'r', label = r'\textcircled{1}+\textcircled{2}+\textcircled{3}+\textcircled{4}+\textcircled{5}')
ax2.legend(prop={"size": 17}, ncol = 2)

#
# FIGURE 3
#

ax3.plot(np.linspace(0,20,10),np.linspace(0,0,10),'-k', lw = 1)
ax3.plot(simu1.time[:-1], simu1.phi_dzbed_sur_dt + simu1.intdVolPhi_sur_dt + simu1.intSlat_PindS + simu1.intS_PindS + simu1.intS2_PindS + simu1.intS3_PindS + simu1.intStop_PindS, color = 'r', label = 'M1' )
ax3.plot(simu2.time[:-1], simu2.phi_dzbed_sur_dt + simu2.intdVolPhi_sur_dt + simu2.intSlat_PindS + simu2.intS_PindS + simu2.intS2_PindS + simu2.intS3_PindS + simu2.intStop_PindS, color = 'b', label = 'M2' )
ax3.plot(simu3.time[:-1], simu3.phi_dzbed_sur_dt + simu3.intdVolPhi_sur_dt + simu3.intSlat_PindS + simu3.intS_PindS + simu3.intS2_PindS + simu3.intS3_PindS + simu3.intStop_PindS, color = 'g', label = 'M3' )
ax3.plot(simu4.time[:-1], simu4.phi_dzbed_sur_dt + simu4.intdVolPhi_sur_dt + simu4.intSlat_PindS + simu4.intS_PindS + simu4.intS2_PindS + simu4.intS3_PindS + simu4.intStop_PindS, color = 'm', label = 'M4' )
ax3.plot(simu5.time[:-1], simu5.phi_dzbed_sur_dt + simu5.intdVolPhi_sur_dt + simu5.intSlat_PindS + simu5.intS_PindS + simu5.intS2_PindS + simu5.intS3_PindS + simu5.intStop_PindS, color = 'c', label = 'M5' )
ax3.legend(prop={"size": 17})

ax3_0.plot(simu1.time[:-1], (simu1.phi_dzbed_sur_dt + simu1.intdVolPhi_sur_dt + simu1.intSlat_PindS + simu1.intS_PindS + simu1.intS2_PindS + simu1.intS3_PindS + simu1.intStop_PindS)/simu1.phi_dzbed_sur_dt*100, color = 'r' )
ax3_0.plot(simu2.time[:-1], (simu2.phi_dzbed_sur_dt + simu2.intdVolPhi_sur_dt + simu2.intSlat_PindS + simu2.intS_PindS + simu2.intS2_PindS + simu2.intS3_PindS + simu2.intStop_PindS)/simu2.phi_dzbed_sur_dt*100, color = 'b' )
ax3_0.plot(simu3.time[:-1], (simu3.phi_dzbed_sur_dt + simu3.intdVolPhi_sur_dt + simu3.intSlat_PindS + simu3.intS_PindS + simu3.intS2_PindS + simu3.intS3_PindS + simu3.intStop_PindS)/simu3.phi_dzbed_sur_dt*100, color = 'g' )
ax3_0.plot(simu4.time[:-1], (simu4.phi_dzbed_sur_dt + simu4.intdVolPhi_sur_dt + simu4.intSlat_PindS + simu4.intS_PindS + simu4.intS2_PindS + simu4.intS3_PindS + simu4.intStop_PindS)/simu4.phi_dzbed_sur_dt*100, color = 'm' )
ax3_0.plot(simu5.time[:-1], (simu5.phi_dzbed_sur_dt + simu5.intdVolPhi_sur_dt + simu5.intSlat_PindS + simu5.intS_PindS + simu5.intS2_PindS + simu5.intS3_PindS + simu5.intStop_PindS)/simu5.phi_dzbed_sur_dt*100, color = 'c' )

#
# FIGURE 4
#

ax4.plot(np.linspace(0,20,10),np.linspace(0,0,10),'-k', lw = 1)
ax4.plot(simu4.time[:-1], simu4.intSlatin_PindS, color = 'g', ls = '--', label = r' \textcircled{3a} Inflow')
ax4.plot(simu4.time[:-1], simu4.intSlatleft_PindS, color = 'g', ls = '-.', label = r' \textcircled{3b} Left outflow')
ax4.plot(simu4.time[:-1], simu4.intSlatright_PindS, color = 'g', ls = ':', label = r' \textcircled{3c} Right outflow')
ax4.plot(simu4.time[:-1], simu4.intS_PindS, color = 'm', ls = '--', label = r' \textcircled{4a} Vertical bed flux')
ax4.plot(simu4.time[:-1], simu4.intS2_PindS, color = 'm', ls = '-.',label = r' \textcircled{4b} Longitudinal bed flux')
ax4.plot(simu4.time[:-1], simu4.intS3_PindS, color = 'm', ls = ':',label = r' \textcircled{4c} Spanwise bed flux')
ax4.legend(prop={"size": 17})

for i in plt.get_fignums():
    if i == 1:
        plt.figure(i).set_size_inches(10,8)
        plt.figure(i).savefig('Figures/3D/'+'massBalance3D'+'.eps', 
                facecolor='w', edgecolor='w', dpi=300, format='eps')
    if i == 2:
        plt.figure(i).set_size_inches(10,8)
        plt.figure(i).savefig('Figures/3D/'+'contrib3D_60_'+caseM4+'.eps', 
                facecolor='w', edgecolor='w', dpi=300, format='eps')
    if i == 3:
        plt.figure(i).set_size_inches(10,8)
        plt.figure(i).savefig('Figures/3D/'+'residual3D_60'+'.eps', 
                facecolor='w', edgecolor='w', dpi=300, format='eps')
    if i == 4:
        plt.figure(i).set_size_inches(10,8)
        plt.figure(i).savefig('Figures/3D/'+'bedLatFlux_60_'+caseM4+'.eps', 
                facecolor='w', edgecolor='w', dpi=300, format='eps') 

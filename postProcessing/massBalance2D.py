#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Filename: massBalance2D.py
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
import fluidfoam
from tqdm import tqdm
import h5py

from func1 import divergence, ddt, create_point_1Dcart, save_1Dpoint, read_1Dpoint

class Twodimsimu(object):
    def __init__(self, basepath, case, nya, nyb, read_point_nc, save_point_nc, read_mesh, read_result, asint):
        self.basepath = basepath
        self.case = case
        self.sol = self.basepath + self.case + "/"
        self.solsav = self.sol + 'constant'
        self.asint = asint      # Solid volume fraction criteria
        self.xmin = 0.05        # left limit of the control volume
        self.xmax = 0.15        # right limit of the control volume
        self.deltaZ = 0.002
        self.read_point_nc = read_point_nc
        self.save_point_nc = save_point_nc
        
        self.createmesh(nya, nyb, read_mesh, read_result)
        
    def createmesh(self, nya, nyb, read_mesh, read_result):
        if not read_mesh:
            X, Y, Z = fluidfoam.readmesh(self.sol, structured=True, precision=13)
            Xb, Yb, Zb = fluidfoam.readmesh(self.sol, precision=13)
            hf = h5py.File(self.solsav+'/read_mesh.h5', 'w')
            hf.create_dataset('X',data=X)
            hf.create_dataset('Y',data=Y)
            hf.create_dataset('Z',data=Z)
            hf.create_dataset('Xb',data=Xb)
            hf.create_dataset('Yb',data=Yb)
            hf.create_dataset('Zb',data=Zb)
            hf.close()
        else :
            hf = h5py.File(self.solsav+'/read_mesh.h5', 'r')
            X = np.array(hf.get('X'))
            Y = np.array(hf.get('Y'))
            Z = np.array(hf.get('Z'))
            Xb = np.array(hf.get('Xb'))
            Yb = np.array(hf.get('Yb'))
            Zb = np.array(hf.get('Zb'))
        
        Nx, Ny, Nz = np.shape(X)
        self.deltaX=1.0/Nx
        ncell = np.size(X)
        dVs = np.zeros(ncell)

        if not self.read_point_nc:
            n1d, ny, pbed, pointer = create_point_1Dcart(Xb, Yb)
            if self.save_point_nc:
                save_1Dpoint(self.solsav, n1d, ny, pbed, Xb, Yb, pointer)
        else:
            n1d, ny, pbed, pointer = read_1Dpoint(self.solsav)

        imin = np.where(X[:, 0, 0] >= self.xmin)[0][0]      # left limit of the control volume
        imax = np.where(X[:, 0, 0] <= self.xmax)[-1][-1]    # right limit of the control volume
        
        dVs[:] = fluidfoam.readscalar(self.sol, "constant", "V", verbose=False, precision=13)
        varorig = dVs[:]
        dV = varorig[pointer].reshape(Nx, Ny, Nz)
        print(np.sum(varorig), np.sum(dV))
  
        #
        # Construct vertical mesh
        #

        # Upper part (water)
        gya = 100 # grading factor
        Lya = 0.15 # Length

        # Lower part (sediment)
        gyb = 0.2461
        Lyb = 0.05

        dya, dyaN = fluidfoam.meshdesign.getdzs(Lya, gya, nya)
        yai, dya, gYa = fluidfoam.meshdesign.getgz(Lya, dya, nya)
        dyb, dybN = fluidfoam.meshdesign.getdzs(Lyb, gyb, nyb)
        ybi, dyb, gYb = fluidfoam.meshdesign.getgz(Lyb, dyb, nyb)
        gYb=1/gYb

        deltaY = np.concatenate((dyb, dya))     # Vector length vertical cells from bottom to top
        yp = np.concatenate((ybi-Lyb, yai[1:])) # Vector elevation each edge cells
        
        self.X = X
        self.Y = Y
        self.Z = Z
        self.imin = imin
        self.imax = imax
        self.ncell = ncell
        self.n1d = n1d
        self.Nx = Nx
        self.Ny = Ny
        self.Nz = Nz
        self.pbed = pbed
        self.pointer = pointer
        self.deltaY = deltaY
        self.yp = yp
        self.dV = dV
    
        self.createtime(read_result)
        
    def createtime(self, read_result):
        try:
            proc = subprocess.Popen(["foamListTimes", "-case", self.sol], stdout=subprocess.PIPE)
        except FileNotFoundError:
            print("foamListTimes : command not found")
            print("Do you have load OpenFoam environement?")
            sys.exit(0)
        output = proc.stdout.read()
        tread = output.decode().rstrip().split("\n")
        if len(tread) != 0 :
            dt = float(tread[1])-float(tread[0])
            Nt = len(tread)
        else :
            tread = range(1,151,1)
            tread = np.array(tread).astype(str)
            dt = 1
            Nt = 150
        time = np.zeros(Nt)
        
        self.time = time
        self.dt = dt
        self.Nt = Nt
        self.tread = tread
                
        self.postProcess(read_result)
            
    def postProcess(self, read_result):
        k = -1
        if not read_result:
            alpha = np.zeros((self.Nx, self.Ny, self.Nz, self.Nt))
            aUa = np.zeros((3, self.Nx, self.Ny, self.Nz, self.Nt))
            zbed, dzbed_sur_dt, dVolPhi_sur_dt, VoldPhi_sur_dt = np.zeros((4, self.Nx, self.Nz, self.Nt))
            phi_dzbed_sur_dt, intdVolPhi_sur_dt, intVolPhi, intVoldPhi_sur_dt = np.zeros((4, self.Nt))
            intS_PindS, intS2_PindS, intSlat_PindS, intStop_PindS = np.zeros((4, self.Nt))
                    
            for t in self.tread:
                print("Reading time: %s s" % t)
                k = k + 1
                alpha[:, :, :, k] = fluidfoam.readscalar(self.sol, t, "alpha.a", structured=True, verbose=False, precision=13)
                aUa[:, :, :, :, k] = fluidfoam.readvector(self.sol, t, "alphaUa", structured=True, verbose=False, precision=13)
                self.time[k] = float(t)
            
            for i in range(self.imin,self.imax+1):
                for k in range(self.Nz):
                    for t in range(self.Nt):
                        bedcondi = np.where(alpha[i,:,k,t] <= self.asint)
                        zbed[i, k, t] = self.yp[bedcondi[0][0]]
                        intVolPhi[t] = intVolPhi[t] + np.sum(alpha[i,bedcondi,k,t]*self.dV[i,bedcondi,k])
                    dzbed_sur_dt[i, k, :] = ddt('eulerFor', zbed[i, k, :], self.dt)
            
            intdVolPhi_sur_dt = ddt('eulerFor', intVolPhi, self.dt)
            
            for i in range(self.imin,self.imax+1):
                for k in range(self.Nz):
                    for t in range(self.Nt):
                        bedcondi = np.where(alpha[i,:,k,t] <= self.asint)
                        bedcondi1 = np.where(alpha[i+1,:,k,t] <= self.asint) 
                        if bedcondi[0][0] == 0 :
                            phi_dzbed_sur_dt[t] = phi_dzbed_sur_dt[t] + alpha[i, bedcondi[0][0], k, t]*dzbed_sur_dt[i, k, t]*self.deltaZ*self.deltaX
                        else:
                            phi_dzbed_sur_dt[t] = phi_dzbed_sur_dt[t] + (alpha[i, bedcondi[0][0], k, t]*self.deltaY[bedcondi[0][0]-1]+alpha[i, bedcondi[0][0]-1, k, t]*self.deltaY[bedcondi[0][0]])/(self.deltaY[bedcondi[0][0]-1]+self.deltaY[bedcondi[0][0]])*dzbed_sur_dt[i, k, t]*self.deltaZ*self.deltaX
                            intS_PindS[t] = intS_PindS[t] - (aUa[1, i, bedcondi[0][0], k, t]*self.deltaY[bedcondi[0][0]-1]+aUa[1, i, bedcondi[0][0]-1, k, t]*self.deltaY[bedcondi[0][0]])/(self.deltaY[bedcondi[0][0]]+self.deltaY[bedcondi[0][0]-1])*self.deltaX*self.deltaZ
                        intStop_PindS[t] = intStop_PindS[t] + aUa[1, i, -1, k, t]*self.deltaX*self.deltaZ
                        if i!=self.imax:
                            dz = (zbed[i+1, k, t] - zbed[i, k, t])
                            match dz:
                                case dz if dz > 0:
                                    S2 = np.where(np.logical_and(alpha[i,:,k,t] <= self.asint, alpha[i+1,:,k,t] > self.asint))[0]
                                    intS2_PindS[t] = intS2_PindS[t] + np.sum((aUa[0, i, S2, k, t]+aUa[0, i+1, S2, k, t])*0.5*self.deltaZ*self.deltaY[S2])
                                case dz if dz < 0:
                                    S2 = np.where(np.logical_and(alpha[i,:,k,t] > self.asint, alpha[i+1,:,k,t] <= self.asint))[0]
                                    intS2_PindS[t] = intS2_PindS[t] - np.sum((aUa[0, i, S2, k, t]+aUa[0, i+1, S2, k, t])*0.5*self.deltaZ*self.deltaY[S2])

            for k in range(self.Nz):
                for t in range(self.Nt):
                    bedcondi = np.where(alpha[self.imin,:,k,t] <= self.asint)
                    for j in range(bedcondi[0][0],self.Ny):
                        intSlat_PindS[t] = intSlat_PindS[t] - (aUa[0, self.imin, j, k, t]+aUa[0, self.imin-1, j, k, t])*0.5*self.deltaZ*self.deltaY[j]

            for k in range(self.Nz):
                for t in range(self.Nt):
                    bedcondi = np.where(alpha[self.imax,:,k,t] <= self.asint)
                    for j in range(bedcondi[0][0],self.Ny):
                        intSlat_PindS[t] = intSlat_PindS[t] + (aUa[0, self.imax, j, k, t] + aUa[0, self.imax +1, j, k, t])*0.5*self.deltaZ*self.deltaY[j]         
            
            hf = h5py.File(self.solsav+'/read_massBalance'+str(self.asint).replace('.','p')+'.h5', 'w')
            hf.create_dataset('phi_dzbed_sur_dt',data=phi_dzbed_sur_dt)
            hf.create_dataset('intdVolPhi_sur_dt',data=intdVolPhi_sur_dt)
            hf.create_dataset('intS_PindS',data=intS_PindS)
            hf.create_dataset('intS2_PindS',data=intS2_PindS)
            hf.create_dataset('intStop_PindS',data=intStop_PindS)
            hf.create_dataset('intSlat_PindS',data=intSlat_PindS)
            hf.close()
        else :
            for t in tqdm(self.tread):
                k = k + 1
                self.time[k] = float(t)
            hf = h5py.File(self.solsav+'/read_massBalance'+str(self.asint).replace('.','p')+'.h5', 'r')
            intVoldPhi_sur_dt = np.array(hf.get('intVoldPhi_sur_dt'))
            phi_dzbed_sur_dt = np.array(hf.get('phi_dzbed_sur_dt'))
            intdVolPhi_sur_dt = np.array(hf.get('intdVolPhi_sur_dt'))
            intS_PindS = np.array(hf.get('intS_PindS'))
            intS2_PindS = np.array(hf.get('intS2_PindS'))
            intStop_PindS = np.array(hf.get('intStop_PindS'))
            intSlat_PindS = np.array(hf.get('intSlat_PindS'))
        
        # Forward [:-1] Backward [1:] for temporal derivative
        self.phi_dzbed_sur_dt = phi_dzbed_sur_dt[:-1]           # term 1 : bed evolution
        self.intdVolPhi_sur_dt = intdVolPhi_sur_dt[:-1]         # term 2 : storage evolution
        self.intSlat_PindS = intSlat_PindS[:-1]                 # term 3 : lateral flux
        self.intS_PindS = intS_PindS[:-1]                       # term 4a : bed flux -> vertical component
        self.intS2_PindS = intS2_PindS[:-1]                     # term 4b : bed flux -> horizontal component
        self.intStop_PindS = intStop_PindS[:-1]                 # term 5 : top flux

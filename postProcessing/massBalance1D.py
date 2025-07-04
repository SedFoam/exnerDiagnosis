#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Filename: massBalance1D.py
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

from func1 import divergence, ddt

class Onedimsimu(object):
    def __init__(self, basepath, case, read_mesh, read_result, asint):
        self.basepath = basepath
        self.case = case
        self.sol = self.basepath + self.case + "/"
        self.solsav = self.sol + 'constant'
        self.asint = asint # Solid volume fraction criteria
        self.Sbh = 0.0006**2 # Bed surface
        
        self.createmesh(read_mesh, read_result)

    def createmesh(self, read_mesh, read_result):
        if not read_mesh:
            X, Y, Z = fluidfoam.readmesh(self.sol)
            Ny = np.size(X)
            hf = h5py.File(self.solsav+'/read_mesh.h5', 'w')
            hf.create_dataset('X',data=X)
            hf.create_dataset('Y',data=Y)
            hf.create_dataset('Z',data=Z)
            hf.close()
        else :
            hf = h5py.File(self.solsav+'/read_mesh.h5', 'r')
            X = np.array(hf.get('X'))
            Y = np.array(hf.get('Y'))
            Z = np.array(hf.get('Z'))

        Ny = np.size(X)
        dV = np.zeros(Ny)
        dV[:] = fluidfoam.readscalar(self.sol, "constant", "V", verbose=False, precision=13)
        
        self.X = X
        self.Y = Y
        self.Z = Z
        self.Ny = Ny
        self.dV = dV
        
        self.createtime(read_mesh, read_result)
        
    def createtime(self, read_mesh, read_result):
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
            tread = range(20,1820,20)
            tread = np.array(tread).astype(str)
            dt = 20
            Nt = 90
        time = np.zeros(Nt)
        
        self.time = time
        self.dt = dt
        self.Nt = Nt
        self.tread = tread
        
        self.postProcess(read_result)
        
    def postProcess(self, read_result):
        k = -1
        if not read_result:
            alpha = np.zeros((self.Ny, self.Nt))
            aUa = np.zeros((3, self.Ny,  self.Nt))
            zbed, intVolPhi, phi_dzbed_sur_dt, intdVolPhi_sur_dt, intS_PindS, intStop_PindS = np.zeros((6, self.Nt))
            
            k = -1
            for t in self.tread:
                print("Reading time: %s s" % t)
                k = k + 1
                alpha[:, k] = fluidfoam.readscalar(self.sol, t, "alpha.a")
                aUa[:, :, k] = fluidfoam.readvector(self.sol, t, "alphaUa")
                self.time[k] = float(t)
            
            for i in range(self.Nt):
                bedcondi = np.where(alpha[:, i] <= self.asint)
                if len(bedcondi[0]) == self.Ny :
                    zbed[i] = 0.0
                else :
                    zbed[i] = 0.5*(self.Y[bedcondi[0][0]]+self.Y[bedcondi[0][0]-1])
                intVolPhi[i] = np.sum(alpha[bedcondi,i]*self.dV[bedcondi])

            intdVolPhi_sur_dt = ddt('eulerFor', intVolPhi, self.dt)
            dzbed_sur_dt = ddt('eulerFor', zbed, self.dt)

            for i in range(0,self.Nt):
                bedcondi = np.where(alpha[:, i] <= self.asint)
                if len(bedcondi[0]) == self.Ny :
                    phi_dzbed_sur_dt[i] = alpha[bedcondi[0][0],i]*dzbed_sur_dt[i]*self.Sbh
                    intS_PindS[i] = 0
                else :
                    phi_dzbed_sur_dt[i] = 0.5*(alpha[bedcondi[0][0],i]+alpha[bedcondi[0][0]-1,i])*dzbed_sur_dt[i]*self.Sbh
                    intS_PindS[i] = -0.5*(aUa[1, bedcondi[0][0],i]+aUa[1, bedcondi[0][0]-1,i])*self.Sbh
            intStop_PindS[i] = intStop_PindS[i] + aUa[1, -1, i]*self.Sbh

            hf = h5py.File(self.solsav+'/read_massBalance'+str(self.asint).replace('.','p')+'.h5', 'w')
            hf.create_dataset('intdVolPhi_sur_dt',data=intdVolPhi_sur_dt)
            hf.create_dataset('phi_dzbed_sur_dt',data=phi_dzbed_sur_dt)
            hf.create_dataset('intS_PindS',data=intS_PindS)
            hf.create_dataset('intStop_PindS',data=intStop_PindS)
            hf.close()
        
        else :
            for t in tqdm(self.tread):
                k = k + 1
                self.time[k] = float(t)
            hf = h5py.File(self.solsav+'/read_massBalance'+str(self.asint).replace('.','p')+'.h5', 'r')
            intdVolPhi_sur_dt = np.array(hf.get('intdVolPhi_sur_dt'))
            phi_dzbed_sur_dt = np.array(hf.get('phi_dzbed_sur_dt'))
            intS_PindS = np.array(hf.get('intS_PindS'))
            intStop_PindS = np.array(hf.get('intStop_PindS'))
        
        # Forward [:-1] Backward [1:] for temporal derivative        
        self.phi_dzbed_sur_dt = phi_dzbed_sur_dt[:-1]           # term 1 : bed evolution
        self.intdVolPhi_sur_dt = intdVolPhi_sur_dt[:-1]         # term 2 : storage evolution
        self.intS_PindS = intS_PindS[:-1]                       # term 4 : bed flux
        self.intStop_PindS = intStop_PindS[:-1]                 # term 5 : top flux

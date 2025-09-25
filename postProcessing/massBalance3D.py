#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Filename: massBalance3D.py
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

from func1 import create_point_2Dcyl, save_point, read_point, uns2cartvec, uns2cart, ddt

class ThreedimSqrsimu(object):
    def __init__(self, basepath, case, g0X, g2X, N0X, N1X, N2X, g0Y, g2Y, N0Y, N1Y, N2Y, gZ, gZ1, nzb, nz1, read_point_nc, save_point_nc, read_mesh, read_result, asint):
        self.basepath = basepath
        self.case = case
        self.sol = self.basepath + self.case + "/"
        self.solsav = self.basepath + self.case + "/constant"
        self.asint = asint      # Solid volume fraction criteria
        self.D = 0.1        
        self.xmin = -0.2        # upstream limit of the control volume
        self.xmax = -0.05       # downstream limit of the control volume
        self.ymin = -0.05       # right limit of the control volume
        self.ymax = 0.05        # left limit of the control volume
        self.read_point_nc = read_point_nc
        self.save_point_nc = save_point_nc
        
        self.createmesh(g0X, g2X, N0X, N1X, N2X, g0Y, g2Y, N0Y, N1Y, N2Y, gZ, gZ1, nzb, nz1, read_mesh, read_result)
        
    def createmesh(self, g0X, g2X, N0X, N1X, N2X, g0Y, g2Y, N0Y, N1Y, N2Y, gZ, gZ1, nzb, nz1, read_mesh, read_result):
        if not read_mesh:
            Xb, Yb, Zb = fluidfoam.readmesh(self.sol, precision = 13)
            hf = h5py.File(self.solsav+'/read_mesh.h5', 'w')
            hf.create_dataset('Xb',data=Xb)
            hf.create_dataset('Yb',data=Yb)
            hf.create_dataset('Zb',data=Zb)
            hf.close()
        else :
            hf = h5py.File(self.solsav+'/read_mesh.h5', 'r')
            Xb = np.array(hf.get('Xb'))
            Yb = np.array(hf.get('Yb'))
            Zb = np.array(hf.get('Zb'))
            
        ncell = np.size(Xb)
        if not self.read_point_nc:
            n2d, nz, pbed, pointer = create_point_2Dcyl(Xb, Yb, Zb)
            if self.save_point_nc:
                save_point(self.solsav, n2d, nz, pbed, Xb, Yb, Zb, pointer)
        else:
            n2d, nz, pbed, pointer = read_point(self.solsav)

        Xb = Xb[pbed]
        Yb = Yb[pbed]
        Zb = Zb[pointer][0,:]
        nx = np.unique(np.round(Xb,6)).shape[0]
        ny = np.unique(np.round(Yb,6)).shape[0]
        nz = np.unique(np.round(Zb,6)).shape[0]
        nvoid = nx*ny*nz - ncell
        indcart = np.zeros((nx,ny))
        xcart = np.unique(np.round(Xb,6))
        ycart = np.unique(np.round(Yb,6))
        indorig = range(n2d)

        i=-1
        j=-1
        for xi in tqdm(xcart):
            j=-1
            i=i+1
            for yi in ycart:
                j=j+1
                if np.logical_and(np.abs(xi) <= 0.05, np.abs(yi) <= 0.05):
                    indcart[i,j] = -1
                else:
                    indcart[i,j]=int(np.argwhere(np.logical_and(np.round(Xb,6)==xi,np.round(Yb,6)==yi))[0][0])
        indcart=indcart.astype(int)
        
        imin = np.where(xcart >= self.xmin)[0][0]        # upstream limit of the control volume
        imax = np.where(xcart <= self.xmax)[-1][-1]      # downstream limit of the control volume

        jmin = np.where(ycart >= self.ymin)[0][0]        # right limit of the control volume
        jmax = np.where(ycart <= self.ymax)[-1][-1]      # left limit of the control volume
        
        #
        # Construct longitudinal mesh
        #
        
        #Length
        L0X = 0.35  # Upstream
        L1X = 0.1   # Cylinder
        L2X = 0.75  # Downstream

        dx2, dx2N = fluidfoam.meshdesign.getdzs(L2X, g2X, N2X)
        x2i, dx2, g2X = fluidfoam.meshdesign.getgz(L2X, dx2, N2X)
        dx1 = np.ones(N1X)*L1X/N1X
        x1i = np.linspace(-0.05,0.05,N1X+1)
        dx0, dx0N = fluidfoam.meshdesign.getdzs(L0X, g0X, N0X)
        x0i, dx0, g0X = fluidfoam.meshdesign.getgz(L0X, dx0, N0X)

        deltaX = np.concatenate((dx0, dx1, dx2))
        xp = np.concatenate((x0i-0.4, x1i[1:], x2i[1:]+0.05))

        #
        # Construct transversal mesh
        #
        
        # Length
        L0Y = 0.35  # Right
        L1Y = 0.1   # Cylinder
        L2Y = 0.35  # Left

        dy2, dy2N = fluidfoam.meshdesign.getdzs(L2Y, g2Y, N2Y)
        y2i, dy2, g2Y = fluidfoam.meshdesign.getgz(L2Y, dy2, N2Y)
        dy1 = np.ones(N1Y)*L1Y/N1Y
        y1i = np.linspace(-0.05,0.05,N1Y+1)
        dy0, dy0N = fluidfoam.meshdesign.getdzs(L0Y, g0Y, N0Y)
        y0i, dy0, g0Y = fluidfoam.meshdesign.getgz(L0Y, dy0, N0Y)

        deltaY = np.concatenate((dy0, dy1, dy2))
        yp = np.concatenate((y0i-0.4, y1i[1:], y2i[1:]+0.05))

        #
        # Construct vertical mesh
        #

        # Length
        Lzb = 0.2   # Lower part (water)
        Lz1 = 0.05  # Part (sediment)

        dzb, dzbN = fluidfoam.meshdesign.getdzs(Lzb, gZ, nzb)
        zbi, dzb, gZ = fluidfoam.meshdesign.getgz(Lzb, dzb, nzb)
        dz1, dz1N = fluidfoam.meshdesign.getdzs(Lz1, gZ1, nz1)
        z1i, dz1, gZ1 = fluidfoam.meshdesign.getgz(Lz1, dz1, nz1)

        deltaZ = np.concatenate((dz1, dzb))
        zp= np.concatenate((z1i-Lz1, zbi[1:]))
        
        self.Xb = Xb
        self.Yb = Yb
        self.Zb = Zb
        self.imin = imin
        self.imax = imax
        self.jmin = jmin
        self.jmax = jmax
        self.xcart = xcart
        self.ycart = ycart
        self.ncell = ncell
        self.n2d = n2d
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.pbed = pbed
        self.pointer = pointer
        self.deltaX = deltaX
        self.deltaY = deltaY
        self.indorig = indorig
        self.indcart = indcart
        self.deltaZ = deltaZ
        self.zp = zp
        
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
        if len(tread) != 1 :
            tread= tread[:80]
            dt = float(tread[1])-float(tread[0])
            Nt = len(tread)
        else :
            tread = np.arange(0.25,20.25,0.25)
            tread = np.array(tread).astype(str)
            dt = 0.25
            Nt = 80
        k=-1
        for s in tread:
            k=k+1
            s=s.replace(".0","")
            tread[k]=s
            
        time = np.zeros(Nt)
        
        self.time = time
        self.dt = dt
        self.Nt = Nt
        self.tread = tread
        
        self.postProcess(read_result)
          
    def postProcess(self, read_result):
        k = -1
        if not read_result:
            dVuns = np.zeros(self.ncell)
            dV = np.zeros((self.nx, self.ny, self.nz))
            print(np.shape(self.xcart), np.shape(self.ycart), np.shape(dV))
            alphauns = np.zeros((self.ncell, self.Nt))
            alpha = np.zeros((self.nx, self.ny, self.nz, self.Nt))
            aUauns = np.zeros((3, self.ncell, self.Nt))
            aUa = np.zeros((3, self.nx, self.ny, self.nz, self.Nt))

            zbed, dzbed_sur_dt, dVolPhi_sur_dt, VoldPhi_sur_dt = np.zeros((4, self.nx, self.ny, self.Nt))
            phi_dzbed_sur_dt, intdVolPhi_sur_dt, intVolPhi = np.zeros((3, self.Nt))
            intS_PindS, intS2_PindS, intS3_PindS, intSlat_PindS, intSlatin_PindS, intSlatleft_PindS, intSlatright_PindS, intSlatout_PindS, intStop_PindS = np.zeros((9, self.Nt))

            dVuns[:] = fluidfoam.readscalar(self.sol, "constant", "V", verbose=False, precision=13)
            varorig = dVuns[:]
            varorig = varorig[self.pointer]
            dV[:,:,:] = uns2cart(varorig, self.nx, self.ny, self.nz, self.indcart, self.indorig)
            print(np.sum(varorig), np.nansum(dV))
            
            for t in tqdm(self.tread):
                print("Reading time: %s s" % t)
                k = k + 1
                alphauns[:, k] = fluidfoam.readscalar(self.sol, t, "alpha.a", verbose=False, precision=13)
                aUauns[:,:,k] = fluidfoam.readvector(self.sol, t, "alphaUa", verbose=False, precision=13)
                self.time[k] = float(t)
                varorig = alphauns[:,k]
                varorig = varorig[self.pointer]
                alpha[:,:,:,k] = uns2cart(varorig, self.nx, self.ny, self.nz, self.indcart, self.indorig)
                varorig = aUauns[:,:,k]
                varorig = varorig[:,self.pointer]
                print(np.shape(varorig))
                aUa[:,:,:,:,k] = uns2cartvec(varorig, self.nx, self.ny, self.nz, self.indcart, self.indorig)
            
            for i in range(self.imin,self.imax+1):
                for j in range(self.jmin,self.jmax+1):
                    for t in range(self.Nt):
                        if np.isnan(alpha[i,j,:,t]).any():
                            pass
                        else:
                            bedcondi = np.where(alpha[i,j,:,t] <= self.asint)
                            zbed[i, j, t] = self.zp[bedcondi[0][0]]
                            intVolPhi[t] = intVolPhi[t] + np.nansum(alpha[i,j,bedcondi,t]*dV[i,j,bedcondi])
                    dzbed_sur_dt[i, j, :] = ddt('eulerFor', zbed[i, j, :], self.dt)
            
            intdVolPhi_sur_dt = ddt('eulerFor', intVolPhi, self.dt)

            for i in range(self.imin,self.imax+1):
                for j in range(self.jmin,self.jmax+1):
                    for t in range(self.Nt):
                        if np.isnan(alpha[i,j,:,t]).any():
                            pass
                        else:
                            bedcondi = np.where(alpha[i,j,:,t] <= self.asint)
                            bedcondi1 = np.where(alpha[i+1,j,:,t] <= self.asint)
                            if bedcondi[0][0] == 0 :
                                phi_dzbed_sur_dt[t] = phi_dzbed_sur_dt[t] + alpha[i, j, bedcondi[0][0], t]*dzbed_sur_dt[i, j, t]*self.deltaY[j]*self.deltaX[i]
                            else:
                                phi_dzbed_sur_dt[t] = phi_dzbed_sur_dt[t] + (alpha[i, j, bedcondi[0][0], t]*self.deltaZ[bedcondi[0][0]-1]+alpha[i, j, bedcondi[0][0]-1, t]*self.deltaZ[bedcondi[0][0]])/(self.deltaZ[bedcondi[0][0]]+self.deltaZ[bedcondi[0][0]-1])*dzbed_sur_dt[i, j, t]*self.deltaY[j]*self.deltaX[i]
                                intS_PindS[t] = intS_PindS[t] - (aUa[2, i, j, bedcondi[0][0], t]*self.deltaZ[bedcondi[0][0]-1]+aUa[2, i, j, bedcondi[0][0]-1, t]*self.deltaZ[bedcondi[0][0]])/(self.deltaZ[bedcondi[0][0]]+self.deltaZ[bedcondi[0][0]-1])*self.deltaY[j]*self.deltaX[i]
                            if i!=self.imax:
                                dz = (zbed[i+1, j, t] - zbed[i, j, t])
                                match dz:
                                    case dz if dz > 0:
                                        S2 = np.where(np.logical_and(alpha[i,j,:,t] <= self.asint, alpha[i+1,j,:,t] > self.asint))[0]
                                        intS2_PindS[t] = intS2_PindS[t] + np.nansum((aUa[0, i, j, S2, t]*self.deltaX[i+1]+aUa[0, i+1, j, S2, t]*self.deltaX[i])/(self.deltaX[i]+self.deltaX[i+1])*self.deltaY[j]*self.deltaZ[S2])
                                    case dz if dz < 0:
                                        S2 = np.where(np.logical_and(alpha[i,j,:,t] > self.asint, alpha[i+1,j,:,t] <= self.asint))[0]
                                        intS2_PindS[t] = intS2_PindS[t] - np.nansum((aUa[0, i, j, S2, t]*self.deltaX[i+1]+aUa[0, i+1, j, S2, t]*self.deltaX[i])/(self.deltaX[i]+self.deltaX[i+1])*self.deltaY[j]*self.deltaZ[S2])
                        intStop_PindS[t] = intStop_PindS[t] + aUa[2, i, j, -1, t]*self.deltaY[j]*self.deltaX[i]

            for i in range(self.imin,self.imax+1):
                for j in range(self.jmin,self.jmax+1):
                    for t in range(self.Nt):
                        if np.isnan(alpha[i,j,:,t]).any():
                            pass
                        else:
                            bedcondi = np.where(alpha[i,j,:,t] <= self.asint)
                            bedcondi1 = np.where(alpha[i,j+1,:,t] <= self.asint)
                            if j!=self.jmax:
                                dz = (zbed[i, j+1, t] - zbed[i, j, t])
                                match dz:
                                    case dz if dz > 0:
                                        S3 = np.where(np.logical_and(alpha[i,j,:,t] <= self.asint, alpha[i,j+1,:,t] > self.asint))[0]
                                        intS3_PindS[t] = intS3_PindS[t] + np.nansum((aUa[1, i, j, S3, t]*self.deltaY[j+1]+aUa[1, i, j+1, S3, t]*self.deltaY[j])/(self.deltaY[j]+self.deltaY[j+1])*self.deltaX[i]*self.deltaZ[S3])
                                    case dz if dz < 0:
                                        S3 = np.where(np.logical_and(alpha[i,j,:,t] > self.asint, alpha[i,j+1,:,t] <= self.asint))[0]
                                        intS3_PindS[t] = intS3_PindS[t] - np.nansum((aUa[1, i, j, S3, t]*self.deltaY[j+1]+aUa[1, i, j+1, S3, t]*self.deltaY[j])/(self.deltaY[j]+self.deltaY[j+1])*self.deltaX[i]*self.deltaZ[S3])
                        
            for j in range(self.jmin,self.jmax+1):
                for t in range(self.Nt):
                    if np.isnan(alpha[self.imin,j,:,t]).any():
                        pass
                    else:
                        bedcondi = np.where(alpha[self.imin,j,:,t] <= self.asint)
                        for k in range(bedcondi[0][0],self.nz):
                            intSlatin_PindS[t] = intSlatin_PindS[t] - (aUa[0, self.imin, j, k, t]*self.deltaX[self.imin-1]+aUa[0, self.imin-1, j, k, t]*self.deltaX[self.imin])/(self.deltaX[self.imin]+self.deltaX[self.imin-1])*self.deltaY[j]*self.deltaZ[k]

            for j in range(self.jmin,self.jmax+1):
                for t in range(self.Nt):
                    if np.isnan(alpha[self.imax,j,:,t]).any():
                        pass
                    else:
                        bedcondi = np.where(alpha[self.imax,j,:,t] <= self.asint)
                        for k in range(bedcondi[0][0],self.nz):
                            intSlatout_PindS[t] = intSlatout_PindS[t] + (aUa[0, self.imax, j, k, t]*self.deltaX[self.imax+1] + aUa[0, self.imax+1, j, k, t]*self.deltaX[self.imax])/(self.deltaX[self.imax]+self.deltaX[self.imax+1])*self.deltaY[j]*self.deltaZ[k]

            for i in range(self.imin,self.imax+1):
                for t in range(self.Nt):
                    if np.isnan(alpha[i,self.jmin,:,t]).any():
                        pass
                    else:
                        bedcondi = np.where(alpha[i,self.jmin,:,t] <= self.asint)
                        for k in range(bedcondi[0][0],self.nz):
                            intSlatright_PindS[t] = intSlatright_PindS[t] - (aUa[1, i, self.jmin, k, t]*self.deltaY[self.jmin-1]+aUa[1, i, self.jmin-1, k, t]*self.deltaY[self.jmin])/(self.deltaY[self.jmin]+self.deltaY[self.jmin-1])*self.deltaX[i]*self.deltaZ[k]
                        
            for i in range(self.imin,self.imax+1):
                for t in range(self.Nt):
                    if np.isnan(alpha[i,self.jmax,:,t]).any():
                        pass
                    else:
                        bedcondi = np.where(alpha[i,self.jmax,:,t] <= self.asint)
                        for k in range(bedcondi[0][0],self.nz):
                            intSlatleft_PindS[t] = intSlatleft_PindS[t] + (aUa[1, i, self.jmax, k, t]*self.deltaY[self.jmax+1]+aUa[1, i, self.jmax+1, k, t]*self.deltaY[self.jmax])/(self.deltaY[self.jmax]+self.deltaY[self.jmax+1])*self.deltaX[i]*self.deltaZ[k]

            for t in range(self.Nt):
                intSlat_PindS[t] = np.nansum([intSlatin_PindS[t],intSlatout_PindS[t],intSlatright_PindS[t],intSlatleft_PindS[t]])
            
            
            hf = h5py.File(self.solsav+'/read_massBalance'+str(self.asint).replace('.','p')+'.h5', 'w')
            hf.create_dataset('phi_dzbed_sur_dt',data=phi_dzbed_sur_dt)
            hf.create_dataset('intdVolPhi_sur_dt',data=intdVolPhi_sur_dt)
            hf.create_dataset('intS_PindS',data=intS_PindS)
            hf.create_dataset('intS2_PindS',data=intS2_PindS)
            hf.create_dataset('intS3_PindS',data=intS3_PindS)
            hf.create_dataset('intStop_PindS',data=intStop_PindS)
            hf.create_dataset('intSlat_PindS',data=intSlat_PindS)
            hf.create_dataset('intSlatleft_PindS',data=intSlatleft_PindS)
            hf.create_dataset('intSlatright_PindS',data=intSlatright_PindS)
            hf.create_dataset('intSlatin_PindS',data=intSlatin_PindS)
            hf.create_dataset('intSlatout_PindS',data=intSlatout_PindS)
            hf.close()
        else :
            for t in tqdm(self.tread):
                k = k + 1
                self.time[k] = float(t)
            with h5py.File(self.solsav+'/read_massBalance'+str(self.asint).replace('.','p')+'.h5', 'r') as hf:
                intVoldPhi_sur_dt = np.array(hf.get('intVoldPhi_sur_dt'))
                phi_dzbed_sur_dt = np.array(hf.get('phi_dzbed_sur_dt'))
                intdVolPhi_sur_dt = np.array(hf.get('intdVolPhi_sur_dt'))
                intS_PindS = np.array(hf.get('intS_PindS'))
                intS2_PindS = np.array(hf.get('intS2_PindS'))
                intS3_PindS = np.array(hf.get('intS3_PindS'))
                intStop_PindS = np.array(hf.get('intStop_PindS'))
                intSlat_PindS = np.array(hf.get('intSlat_PindS'))
                intSlatleft_PindS = np.array(hf.get('intSlatleft_PindS'))
                intSlatright_PindS = np.array(hf.get('intSlatright_PindS'))
                intSlatin_PindS = np.array(hf.get('intSlatin_PindS'))
                intSlatout_PindS = np.array(hf.get('intSlatout_PindS'))

        # Forward [:-1] Backward [1:] for temporal derivative
        self.phi_dzbed_sur_dt = phi_dzbed_sur_dt[:-1]           # term 1 : bed evolution
        self.intdVolPhi_sur_dt = intdVolPhi_sur_dt[:-1]         # term 2 : storage evolution
        self.intSlat_PindS = intSlat_PindS[:-1]                 # term 3 : lateral flux
        self.intSlatleft_PindS = intSlatleft_PindS[:-1]         # term 3a : lateral flux -> left face
        self.intSlatright_PindS = intSlatright_PindS[:-1]       # term 3b : lateral flux -> right face
        self.intSlatin_PindS = intSlatin_PindS[:-1]             # term 3c : lateral flux -> upstream face
        self.intSlatout_PindS = intSlatout_PindS[:-1]           # term 3d : lateral flux -> downstream face
        self.intS_PindS = intS_PindS[:-1]                       # term 4a : bed flux -> vertical component
        self.intS2_PindS = intS2_PindS[:-1]                     # term 4b : bed flux -> streamwise component
        self.intS3_PindS = intS3_PindS[:-1]                     # term 4c : bed flux -> spanwise component
        self.intStop_PindS = intStop_PindS[:-1]                 # term 5 : top flux

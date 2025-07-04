#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Filename: func1.py
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

def divergence(f,x):
    import numpy as np
    num_dims = len(f)
    return np.ufunc.reduce(np.add, [np.gradient(f[i], x[i], axis=i) for i in range(num_dims)])

def ddt(arg, f, dt):
    import numpy as np
    if arg == 'eulerBack':
        df_dt = np.zeros(len(f))
        for i in range(len(f)-1):
            df_dt[i+1] = (f[i+1]-f[i])/dt
        return df_dt
    if arg == 'eulerFor':
        df_dt = np.zeros(len(f))
        for i in range(len(f)-1):
            df_dt[i] = (f[i+1]-f[i])/dt
        return df_dt
    if arg == 'trapBack':
        df_dt = np.zeros(len(f))
        for i in range(len(f)-1):
            df_dt[i+1] = ((f[i+1]+f[i])/2-f[i])/(dt/2)
        return df_dt
    if arg == 'trapFor':
        df_dt = np.zeros(len(f))
        for i in range(len(f)-1):
            df_dt[i] = ((f[i+1]+f[i])/2-f[i])/(dt/2)
        return df_dt
    if arg == 'gradient':
        return np.gradient(f,dt)
    
def create_point_2Dcyl(X, Y, Z):
    import numpy as np
    from tqdm import tqdm

    Ncell = np.size(X)
    #
    #---create a list of (x,y) positions correpsonding to Zb = min(Zb)
    #
    print("Generate pointer list")
    # BONAMY TU EN LOUPES PEUT ETRE A CAUSE DE LA PRECISION
    pbed = np.where(Z == np.amin(Z))[0]

    n2d = np.size(pbed)
    nz = 0
    profileList = []
    for i in tqdm(range(n2d)):
        #Extract the list of points having the same (Xb, Yb) values
        indices = np.where(np.logical_and(np.round(X,10) == np.round(X[pbed[i]],10),
                                          np.round(Y,10) == np.round(Y[pbed[i]],10)))[0]
        nz = max(nz, np.size(indices))
        profile = list(indices)
        #Sort the list of points by increasing Zb values (profile)
        profileSorted = [x for _, x in sorted(zip(Z[profile], profile))]
        #Append the list to the profileList
        profileList.append(profileSorted)

    # Convert profileList to numpy.array
    pointer = np.array(profileList)
    return n2d, nz, pbed, pointer
    
def create_point_1Dcart(X, Y):
    import numpy as np
    from tqdm import tqdm

    Ncell = np.size(X)
    #
    #---create a list of (x) positions correpsonding to Zb = min(Zb)
    #
    print("Generate pointer list")
    # BONAMY TU EN LOUPES PEUT ETRE A CAUSE DE LA PRECISION
    pbed = np.where(Y == np.amin(Y))[0]

    n1d = np.size(pbed)
    ny = 0
    profileList = []
    for i in tqdm(range(n1d)):
        #Extract the list of points having the same (Xb) values
        indices = np.where(X == X[pbed[i]])[0]
        ny = max(ny, np.size(indices))
        profile = list(indices)
        #Sort the list of points by increasing Yb values (profile)
        profileSorted = [x for _, x in sorted(zip(Y[profile], profile))]
        #Append the list to the profileList
        profileList.append(profileSorted)

    # Convert profileList to numpy.array
    pointer = np.array(profileList)
    return n1d, ny, pbed, pointer

def save_1Dpoint(solsav, n1d, ny, pbed, X, Y, pointer):
    import os
    import numpy as np
    from netCDF4 import Dataset

    ###############################################################################
    # save NetCDF files in constant and timeStep folders
    # --------------------------------------
    #
    postProcFile = os.path.join(solsav, 'pointerpostproc.nc')
    print ("save postprocFile in:",postProcFile)

    # NetCDF file creation
    rootgrp = Dataset(postProcFile, 'w')

    # Dimensions creation
    rootgrp.createDimension('X', n1d)
    rootgrp.createDimension('Y', ny)

    # Variables creation
    pb_file = rootgrp.createVariable('pbed', np.int64, 'X')
    # x_file  = rootgrp.createVariable('x', np.float64, 'X')
    # y_file  = rootgrp.createVariable('y', np.float64, ('X','Y'))
    p_file  = rootgrp.createVariable('p', np.int64, ('X','Y'))

    # Writing variables
    pb_file[:]  = pbed
#    x_file[:]   = X[pbed]
#    y_file[:]   = Y[pbed]
    p_file[:,:] = pointer[:,:]

    # File closing
    rootgrp.close()

def read_1Dpoint(solsav):
    import os
    import numpy as np
    from netCDF4 import Dataset

    ###############################################################################
    # read NetCDF files in constant and timeStep folders
    # --------------------------------------
    #
    postProcFile = os.path.join(solsav, 'pointerpostproc.nc')
    print ("Read postprocFile in:",postProcFile)

    # NetCDF file read
    pointer_file = Dataset(postProcFile)
    pbed = pointer_file.variables['pbed'][:]
    pointer = pointer_file.variables['p'][:,:]
    n1d = np.size(pbed)
    ny = np.shape(pointer)[1]
    return n1d, ny, pbed, pointer

def save_point(solsav, n2d, nz, pbed, X, Y, Z, pointer):
    import os
    import numpy as np
    from netCDF4 import Dataset

    ###############################################################################
    # save NetCDF files in constant and timeStep folders
    # --------------------------------------
    #
    postProcFile = os.path.join(solsav, 'pointerpostproc.nc')
    print ("save postprocFile in:",postProcFile)

    # NetCDF file creation
    rootgrp = Dataset(postProcFile, 'w')

    # Dimensions creation
    rootgrp.createDimension('XY', n2d)
    rootgrp.createDimension('Z', nz)

    # Variables creation
    pb_file = rootgrp.createVariable('pbed', np.int64, 'XY')
#    x_file  = rootgrp.createVariable('x', np.float64, 'XY')
#    y_file  = rootgrp.createVariable('y', np.float64, 'XY')
#    z_file  = rootgrp.createVariable('z', np.float64, ('XY','Z'))
    p_file  = rootgrp.createVariable('p', np.int64, ('XY','Z'))

    # Writing variables
    pb_file[:]  = pbed
#    x_file[:]   = X[pbed]
#    y_file[:]   = Y[pbed]
#    z_file[:,:] = Z[pointer[:,:]]
    p_file[:,:] = pointer[:,:]

    # File closing
    rootgrp.close()


def read_point(solsav):
    import os
    import numpy as np
    from netCDF4 import Dataset

    ###############################################################################
    # read NetCDF files in constant and timeStep folders
    # --------------------------------------
    #
    postProcFile = os.path.join(solsav, 'pointerpostproc.nc')
    print ("Read postprocFile in:",postProcFile)

    # NetCDF file read
    pointer_file = Dataset(postProcFile)
    pbed = pointer_file.variables['pbed'][:]
    pointer = pointer_file.variables['p'][:,:]
    n2d = np.size(pbed)
    nz = np.shape(pointer)[1]
    return n2d, nz, pbed, pointer


def uns2cylvec(array, nr, ntheta, nz, indcyl, indcart):
    import numpy as np
    varcyl = np.zeros((3, nr, ntheta, nz))
    for l in range(nr):
        for m in range(ntheta):
            ind=indcyl[l,m]
            #varcyl[:,l,m,:]=array[:,np.where(ind==indcart)[0][0],:]
            varcyl[:,l,m,:]=array[:,ind,:]
    return varcyl


def uns2cyl(array, nr, ntheta, nz, indcyl, indcart):
    import numpy as np
    varcyl = np.zeros((nr, ntheta, nz))
    for l in range(nr):
        for m in range(ntheta):
            ind=indcyl[l,m]
            #varcyl[l,m,:]=array[np.where(ind==indcart)[0][0],:]
            varcyl[l,m,:]=array[ind,:]
    return varcyl

def uns2cart(array, nx, ny, nz, indcart, indorig):
    import numpy as np
    varcart = np.zeros((nx, ny, nz))
    for xi in range(nx):
        for yi in range(ny):
            ind=indcart[xi,yi]
            if ind < 0:
                varcart[xi,yi,:]= np.nan
            else:
                #varcyl[l,m,:]=array[np.where(ind==indcart)[0][0],:]
                varcart[xi,yi,:]=array[ind,:]
    return varcart

def uns2cartvec(array, nx, ny, nz, indcart, indorig):
    import numpy as np
    varcart = np.zeros((3, nx, ny, nz))
    for xi in range(nx):
        for yi in range(ny):
            ind=indcart[xi,yi]
            if ind < 0:
                varcart[:,xi,yi,:]= np.nan
            else:
                #varcyl[:,l,m,:]=array[:,np.where(ind==indcart)[0][0],:]
                varcart[:,xi,yi,:]=array[:,ind,:]
    return varcart

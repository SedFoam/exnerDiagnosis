# README #

This repository provides the python scripts reproducing the figures of the article for sedFOAM solver (version python 3.12.0 and version OpenFOAM classic). You may need to install two python package: tqdm, h5py. You may simply use "pip install --user tqdm" and "pip install --user h5py"

## Name

Exner diagnosis method for two-fluid morphodynamics simulations

## Description

The objective is to reproduce the time evolution of each term part the Exner equation for the two-phase flow model.

## Usage

Desription of each python scripts:

Python scripts creating class and produce time evolution of each term :

massBalance1D.py -> 1D sedimentation configuration
massBalance2D.py -> 2D scour downstream an apron configuration
massBalance3D.py -> 3D scour around a square cylinder configuration

Name of each term : 

phi_dzbed_sur_dt   -> term 1  : bed evolution
intdVolPhi_sur_dt  -> term 2  : storage evolution
intSlat_PindS      -> term 3  : lateral flux
intS_PindS         -> term 4a : bed flux vertical component   ( available for 1D, 2D, 3D configurations )
intS2_PindS        -> term 4b : bed flux streamwise component ( available for 2D, 3D configurations     )
intS3_PindS        -> term 4c : bed flux spanwise component   ( available for 3D configuration          )
intStop_PindS      -> term 5 : top flux                       ( available for 1D, 2D, 3D configurations )

Python scripts reproducing the figure of the article:

figures1D.py -> reproduce figures 6, 7a and 7b
figures2D.py -> reproduce figures 9a, 9b, 10 (result output at time 50 required), 12a and 12b
figures3D.py -> reproduce figures 14a, 14b, 16a and 16b

Results are already saved in Configuration folder, to reproduce the figures, simply run figures*D.py as it is coded. It is possible to reconstruct the data from scratch with the following steps:

1D configuration :

Step 1 : run Allrun and wait for the end of the simulation
Step 2 : switch to False the variable meshCreated l.24 in figures1D.py to produce the read_mesh.h5 file
Step 3 : switch to False the variable resultCreated l.25 in figures1D.py to produce the read_massBalance*.h5 files
Step 4 : run figures1D.py

2D configuration :

Step 1 : run Allrun.pre and wait for the end of the 1D simulation required for the inlet boundary condition
Step 2 : run Allrun and wait for the end of the 2D simulation
Step 3 : switch to False the variable meshCreated l.25 in figures2D.py to produce the read_mesh.h5 file
Step 4 : switch to False the variable resultCreated l.26 in figures2D.py to produce the read_massBalance*.h5 files
Step 5 : switch to False the variable read_point_nc l.30 in figures2D.py to produce pointerpostproc.nc file
Step 6 : switch to True the variable save_point_nc l.31 in figures2D.py to save pointerpostproc.nc file
Step 7 : run figures2D.py

3D configuration :

Step 1 : run Allrun.pre and wait for the end of the 1D simulation required for the inlet boundary condition
Step 2 : run Allrun and wait for the end of the 3D simulation
Step 3 : switch to False the variable meshCreated l.40 in figures3D.py to produce the read_mesh.h5 file
Step 4 : switch to False the variable resultCreated l.41 in figures3D.py to produce the read_massBalance*.h5 files
Step 5 : switch to False the variable read_point_nc l.38 in figures3D.py to produce pointerpostproc.nc file
Step 6 : switch to True the variable save_point_nc l.39 in figures3D.py to save pointerpostproc.nc file
Step 7 : run figures3D.py

It is also possible to construct new results in a different zone by modifying the massBalance*D.py file described as follows:

1D configuration :

Only 1 zone possible, the one above the bed.

2D configuration :

Step 1 : modify the upstream limit with the variable xmin l.21 in massBalance2D.py
Step 2 : modify the downstream limit with the variable xmin l.22 in massBalance2D.py
Step 3 : switch to False the variable resultCreated l.26 in figures2D.py to produce the new read_massBalance*.h5 files
Step 4 : run figures2D.py (no need to do steps 3, 5 and 6 previously described if read_mesh.h5 and pointerpostproc.nc already exist)

3D configuration :

Step 1 : modify the upstream limit with the variable xmin l.22 in massBalance3D.py
Step 2 : modify the downstream limit with the variable xmin l.23 in massBalance3D.py
Step 3 : modify the right limit with the variable ymin l.24 in massBalance3D.py
Step 4 : modify the left limit with the variable ymin l.25 in massBalance3D.py
Step 5 : switch to False the variable resultCreated l.41 in figures3D.py to produce the new read_massBalance*.h5 files
Step 6 : run figures3D.py (no need to do steps 3, 5 and 6 previously described if read_mesh.h5 and pointerpostproc.nc already exist)

## Authors

Cyrille Bonamy (LEGI UMR 5519 - Grenoble - France)
Cyrille.Bonamy (A) univ-grenoble-alpes.fr
Julien Chauchat (LEGI UMR 5519 - Grenoble - France)
Julien.Chauchat (A) univ-grenoble-alpes.fr
Alban Gilletta (LEGI UMR 5519 - Grenoble - France)
Alban.Gilletta-de-Saint-Joseph (A) univ-grenoble-alpes.fr

## License
                    GNU GENERAL PUBLIC LICENSE
                       Version 3, 29 June 2007

 Copyright (C) 2007 Free Software Foundation, Inc. <https://fsf.org/>
 Everyone is permitted to copy and distribute verbatim copies
 of this license document, but changing it is not allowed.

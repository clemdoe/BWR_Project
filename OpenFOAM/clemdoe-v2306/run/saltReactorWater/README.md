# OpenFoam v2306 Models for THMSAL validation.

* Author: Cristian Garrido Tamm
* e-mail: cristian.garrido@idom.com
* Version: v0
* Date: 03/03/2024

## PREREQUISITES

To generate the mesh:
* gmsh Python API 

To run the simulation:
* OpenFoam v2306
* PyFoam (Recommended)

## MODEL DESCRIPTION

### GEOMETRY

This model represent half single pin with three regions (fuel, cladding and coolant) in cartesian geometry. The mesh for the 3 regions is generated using the Python script cartesian_3region.py inside the gmsh folder. At the top of the script there is parameters to modify the side, height and radius as well as to control the mesh resolution. This script will generate the file cartesian_3region.msh. This file (or a symbolic link to it) needs to be located in a gmsh folder in the case directory. The Allrun* scripts will then use the gmshToFoam utility to convert the mesh to OpenFoam format. 
For convinience, two meshes for 1 cm and 100 cm height rod models are provided in the gmsh folder.

### BOUNDARY CONDITIONS

Boundary conditions and initial fields are set in the changeDictionaryDict files in the system/<region> folders. The boundary conditions and initial fields in the 0.orig directory are dummy and will be overwriten before running.

### HEAT SOURCE

The heat source is given in the constant/fuel/fvOptions file.

### THERMOPHYSICAL PROPERTIES

The thermophysical properties of the coolant (FLiBe) are obtained from Oak Ridge's Molten Salt Thermophysical Properties Database, for cladding are extracted from the correlations used in DONJON source code. Polynomial fits are used when the correlation is not polynomial.

## RUNNING THE SIMULATION

The Allrun script will run the simulation using a single cpu with pyFoam.
The Allrun-parallell script will run the simulation using 40 cpus using mpi (for calculation clusters).

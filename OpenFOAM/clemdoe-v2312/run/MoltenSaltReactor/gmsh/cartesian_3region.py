#!/usr/bin/env python
# -*- coding: utf-8 -*-
# /home/cgt/OpenFOAM/cgt-v2306/gmsh/cartesian_3region.py
# Based on tutorial_gmsh.py (https://jsdokken.com/src/tutorial_gmsh.html)
import logging
import colorlog
from docopt import docopt
import sys
import gmsh

__AUTHOR__ = 'Cristian Garrido Tamm'
__EMAIL__ = 'cristian.garrido@idom.com'
__CREATION_DATE__ = 'Sunday Jul 09, 2023 16:25:04 CEST'
__VERSION__ = '1.0.1'
__LAST_MODIFIED__ = 'Sun 03 Mar 2024 05:04:34 PM CET'

# Definition of Geometrical parameters
fuelRadius = 0.455063*1E-2
innerCladRadius = 0.460169*1E-2
outerCladRadius = 0.500184*1E-2
side = 1.26*1E-2
#  coreHeight = 1.0E-2
coreHeight = 1.0

# Mesh size control
lc = side/15        # Maximum mesh size
factor = 0.95       # Factor defining refinment region arround cladding surface
div = 10            # Factor for minimum mesh size

def logger_setup():
    """ 
        Setup logger
    """
    # Set up logger
    global logger
    handler = colorlog.StreamHandler(sys.stdout)
    handler.setFormatter(colorlog.ColoredFormatter('%(log_color)s[ %(levelname)s ] %(log_color)s%(asctime)s -%(name)s- : %(message)s'))
    logger = logging.getLogger()
    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG)

def process(args):
    gmsh.initialize()
    gmsh.model.add("subChannel")

    # Points describing the problem
    p0 = gmsh.model.occ.addPoint(0, 0, 0, lc)
    p1 = gmsh.model.occ.addPoint(fuelRadius, 0, 0, lc)
    p2 = gmsh.model.occ.addPoint(outerCladRadius, 0, 0, lc)
    p3 = gmsh.model.occ.addPoint(side/2, 0, 0, lc)
    p4 = gmsh.model.occ.addPoint(side/2, side/2, 0, lc)
    p5 = gmsh.model.occ.addPoint(-side/2, side/2, 0, lc)
    p6 = gmsh.model.occ.addPoint(-side/2, 0, 0, lc)
    p7 = gmsh.model.occ.addPoint(-outerCladRadius, 0, 0, lc)
    p8 = gmsh.model.occ.addPoint(-fuelRadius, 0, 0, lc)

    # Add Lines
    l0_1 = gmsh.model.occ.addLine(p0, p1)
    l1_2 = gmsh.model.occ.addLine(p1, p2)
    l2_3 = gmsh.model.occ.addLine(p2, p3)
    l3_4 = gmsh.model.occ.addLine(p3, p4)
    l4_5 = gmsh.model.occ.addLine(p4, p5)
    l5_6 = gmsh.model.occ.addLine(p5, p6)
    l6_7 = gmsh.model.occ.addLine(p6, p7)
    l7_8 = gmsh.model.occ.addLine(p7, p8)
    l8_0 = gmsh.model.occ.addLine(p8, p0)

    # Add Arcs
    lfuel = gmsh.model.occ.addCircleArc(p8, p0, p1)
    lcladding = gmsh.model.occ.addCircleArc(p7, p0, p2)

    # Add a curve loop and a plane surface for the coolant region
    coolLoop = gmsh.model.occ.addCurveLoop([l2_3, l3_4, l4_5, l5_6, l6_7, lcladding])

    # Add a curve loop and a plane surface for the cladding region
    cladLoop = gmsh.model.occ.addCurveLoop([l1_2, -lcladding, l7_8, lfuel])
   
    # Add a curve loop and a plane surface for the cladding region
    fuelLoop = gmsh.model.occ.addCurveLoop([l0_1, -lfuel, l8_0])

    # Add regions
    fuelRegion = gmsh.model.occ.addPlaneSurface([fuelLoop])
    cladRegion = gmsh.model.occ.addPlaneSurface([cladLoop])
    coolRegion = gmsh.model.occ.addPlaneSurface([coolLoop])

    # Generate the cells by extrusion of the 2D regions
    coolCell = gmsh.model.occ.extrude([(2, coolRegion)], 0, 0, coreHeight, recombine=True)
    fuelCell = gmsh.model.occ.extrude([(2, fuelRegion)], 0, 0, coreHeight, recombine=True)
    cladCell = gmsh.model.occ.extrude([(2, cladRegion)], 0, 0, coreHeight, recombine=True)
    gmsh.model.occ.synchronize()
    
    # Generate Physical group for each material cell
    coolVolTag = [j for (i, j) in coolCell if i == 3]
    gmsh.model.addPhysicalGroup(3, coolVolTag, name='coolant')
    cladVolTag = [j for (i, j) in cladCell if i == 3] 
    gmsh.model.addPhysicalGroup(3, cladVolTag, name='cladding')
    fuelVolTag = [j for (i, j) in fuelCell if i == 3] 
    gmsh.model.addPhysicalGroup(3, fuelVolTag, name='fuel')

    # Set boundaries
    right = [coolCell[3][1]]
    back = [coolCell[4][1]]
    left = [coolCell[5][1]]
    front = [coolCell[2][1], coolCell[6][1], cladCell[2][1], cladCell[4][1], fuelCell[2][1], fuelCell[4][1]]
    inlet = [coolRegion]
    bottom = [fuelRegion, cladRegion]
    outlet = [coolCell[0][1]]
    top = [fuelCell[0][1], cladCell[0][1]]
    coolant_to_cladding = [coolCell[7][1]]
    cladding_to_coolant = [cladCell[3][1]]
    cladding_to_fuel = [cladCell[5][1]]
    fuel_to_cladding = [fuelCell[3][1]]

    gmsh.model.addPhysicalGroup(2, left, name='left')
    gmsh.model.addPhysicalGroup(2, front, name='front')
    gmsh.model.addPhysicalGroup(2, right, name='right')
    gmsh.model.addPhysicalGroup(2, back, name='back')
    gmsh.model.addPhysicalGroup(2, inlet, name='inlet')
    gmsh.model.addPhysicalGroup(2, bottom, name='bottom')
    gmsh.model.addPhysicalGroup(2, outlet, name='outlet')
    gmsh.model.addPhysicalGroup(2, top, name='top')
    gmsh.model.addPhysicalGroup(2, coolant_to_cladding, name='coolant_to_cladding')
    gmsh.model.addPhysicalGroup(2, cladding_to_coolant, name='cladding_to_coolant')
    gmsh.model.addPhysicalGroup(2, cladding_to_fuel, name='cladding_to_fuel')
    gmsh.model.addPhysicalGroup(2, fuel_to_cladding, name='fuel_to_cladding')

	# Add field for meshing control
    gmsh.model.occ.synchronize()
    distance = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(distance, "FacesList", [coolant_to_cladding[0]])
    gmsh.model.mesh.field.setNumber(distance, "Sampling", 360)
    threshold = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(threshold, "InField", distance)
    gmsh.model.mesh.field.setNumber(threshold, "SizeMin", lc/div)
    gmsh.model.mesh.field.setNumber(threshold, "SizeMax", lc)
    gmsh.model.mesh.field.setNumber(threshold, "DistMin", (1-factor)*outerCladRadius)
    gmsh.model.mesh.field.setNumber(threshold, "DistMax", factor*outerCladRadius)
	
    # Refine mesh near the inlet
    #  inlet_dist = gmsh.model.mesh.field.add("Distance")
    #  gmsh.model.mesh.field.setNumbers(inlet_dist, "FacesList", inlet)
    #  inlet_thre = gmsh.model.mesh.field.add("Threshold")
    #  gmsh.model.mesh.field.setNumber(inlet_thre, "InField", inlet_dist)
    #  gmsh.model.mesh.field.setNumber(inlet_thre, "SizeMin", lc/div)
    #  gmsh.model.mesh.field.setNumber(inlet_thre, "SizeMax", lc)
    #  gmsh.model.mesh.field.setNumber(inlet_thre, "DistMin", (1-factor2)*coreHeight)
    #  gmsh.model.mesh.field.setNumber(inlet_thre, "DistMax", factor2*coreHeight)

    minimum = gmsh.model.mesh.field.add("Min")
    #  gmsh.model.mesh.field.setNumbers(minimum, "FieldsList", [threshold, inlet_thre])
    gmsh.model.mesh.field.setNumbers(minimum, "FieldsList", [threshold])
    #  gmsh.model.mesh.field.setNumbers(minimum, "FieldsList", [inlet_thre])
    gmsh.model.mesh.field.setAsBackgroundMesh(minimum)

    # Change meshing options
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    #  gmsh.option.setNumber("Mesh.Recombine3DAll", 1)    # Apply recombination3D algorithm to all volumes, ignoring per-volume spec (experimental)  
    #  gmsh.option.setNumber("Mesh.Recombine3DLevel", 0)  # 3d recombination level (0: hex, 1: hex+prisms, 2: hex+prism+pyramids) (experimental)
    #  gmsh.option.setNumber("Mesh.Recombine3DConformity", 0)  # 3d recombination conformity type (0: nonconforming, 1: trihedra, 2: pyramids+trihedra, 3:pyramids+hexSplit+trihedra, 4:hexSplit+trihedra)(experimental)
    #  gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 0)  # Mesh recombination algorithm (0: simple, 1: blossom, 2: simple full-quad, 3: blossom full-quad)
    gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1) # Mesh subdivision algorithm (0: none, 1: all quadrangles, 2: all hexahedra, 3: barycentric)
    gmsh.option.setNumber("Mesh.Algorithm", 5) # 2D mesh algorithm (1: MeshAdapt, 2: Automatic, 3: Initial mesh only, 5: Delaunay, 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms, 11: Quasi-structured Quad)
    gmsh.option.setNumber("Mesh.Algorithm3D", 10) # 3D mesh algorithm (1: Delaunay, 3: Initial mesh only, 4: Frontal, 7: MMG3D, 9: R-tree, 10: HXT)
    #  gmsh.option.setNumber("General.NumThreads", 10) # Number of OpenMP threads
    #  gmsh.option.setNumber("General.Verbosity", 100) # Print debug messages

    # Before meshing the model, we need to use the syncronize command
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(3)

    # Remove duplicates
    gmsh.model.mesh.removeDuplicateElements([])
    
    # We can write the mesh to msh to be visualized with gmsh with the following command
    gmsh.write("cartesian_3region.msh")

    #  gmsh.fltk.run()

def main():
    logger_setup()

    DESCRIPTION = """
    AUTHOR/S: {author}

    e-MAIL/S: {email}

    DATE: {date}
    VERSION: {version}
    
    Generate mesh for half a single cartesian pincell with fuel, cladding and coolant regions

    Usage:
        cartesian_3region.py [-h | --help] [-d | --debug]

    Options:
        -h, --help      Show help message and exit.
        -d, --debug     Display debug messages.

    """.format(author=__AUTHOR__, email=__EMAIL__, date=__LAST_MODIFIED__, version=__VERSION__)
    # Parse arguments
    args = docopt(DESCRIPTION)
    logger.info(DESCRIPTION)
    if args['--debug']:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    logger.debug('args = {}'.format(args))
    process(args)

if __name__ == '__main__':
    main()
else:
    logger = logging.getLogger(__name__)



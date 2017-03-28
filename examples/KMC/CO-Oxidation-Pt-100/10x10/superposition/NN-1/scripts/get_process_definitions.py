#!/usr/bin/env python
# -*- coding: utf-8 -*-

from operator import add
import sys

from kynetix import PY2
if not PY2:
    from functools import reduce

from kynetix.utilities.coordinate_utilities import CoordsGroup

def generate_process_dict(origin, ori_coord, coord, rxn_expression, coord_groups):
    """
    Function to get kMC process definition dict.

    Parameters:
    -----------
    origin: The CoordsGroup object of the original origin coordinates group.
    ori_coord: The original coordinates of origin, list of floats.
    coord: The coordinates after moving, list of floats.
    rxn_expression: Reaction expression, str.
    coord_groups: Original coordinates groups, list of CoordsGroup.
    """
    offset = [i - j for i, j in zip(coord, ori_coord)]
    coord_group = reduce(add, coord_groups) + origin.move(offset)
    elements_after = ["V"]*len(coord_group.elements())
    process_dict = dict(reaction=rxn_expression,
                        coordinates_group=[coord_group.coordinates()],
                        elements_before=coord_group.elements(),
                        elements_after=elements_after,
                        basis_sites=[0])

    return process_dict


def tofile(process_dicts, filename):
    with open(filename, "w") as f:
        process_str = "processes = [\n"
        for process_dict in process_dicts:
            process_str += "\n" + str(process_dict) + ",\n"
        process_str += "\n]\n"
        f.write(process_str)


# CoordGroup objects of origins.
# top
c = [[0.0, 0.0, 0.0],   
     [-0.5, 0.0, 0.0],  
     [-0.5, 0.5, 0.0],  
     [0.0, 0.5, 0.0],   
     [0.5, 0.5, 0.0],   
     [0.5, 0.0, 0.0],   
     [0.5, -0.5, 0.0],  
     [0.0, -0.5, 0.0],  
     [-0.5, -0.5, 0.0]]
e = ["C", "V", "V", "V", "V", "V", "V", "V", "V"]
top = CoordsGroup(c, e)

# bri1
c = [[0.0, 0.0, 0.0],  
     [0.5, 0.5, 0.0],  
     [1.0, 0.0, 0.0],  
     [0.5, -0.5, 0.0], 
     [0.5, 0.0, 0.0]]
e = ["V", "V", "V", "V", "C"]
bri1 = CoordsGroup(c, e)

# bri2.
c = [[0.0, 0.0, 0.0],  
     [-0.5, 0.5, 0.0], 
     [0.0, 1.0, 0.0],  
     [0.5, 0.5, 0.0],  
     [0.0, 0.5, 0.0]]
e = ["V", "V", "V", "V", "C"]
bri2 = CoordsGroup(c, e)

origins = [top, bri1, bri2]

if "__main__" == __name__:

    if len(sys.argv) != 2:
        print "Usage: python get_process_definitions.py inputfile"
        sys.exit(1)
    else:
        inputfile = sys.argv[1]

    # Get parameters in input files.
    glob, locs = {}, {}
    exec(open(inputfile, "rb").read(), glob, locs)

    coords_indices = locs["coords_indices"]
    rxn_expression = locs["rxn_expression"]
    ori_coords = locs["ori_coords"]
    coord_groups = locs["coord_groups"]

    # Get all processes.
    process_dicts = []
    for coord_group in coord_groups:
        for coord, idx in coords_indices:
            origin = origins[idx]
            ori_coord = ori_coords[idx]
            process_dict = generate_process_dict(origin,
                                                 ori_coord,
                                                 coord,
                                                 rxn_expression,
                                                 coord_group)
            process_dicts.append(process_dict)

    # Write to file.
    filename = inputfile.split(".")[0] + "_out.py"
    tofile(process_dicts, filename)


#!/usr/bin/env python
# -*- coding: utf-8 -*-

from operator import add

from kynetix.utilities.coordinate_utilities import CoordsGroup

def generate_process_dict(origin, ori_coord, coord, rxn_expression, coord_groups):
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
        for process_dict in process_dicts:
            process_str = "\n" + str(process_dict) + ",\n"
            f.write(process_str)


coords_indices = [
    ([0.0, 0.0, 0.0], 0), ([0.5, 0.0, 0.0], 1), ([1.0, 0.0, 0.0], 0), ([1.5, 0.0, 0.0], 1),
    ([2.0, 0.0, 0.0], 0), ([2.5, 0.0, 0.0], 1), ([3.0, 0.0, 0.0], 0), ([0.0, 0.5, 0.0], 2),
    ([1.0, 0.5, 0.0], 2), ([2.0, 0.5, 0.0], 2), ([3.0, 0.5, 0.0], 2), ([0.0, 1.0, 0.0], 0),
    ([0.5, 1.0, 0.0], 1), ([1.0, 1.0, 0.0], 0), ([1.5, 1.0, 0.0], 1), ([2.0, 1.0, 0.0], 0),
    ([2.5, 1.0, 0.0], 1), ([3.0, 1.0, 0.0], 0), ([0.0, 1.5, 0.0], 2), ([3.0, 1.5, 0.0], 2),
    ([0.0, 2.0, 0.0], 0), ([0.5, 2.0, 0.0], 1), ([2.5, 2.0, 0.0], 1), ([3.0, 2.0, 0.0], 0),
    ([0.0, 2.5, 0.0], 2), ([3.0, 2.5, 0.0], 2), ([0.0, 3.0, 0.0], 0), ([0.5, 3.0, 0.0], 1),
    ([1.0, 3.0, 0.0], 0), ([1.5, 3.0, 0.0], 1), ([2.0, 3.0, 0.0], 0), ([2.5, 3.0, 0.0], 1),
    ([3.0, 3.0, 0.0], 0), ([0.0, 3.5, 0.0], 2), ([1.0, 3.5, 0.0], 2), ([2.0, 3.5, 0.0], 2),
    ([3.0, 3.5, 0.0], 2), ([0.0, 4.0, 0.0], 0), ([0.5, 4.0, 0.0], 1), ([1.0, 4.0, 0.0], 0),
    ([1.5, 4.0, 0.0], 1), ([2.0, 4.0, 0.0], 0), ([2.5, 4.0, 0.0], 1), ([3.0, 4.0, 0.0], 0),
]


if "__main__" == __name__:
    ori_coord = [1.0, 2.0, 0.0]
    coord = [0.0, 4.0, 0.0]
    rxn_expression = "CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b"

    # Coordinates of origin.
    ori_coords = [[1.0, 2.0, 0.0], [1.5, 2.0, 0.0], [1.0, 2.5, 0.0]]

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

    # Get fixed coord_groups.
    c = [[0.0, 0.0, 0.0],  
         [-0.5, 0.5, 0.0], 
         [0.0, 1.0, 0.0],  
         [0.5, 0.5, 0.0],  
         [0.0, 0.5, 0.0]]
    e = ["V", "V", "V", "V", "V"]
    free1 = CoordsGroup(c, e)

    free2 = free1.move([1, 0, 0])
    free3 = free1.move([0, -1, 0])
    free4 = free2.move([0, -1, 0])

    c = [[0.0, 0.0, 0.0],  
         [0.5, 0.5, 0.0],  
         [1.0, 0.0, 0.0],  
         [0.5, -0.5, 0.0], 
         [0.5, 0.0, 0.0]]
    e = ["V", "V", "V", "V", "O_s"]
    O_s = CoordsGroup(c, e)

    coord_groups = [[O_s, free1], [O_s, free2], [O_s, free3], [O_s, free4]]

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
    filename = "kmc_process_1.py"
    tofile(process_dicts, filename)


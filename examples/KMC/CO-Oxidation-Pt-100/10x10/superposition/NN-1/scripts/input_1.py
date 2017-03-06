#!/bin/env/python

from kynetix.utilities.coordinate_utilities import CoordsGroup

coords_indices = [
    # NN-1
    ([-0.5, 1.0, 0.0], 1), ([0.0, 1.0, 0.0], 0), ([0.5, 1.0, 0.0], 1), ([1.0, 1.0, 0.0], 0),
    ([1.5, 1.0, 0.0], 1), ([1.5, 0.0, 0.0], 1), ([1.5, -1.0, 0.0], 1), ([1,0, -1.0, 0.0], 0),
    ([0.5, -1.0, 0.0], 1), ([0.0, -1.0, 0.0], 0), ([-0.5, -1.0, 0.0], 1), ([-0.5, 0.0, 0.0], 0),
]

rxn_expression = "CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b"

# Coordinates of origin.
ori_coords = [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [0.0, 0.5, 0.0]]

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


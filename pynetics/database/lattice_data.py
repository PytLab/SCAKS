'''
    Module to store lattice related data.
'''

grid_neighbor_offsets = {

    'square': ((0.0, -1.0, 0.0), (-1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (1.0, 0.0, 0.0)),  # from left, clockwise

    'hexagonal': ((-1.0, 0.0, 0.0), (-1.0, 1.0, 0.0), (0.0, 1.0, 0.0),            # from left, clockwise
                  (1.0, 0.0, 0.0), (1.0, -1.0, 0.0), (0.0, -1.0, 0.0))
}

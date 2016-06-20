import numpy as np


def collect_coverage(types, possible_types, ncvgs):
    '''
    Function to get current coverages of possible types.

    Parameters:
    -----------
    types: The site types at the lattice points as a list, list of str.

    possible_types: possible species type in grid.

    ncvgs: number of possible species type.

    Returns:
    --------
    cvgs: coverages of possible types, numpy.array int
    '''
    # total number of sites
    nsite = len(types)
    # numbers of different types
    ntypes = [0]*len(possible_types)

    # go through all site to get species numbers
    for sp in types:
        idx = possible_types.index(sp)
        ntypes[idx] += 1

    cvgs = np.array([float(ntype)/nsite for ntype in ntypes])

    return cvgs


def match_elements_list(types,
                        nrow, ncol,
                        stripped_elements_list,
                        stripped_coordinates_list,
                        grid_shape):
    '''
    Function to get total matching success number for,
    a list of stripped elements list and coordinates.

    Parameters:
    -----------
    types: The site types at the lattice points as a list, list of str.

    nrow, ncol: numbers of row and column of lattice grid.

    stripped_elements_list: a list of stripped_elements, a **1D** string list.

    stripped_coordinates_list: a list of stripped coordinates, 2D string list.

    grid_shape: shape of grid, tuple of int.

    Returns:
    --------
    total_nsuccess: total number of successful matching, int
    '''
    # elements_list reshape
    stripped_elements_list = np.array(stripped_elements_list)
    stripped_elements_list.shape = (nrow, ncol)

    total_nsuccess = 0
    for stripped_elements, stripped_coordinates in \
            zip(stripped_elements_list, stripped_coordinates_list):
        n_success = match_elements(types,
                                   stripped_elements,
                                   stripped_coordinates,
                                   grid_shape)
        total_nsuccess += n_success

    return total_nsuccess


def match_elements(types, stripped_elements, stripped_coordinates, grid_shape):
    '''
    Function to go through grid to match elements local configuration.

    Parameters:
    -----------
    types: The site types at the lattice points as a list, list of str.

    stripped_elements: stripped elements list(without wildcards),
                       numpy.array of str.

    stripped_coordinates: stripped relative coordinates list(without wildcards),
                          2d numpy.array of float.

    grid_shape: shape of grid, tuple of int.

    Returns:
    --------
    n_success: number of successful matching, int

    '''
    m, n = grid_shape
    # loop through types
    n_success = 0
    for i in xrange(m):
        for j in xrange(n):
            # check all elements
            match_fail = False
            for element, coordinate in zip(stripped_elements, stripped_coordinates):
                x_offset, y_offset = coordinate[: 2]
                # do pdc check
                x, y = i + int(x_offset), j + int(y_offset)
                # check x
                if x < 0:
                    x = m - 1
                elif x > m-1:
                    x = 0
                # check y
                if y < 0:
                    y = n - 1
                elif y > n-1:
                    y = 0

                idx = x*n + y

                if types[idx] != element:
                    match_fail = True
                    break

            if not match_fail:
                n_success += 1

    return n_success

import numpy as np


def match_elements_list(types,
                        nrow, ncol,
                        stripped_elements_list,
                        stripped_coordinates_list,
                        grid_shape):
    '''
    Function to get total matching success number for,
    a list of stripped elements list and coordinates.
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

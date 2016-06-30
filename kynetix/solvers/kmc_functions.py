import numpy as np


def collect_coverages(types, possible_types):
    '''
    Function to get current coverages of possible types.

    Parameters:
    -----------
    types: The site types at the lattice points as a list, list of str.

    possible_types: possible species type in grid.

    Returns:
    --------
    cvgs: coverages of possible types, numpy.array int
    '''
    # total number of sites
    nsite = len(types)

    # numbers of different types
    ntypes = [0]*len(possible_types)

    # Loop to collect coverages.
    for element in types:
        idx = possible_types.index(element)
        ntypes[idx] += 1

    cvgs = [float(ntype)/nsite for ntype in ntypes]

    return cvgs


import numpy as np


def collect_coverages(types, possible_types, coverage_ratios):
    """
    Function to get current coverages of possible types.

    Parameters:
    -----------
    types: The site types at the lattice points as a list, list of str.

    possible_types: possible species type in grid.

    coverage_ratios: The coverages ratios for all basis sites.

    Returns:
    --------
    cvgs: coverages of possible types, numpy.array int
    """
    # Total number of sites.
    nsite = len(types)/len(coverage_ratios)

    # Numbers of different types.
    ntypes = [0]*len(possible_types)

    # Loop to collect coverages.
    all_ratios = coverage_ratios*nsite
    for ratio, element in zip(all_ratios, types):
        # Ignore the empty type.
        if element not in possible_types:
            continue

        idx = possible_types.index(element)
        ntypes[idx] += ratio

    cvgs = [float(ntype)/nsite for ntype in ntypes]

    return cvgs


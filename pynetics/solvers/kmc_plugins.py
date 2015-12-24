import logging

import numpy as np
try:
    from KMCLib import KMCAnalysisPlugin
except ImportError:
    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    print "!!!                                                   !!!"
    print "!!!          WARNING: KMCLib is not installed         !!!"
    print "!!! Any kMC calculation using KMCLib will be disabled !!!"
    print "!!!                                                   !!!"
    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

from ..functions import get_list_string
from .. import file_header


def collect_coverage(types, possible_types):
    '''
    Function to get current coverages of possible types.

    '''
    # total number of sites
    nsite = len(types)
    # numbers of different types
    ntypes = [0]*len(possible_types)

    # go through all site to get species numbers
    for sp in types:
        idx = possible_types.index(sp)
        ntypes[idx] += 1

    cvgs = [float(ntype)/nsite for ntype in ntypes]

    return cvgs


class CoveragesAnalysis(KMCAnalysisPlugin):
    '''
    Sub-class of KMCLib KMCAnalysisPlugin to do on-the-fly coverage analysis.
    '''
    def __init__(self, kinetic_model):
        # get possible elements type
        adsorbate_names = [ads.split('_')[0] for ads in
                           kinetic_model.adsorbate_names]
        possible_types = adsorbate_names + ['Vac']
        self.possible_types = possible_types

        # initialize recorder variables
        self.steps = []
        self.times = []

        # initialize surface species cvg
        self.species_cvgs = []

        # set logger
        self.logger = logging.getLogger('model.solvers.CoveragesAnalysis')

    def setup(self, step, time, configuration):
        '''
        Function called right before the start of the KMC loop to allow for
        custom setup of the analysis object.

        :param step: The simulation step number.
        :type step: int

        :param time: The simulation time.
        :type time: float

        :param configuration: The up to date configuration of the simulation.
        :type configuration: KMCConfiguration
        '''
        self.steps.append(step)
        self.times.append(time)

        # get species coverages
        types = configuration.types()
        cvgs = collect_coverage(types, self.possible_types)
        self.species_cvgs.append(cvgs)

    def registerStep(self, step, time, configuration):
        '''
        Called from the KMC loop after each step.

        :param step: The simulation step number.
        :type step: int

        :param time: The simulation time.
        :type time: float

        :param configuration: The up to date configuration of the simulation.
        :type configuration: KMCConfiguration
        '''
        self.steps.append(step)
        self.times.append(time)

        # get species coverages
        types = configuration.types()
        cvgs = collect_coverage(types, self.possible_types)
        self.species_cvgs.append(cvgs)

    def finalize(self):
        '''
        Called after the KMC loop to allow for custom finalization and
        post-processing of collected data.
        '''
        # get variables strings
        cvgs = zip(*self.species_cvgs)
        cvgs_str = get_list_string('coverages', cvgs)
        times_str = get_list_string('times', self.times)
        steps_str = get_list_string('steps', self.steps)
        # get possible types
        possible_types_str = get_list_string('possible_types',
                                             self.possible_types)

        # write to file
        content = (file_header + times_str + steps_str +
                   cvgs_str + possible_types_str)
        with open('auto_coverages.py', 'w') as f:
            f.write(content)

        self.logger.info('coverages info are wirtten to auto_coverages.py.')

        return


class TOFAnalysis(KMCAnalysisPlugin):
    '''
    Sub-class of KMCLib KMCAnalysisPlugin to do on-the-fly TOF analysis.
    '''
    def __init__(self, kinetic_model, rxn_indices):
        self.kinetic_model = kinetic_model
        self.rxn_indices = rxn_indices

    def setup(self, step, time, configuration):
        '''
        Function called right before the start of the KMC loop to allow for
        custom setup of the analysis object.

        :param step: The simulation step number.
        :type step: int

        :param time: The simulation time.
        :type time: float

        :param configuration: The up to date configuration of the simulation.
        :type configuration: KMCConfiguration
        '''
        # get elementary reaction rates and elements changes
        rates_list = []
        elements_changes_list = []
        for idx in self.rxn_indices:
            rxn_list = self.kinetic_model.elementary_rxns_list[idx]
            # elementary forward and reversed rates
            rates = self.kinetic_model.parser.get_elementarty_rate(rxn_list)
            rates_list.append(rates)
            # elementary reaction elements changes
            elements_changes = self.kinetic_model.parser.get_elementary_elements_changes(rxn_list)
            elements_changes_list.append(elements_changes)

        self.rates_list = rates_list
        self.elements_changes_list = elements_changes_list

    def registerStep(self, step, time, configuration):
        '''
        Called from the KMC loop after each step.

        :param step: The simulation step number.
        :type step: int

        :param time: The simulation time.
        :type time: float

        :param configuration: The up to date configuration of the simulation.
        :type configuration: KMCConfiguration
        '''

    def finalize(self):
        '''
        Called after the KMC loop to allow for custom finalization and
        post-processing of collected data.
        '''


def strip_local_configuration(elements, coordinates):
    '''
    Function to strip wild cards off from elements and corresponding coordinates.

    Parameters:
    -----------
    elements: local elements around a center elements, numpy.array.

    coordinates: corresponding coordinates of elements, 2d numpy.array.

    Example:
    --------
    >>> e = np.array(['CO', '*', '*', 'O', '*'])
    >>> c = array([[ 0., -1.,  0.],
                   [-1.,  0.,  0.],
                   [ 0.,  1.,  0.],
                   [ 1.,  0.,  0.]])

    >>> (array(['CO', 'O'], dtype='|S2'),
         array([[ 0.,  0.,  0.], [ 0.,  1.,  0.]]))

    '''
    indices = [idx for idx, elem in enumerate(elements) if elem != '*']
    # mask elements
    stripped_elements = np.array([elements[idx] for idx in indices])
    # mask coordinates
    stripped_coordinates = np.array([coordinates[idx] for idx in indices])

    return stripped_elements, stripped_coordinates


def match_elements(types, stripped_elements, stripped_coordinates, grid_shape):
    '''
    Function to go through grid to match elements local configuration.

    Parameters:
    -----------
    types: The site types at the lattice points as a list, list of str.

    stripped_elements: stripped elements list(without wildcards),
                       numpy.array of str.

    stripped_coordinates: stripped coordinate list(without wildcards),
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

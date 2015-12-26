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
    def __init__(self, kinetic_model):
        self.kinetic_model = kinetic_model
        self.rxn_indices = kinetic_model.TOF_indices
        # get total site number
        grid_shape = self.kinetic_model.grid_shape
        self.Ntot = reduce(lambda x, y: x*y, grid_shape)

        # set logger
        self.logger = logging.getLogger('model.solvers.TOFAnalysis')

    def setup(self, step, time, configuration):
        '''
        Function called right before the start of the KMC loop to allow for
        custom setup of the analysis object.

        Get forward and reversed rates for each elementary reaction and get local
        configuration changed list for each elementary reaction.

        :param step: The simulation step number.
        :type step: int

        :param time: The simulation time.
        :type time: float

        :param configuration: The up to date configuration of the simulation.
        :type configuration: KMCConfiguration
        '''
        self.logger.info(' ')
        self.logger.info('TOFAnalysis setting up...')
        parser = self.kinetic_model.parser

        # get elementary reaction rates and elements changes
        rates_list = []
        elements_changes_list = []
        for idx in self.rxn_indices:
            # output rxn expression info
            self.logger.info('-'*20)
            rxn_expression = self.kinetic_model.rxn_expressions[idx]
            self.logger.info('[ %s ]', rxn_expression)

            rxn_list = self.kinetic_model.elementary_rxns_list[idx]

            # elementary forward and reversed rates
            rates = self.kinetic_model.solver.get_elementary_rate(rxn_list)
            rates_list.append(rates)

            # elementary reaction elements changes
            elements_changes = parser.get_elementary_elements_changes(rxn_list)
            elements_changes_list.append(elements_changes)
            # elements change info output
            self.logger.info('elements changes =')
            for elements_change in elements_changes:
                self.logger.info('    %s', str(elements_change))

        # set attrs
        self.rates_list = rates_list
        self.elements_changes_list = elements_changes_list

        # get local relative coordinates
        self.coordinates = self.kinetic_model.solver.get_coordinates()
        # set delta t
        self.last_time = time
        # set temp total rate
        self.total_rates = [[0.0, 0.0] for i in xrange(len(rates_list))]
        # set counter
        self.analysis_counter = 0
        # step and time collection ccontainers
        self.steps, self.times = [], []
        # TOF collection
        self.TOFs = [[] for i in xrange(len(rates_list))]

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
        # collect step and time
        self.steps.append(step)
        self.times.append(time)

        # info output
        self.analysis_counter += 1
        self.logger.info('-------- Entering TOF analysis ( %d ) --------',
                         self.analysis_counter)

        # get current types
        types = configuration.types()
        grid_shape = self.kinetic_model.grid_shape

        # get elapsed time
        current_time = time
        dt = current_time - self.last_time
        self.last_time = current_time
        self.end_time = time

        # count available reaction number in the current configuration
        for idx, (rates, elements_changes) in \
                enumerate(zip(self.rates_list, self.elements_changes_list)):
            # output rxn expression info
            rxn_expression = self.kinetic_model.rxn_expressions[self.rxn_indices[idx]]

            # get rate constants
            Rf, Rr = rates

            # get elements changes for both direction
            reactant_configs, product_configs = zip(*elements_changes)
            # strip reactant elements and coordinates
            stripped_reactants = [strip_local_configuration(reactant_config,
                                                            self.coordinates)
                                  for reactant_config in reactant_configs]
            stripped_reactant_elems, stripped_reactant_coords = zip(*stripped_reactants)
            # strip product elements and coordinates
            stripped_products = [strip_local_configuration(product_config,
                                                           self.coordinates)
                                 for product_config in product_configs]
            stripped_product_elems, stripped_product_coords = zip(*stripped_products)

            # number of successfule matching for forward
            n_forward = match_elements_list(types,
                                            stripped_reactant_elems,
                                            stripped_reactant_coords,
                                            grid_shape)
            total_forward_rate = Rf*n_forward*dt

            # number of successfule matching for reversed
            n_reversed = match_elements_list(types,
                                             stripped_product_elems,
                                             stripped_product_coords,
                                             grid_shape)
            total_reversed_rate = Rr*n_reversed*dt

            # add total rate to total rates list
            self.total_rates[idx][0] += total_forward_rate
            self.total_rates[idx][-1] += total_reversed_rate

            # info output
            self.logger.info('[ %s ] n_forward = %d, n_reversed = %d, dt = %e',
                             rxn_expression, n_forward, n_reversed, dt)

        # collect TOFs for this step
        current_tofs = self.append_TOFs()
        self.logger.info('current TOFs = ')
        for tof_list in current_tofs:
            self.logger.info('    %s', str(tof_list))
        self.logger.info(' ')

    def finalize(self):
        '''
        Called after the KMC loop to allow for custom finalization and
        post-processing of collected data.

        Write steps, time, tof trajectories to file.
        '''
        times_str = get_list_string('times', self.times)
        steps_str = get_list_string('steps', self.steps)
        # get TOFs string
        TOFs_str = ''
        for idx, tofs in enumerate(self.TOFs):
            elementary_str = get_list_string('TOFs_'+str(idx), tofs)
            TOFs_str += elementary_str

        content = file_header + times_str + steps_str + TOFs_str

        with open('auto_TOFs.py', 'w') as f:
            f.write(content)

        self.logger.info('TOFs info are written to auto_TOFs.py.')

    def append_TOFs(self):
        '''
        Function to calculate TOF values from statistic number, and
        append them to self.TOFs, return current_tofs list.
        '''
        current_tofs = []
        for idx, rate_list in enumerate(self.total_rates):
            tof_list = [rate/(self.Ntot*self.end_time) for rate in rate_list]
            self.TOFs[idx].append(tof_list)
            current_tofs.append(tof_list)

        return current_tofs


def match_elements_list(types,
                        stripped_elements_list,
                        stripped_coordinates_list,
                        grid_shape):
    '''
    Function to get total matching success number for,
    a list of stripped elements list and coordinates.
    '''
    total_nsuccess = 0
    for stripped_elements, stripped_coordinates in \
            zip(stripped_elements_list, stripped_coordinates_list):
        n_success = match_elements(types,
                                   stripped_elements,
                                   stripped_coordinates,
                                   grid_shape)
        total_nsuccess += n_success

    return total_nsuccess


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
    >>> c = array([[ 0.,  0.,  0.],
                   [ 0., -1.,  0.],
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

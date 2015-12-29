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
try:
    from .plugin_backends.kmc_functions import *
except ImportError:
    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    print "!!!   WARNING: plugin backends extension not found.   !!!"
    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    from kmc_functions import *


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
        self.ncvgs = len(possible_types)

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
        cvgs = collect_coverage(types, self.possible_types, self.ncvgs)
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
        cvgs = collect_coverage(types, self.possible_types, self.ncvgs)
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
        self.logger.info('collect statistics about possible reaction:')

        # get current types
        types = configuration.types()
        grid_shape = self.kinetic_model.grid_shape

        # get elapsed time
        current_time = time
        dt = current_time - self.last_time
        self.last_time = current_time
        self.end_time = time

        def collect_statistics(species_configs):
            "collect all possible reaction number for a list of configurations."
            # strip reactant elements and coordinates
            stripped_species = [strip_local_configuration(species_config,
                                                          self.coordinates)
                                for species_config in species_configs]
            stripped_species_elems, stripped_species_coords = zip(*stripped_species)

            # datatype conversion before statistics
            stripped_species_elems = np.array(stripped_species_elems)
            nrow, ncol = stripped_species_elems.shape
            # elements_list to 1D list
            stripped_species_elems, = stripped_species_elems.reshape(1, -1).tolist()
            # coordinates list to 3D numpy.array
            stripped_species_coords = np.array(stripped_species_coords)

            # number of successfule matching for forward
            n = match_elements_list(types,
                                    nrow, ncol,
                                    stripped_species_elems,
                                    stripped_species_coords,
                                    grid_shape)

            return n

        # count available reaction number in the current configuration
        for idx, (rates, elements_changes) in \
                enumerate(zip(self.rates_list, self.elements_changes_list)):
            # output rxn expression info
            rxn_expression = self.kinetic_model.rxn_expressions[self.rxn_indices[idx]]

            # get rate constants
            Rf, Rr = rates

            # get elements changes for both direction
            reactant_configs, product_configs = zip(*elements_changes)

            # get total sum of ka of forward reaction
            n_forward = collect_statistics(reactant_configs)
            total_forward_rate = Rf*n_forward*dt

            # get total sum of ka of reversed reaction
            n_reversed = collect_statistics(reactant_configs)
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
            self.logger.info('   %17.10e(+), %e(-)', *tof_list)
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


def strip_local_configuration(elements, coordinates):
    '''
    Function to strip wild cards off from elements and corresponding coordinates.

    Parameters:
    -----------
    elements: local elements around a center elements, list of strings.

    coordinates: corresponding coordinates of elements, 2d list of float.

    Example:
    --------
    >>> e = ['CO', '*', '*', 'O', '*']
    >>> c = array([[ 0.,  0.,  0.],
                   [ 0., -1.,  0.],
                   [-1.,  0.,  0.],
                   [ 0.,  1.,  0.],
                   [ 1.,  0.,  0.]])

    >>> (['CO', 'O'], [[ 0.,  0.,  0.], [ 0.,  1.,  0.]])

    '''
    indices = [idx for idx, elem in enumerate(elements) if elem != '*']
    # mask elements
    stripped_elements = [elements[idx] for idx in indices]
    # mask coordinates
    stripped_coordinates = [coordinates[idx] for idx in indices]

    return stripped_elements, stripped_coordinates

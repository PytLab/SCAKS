#####################################################################
#           KMC plugins for coverages and TOF analysis
#
# * CoveragesAnalysis and TOFAnalysis are sub-classes of KynetixPlugin,
#   and support continuous job calculations.
#
# * TOFCoveragesAnalysis is sub-class of TOFAnalysis,
#   dose not support continuous job calculations.
#####################################################################
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


def reset_step_and_time(func):
    '''
    Decorator function to decorate analysis plugin methods.
    '''
    def wrapped_func(self, step, time, configuration):
        step = step + self.step_base
        time = time + self.time_base
        func(self, step, time, configuration)

    return wrapped_func


class KynetixPlugin(KMCAnalysisPlugin):
    '''
    Base class from which analysis classes in kynetix inherit.
    '''
    def __init__(self, kinetic_model):
        self.kinetic_model = kinetic_model

        # ---------------------------------------------------------
        # base point is used to get correct time and step value
        # when the job is a continuous job
        #
        # start_time and start_step are got in load() in /model.py
        # ---------------------------------------------------------
        # set step base point
        if hasattr(self.kinetic_model, 'start_step'):
            self.step_base = self.kinetic_model.start_step
        else:
            self.step_base = 0

        # set time base point
        if hasattr(self.kinetic_model, 'start_time'):
            self.time_base = self.kinetic_model.start_time
        else:
            self.time_base = 0.0

        # total step number for kMC loop
        self.total_step = self.kinetic_model.nstep

        # condition for continuous job or not
        self.continuous_job = (hasattr(self.kinetic_model, 'kmc_continue') and
                               self.kinetic_model.kmc_continue)

        # dump interval
        if hasattr(self.kinetic_model, 'analysis_dump_interval'):
            self.dump_interval = self.kinetic_model.analysis_dump_interval
        else:
            self.dump_interval = 1

        # set counter
        self.analysis_counter = 0

    def store_last_types(self):
        '''
        Archive self.last_types by writting it to file.
        '''
        last_types_str = get_list_string('types', self.last_types)
        content = file_header + last_types_str

        with open('auto_last_types.py', 'w') as f:
            f.write(content)
        self.logger.info('types info written into auto_last_types.py')


class CoveragesAnalysis(KynetixPlugin):
    '''
    Sub-class of KMCLib KynetixPlugin to do on-the-fly coverage analysis.
    '''
    def __init__(self, kinetic_model):
        KynetixPlugin.__init__(self, kinetic_model)
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

    @reset_step_and_time
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

    @reset_step_and_time
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
        self.last_types = configuration.types()

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

        # write to coverages file
        content = (file_header + times_str + steps_str +
                   cvgs_str + possible_types_str)
        if self.continuous_job:
            filename = 'auto_coverages_ctn.py'
        else:
            filename = 'auto_coverages.py'
        with open(filename, 'w') as f:
            f.write(content)
        self.logger.info('coverages info are written to ' + filename + '.')

        # write to last types file
        self.store_last_types()

        return


class TOFAnalysis(KynetixPlugin):
    '''
    Sub-class of KMCLib KynetixPlugin to do on-the-fly TOF analysis.
    '''
    def __init__(self, kinetic_model):
        KynetixPlugin.__init__(self, kinetic_model)

        self.rxn_indices = kinetic_model.TOF_indices
        # get total site number
        grid_shape = self.kinetic_model.grid_shape
        self.Ntot = reduce(lambda x, y: x*y, grid_shape)

        # set logger
        self.logger = logging.getLogger('model.solvers.TOFAnalysis')

    @reset_step_and_time
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
        # set total rate
        if hasattr(self.kinetic_model, 'total_rates'):
            self.total_rates = self.kinetic_model.total_rates
        else:
            self.total_rates = [[0.0, 0.0] for i in xrange(len(rates_list))]
        # step and time collection ccontainers
        self.steps, self.times = [], []
        # TOF collection
        self.TOFs = [[] for i in xrange(len(rates_list))]

        # variables for TOF analysis delay
        # if it's a continuous job, start at the first step
        if self.continuous_job:
            self.start = True
        else:
            self.start = False

        # if it's a continuous job, taf start time is the same as
        # that in last kmc job
        if hasattr(self.kinetic_model, 'tof_start_time'):
            self.start_time = self.kinetic_model.tof_start_time
        else:
            self.start_time = 0.0

        # use TOF_start_step value if set in setup file
        # or if it's a continuous job, start at the first step(set to 0)
        if self.continuous_job or not hasattr(self.kinetic_model, 'TOF_start_step'):
            self.start_step = 0
        else:
            self.start_step = self.kinetic_model.TOF_start_step

        # archive current types
        self.last_types = configuration.types()

    @reset_step_and_time
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
        # process percentage
        percentage = int(float(step - self.step_base)/float(self.total_step)*100)

        # TOF analysis start when coverage is converged
        if step < self.start_step:
            self.logger.info("[ %d%% ] TOF analysis skipped ( step: %d, time: %e s )",
                             percentage, step, time)
            return

        # TOF start stamp
        if not self.start:
            self.start = True
            self.start_time = time
            self.logger.info('[ %d%% ] TOF analysis starts (%d) at time = %e s\n',
                             percentage, step, self.start_time)
            return

        # collect step and time
        self.steps.append(step)
        self.times.append(time)

        # info output
        self.analysis_counter += 1
        if self.analysis_counter % self.dump_interval == 0:
            self.logger.info('[ %d%% ] Entering TOF analysis %d ( step: %d, time: %e s ) ',
                             percentage, self.analysis_counter, step, time)
            self.logger.info('collect statistics about possible reactions:')

        # get current types
        types = configuration.types()
        grid_shape = self.kinetic_model.grid_shape

        # get elapsed time
        dt = time - self.last_time
        self.last_time = time
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
            n_reversed = collect_statistics(product_configs)
            total_reversed_rate = Rr*n_reversed*dt

            # add total rate to total rates list
            self.total_rates[idx][0] += total_forward_rate
            self.total_rates[idx][-1] += total_reversed_rate

            # info output

            if self.analysis_counter % self.dump_interval == 0:
                self.logger.info('[ %s ] n_forward = %d, n_reversed = %d, dt = %e',
                                 rxn_expression, n_forward, n_reversed, dt)

        # collect TOFs for this step
        current_tofs = self.append_TOFs()

        if self.analysis_counter % self.dump_interval == 0:
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
        start_time_str = '\nstart_time = %e\n\n' % self.start_time  # time when TOF analysis begins
        times_str = get_list_string('times', self.times)
        steps_str = get_list_string('steps', self.steps)
        total_rates_str = get_list_string('total_rates', self.total_rates)
        # get TOFs string
        TOFs_str = ''
        for idx, tofs in enumerate(self.TOFs):
            elementary_str = get_list_string('TOFs_'+str(idx), tofs)
            TOFs_str += elementary_str

        content = (file_header + start_time_str + times_str +
                   steps_str + TOFs_str + total_rates_str)

        if self.continuous_job:
            filename = 'auto_TOFs_ctn.py'
        else:
            filename = 'auto_TOFs.py'

        with open(filename, 'w') as f:
            f.write(content)

        self.logger.info('TOFs info are written to ' + filename + '.')

        # write to last types file
        self.store_last_types()

    def append_TOFs(self):
        '''
        Function to calculate TOF values from statistic number, and
        append them to self.TOFs, return current_tofs list.
        '''
        current_tofs = []
        time_span = self.end_time - self.start_time
        for idx, rate_list in enumerate(self.total_rates):
            tof_list = [rate/(self.Ntot*time_span) for rate in rate_list]
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


class TOFCoveragesAnalysis(TOFAnalysis):
    '''
    Sub-class of TOFAnalysis to do on-the-fly TOF and coverages analysis.
    '''
    def __init__(self, kinetic_model):
        self.kinetic_model = kinetic_model
        # get total site number
        grid_shape = self.kinetic_model.grid_shape
        self.Ntot = reduce(lambda x, y: x*y, grid_shape)

        # step and time collection ccontainers
        self.cvg_steps, self.cvg_times = [], []
        self.tof_steps, self.tof_times = [], []

        # set attrs for TOF analysis
        self.rxn_indices = kinetic_model.TOF_indices

        # set attrs for coverages analysis
        adsorbate_names = [ads.split('_')[0] for ads in
                           kinetic_model.adsorbate_names]
        possible_types = adsorbate_names + ['Vac']
        self.possible_types = possible_types
        self.ncvgs = len(possible_types)
        # initialize surface species cvg
        self.species_cvgs = []

        self.total_step = self.kinetic_model.nstep

        # dump interval
        if hasattr(self.kinetic_model, 'analysis_dump_interval'):
            self.dump_interval = self.kinetic_model.analysis_dump_interval
        else:
            self.dump_interval = 1

        # set logger
        self.logger = logging.getLogger('model.solvers.TOFCoveragesAnalysis')

    def setup(self, step, time, configuration):
        '''
        Function called right before the start of the KMC loop to allow for
        custom setup of the analysis object.

        Get forward and reversed rates for each elementary reaction and get local
        configuration changed list for each elementary reaction.

        Get initial species coverages on surface.

        :param step: The simulation step number.
        :type step: int

        :param time: The simulation time.
        :type time: float

        :param configuration: The up to date configuration of the simulation.
        :type configuration: KMCConfiguration
        '''
        #-----------------------------------------
        #                TOF setup
        #-----------------------------------------
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
        self.tof_analysis_counter = 0
        # TOF collection
        self.TOFs = [[] for i in xrange(len(rates_list))]

        # variables for TOF analysis delay
        self.start = False
        self.start_time = 0.0
        if not hasattr(self.kinetic_model, 'TOF_start_step'):
            self.start_step = 0
        else:
            self.start_step = self.kinetic_model.TOF_start_step

        #----------------------------------------------
        #            coverages analysis setup
        #----------------------------------------------
        self.cvg_steps.append(step)
        self.cvg_times.append(time)
        # get species coverages
        types = configuration.types()
        cvgs = collect_coverage(types, self.possible_types, self.ncvgs)
        self.species_cvgs.append(cvgs)

    def registerStep(self, step, time, configuration):
        '''
        Called from the KMC loop after each step.

        Collect statistics about coverages and TOFs of the configuration.

        :param step: The simulation step number.
        :type step: int

        :param time: The simulation time.
        :type time: float

        :param configuration: The up to date configuration of the simulation.
        :type configuration: KMCConfiguration
        '''
        percentage = int(float(step)/float(self.total_step)*100)
        self.last_types = configuration.types()

        #-----------------------------------------------
        #    do coverage analysis before tof analysis
        #-----------------------------------------------
        self.cvg_steps.append(step)
        self.cvg_times.append(time)
        # get species coverages
        types = configuration.types()
        current_cvgs = collect_coverage(types, self.possible_types, self.ncvgs)
        old_cvgs = self.species_cvgs[-1]
        self.species_cvgs.append(current_cvgs)

        # check covergaes convergency
        old = np.array(old_cvgs)
        current = np.array(current_cvgs)
        deltas = np.abs(current - old)
        converged = (deltas < self.kinetic_model.cvgs_convergency).all()

        # TOF analysis start when coverage is converged
        if not converged and not self.start:
            self.logger.info("[ %d%% ] Coverages not converged, TOF analysis skipped " +
                             "( step: %d, time: %e s )", percentage, step, time)
            return

        # -----------------------------------------------
        # TOF analysis starts when coverages converged
        # -----------------------------------------------
        # TOF start stamp
        if not self.start:
            self.start = True
            self.start_time = time
            self.logger.info('[ %d%% ] TOF analysis starts (%d) at time = %e s\n',
                             percentage, step, self.start_time)
            return

        # collect step and time
        self.tof_steps.append(step)
        self.tof_times.append(time)

        # info output
        self.tof_analysis_counter += 1

        if self.tof_analysis_counter % self.dump_interval == 0:
            self.logger.info('[ %d%% ] Entering TOF analysis %d ( step: %d, time: %e s )',
                             percentage, self.tof_analysis_counter, step, time)
            self.logger.info('collect statistics about possible reaction:')

        # get current types
        types = configuration.types()
        grid_shape = self.kinetic_model.grid_shape

        # get elapsed time
        dt = time - self.last_time
        self.last_time = time
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
            n_reversed = collect_statistics(product_configs)
            total_reversed_rate = Rr*n_reversed*dt

            # add total rate to total rates list
            self.total_rates[idx][0] += total_forward_rate
            self.total_rates[idx][-1] += total_reversed_rate

            # info output

            if self.tof_analysis_counter % self.dump_interval == 0:
                self.logger.info('[ %s ] n_forward = %d, n_reversed = %d, dt = %e',
                                 rxn_expression, n_forward, n_reversed, dt)

        # collect TOFs for this step
        current_tofs = self.append_TOFs()

        if self.tof_analysis_counter % self.dump_interval == 0:
            self.logger.info('current TOFs = ')
            for tof_list in current_tofs:
                self.logger.info('   %17.10e(+), %e(-)', *tof_list)
            self.logger.info(' ')

    def finalize(self):
        '''
        Called after the KMC loop to allow for custom finalization and
        post-processing of collected data.

        Write steps, time, tof and coverages trajectories to files.
        '''
        #----------------------------------------------
        #          generate TOF analysis file
        #----------------------------------------------
        start_time_str = '\nstart_time = %e\n\n' % self.start_time  # time when TOF analysis begins
        times_str = get_list_string('times', self.tof_times)
        steps_str = get_list_string('steps', self.tof_steps)
        total_rates_str = get_list_string('total_rates', self.total_rates)
        # get TOFs string
        TOFs_str = ''
        for idx, tofs in enumerate(self.TOFs):
            elementary_str = get_list_string('TOFs_'+str(idx), tofs)
            TOFs_str += elementary_str

        content = (file_header + start_time_str + times_str +
                   steps_str + TOFs_str + total_rates_str)

        with open('auto_TOFs.py', 'w') as f:
            f.write(content)

        self.logger.info('TOFs info are written to auto_TOFs.py.')

        #----------------------------------------------
        #          generate coverages analysis file
        #----------------------------------------------
        cvgs = zip(*self.species_cvgs)
        cvgs_str = get_list_string('coverages', cvgs)
        times_str = get_list_string('times', self.cvg_times)
        steps_str = get_list_string('steps', self.cvg_steps)
        # get possible types
        possible_types_str = get_list_string('possible_types',
                                             self.possible_types)

        # write to file
        content = (file_header + times_str + steps_str +
                   cvgs_str + possible_types_str)
        with open('auto_coverages.py', 'w') as f:
            f.write(content)

        self.logger.info('coverages info are written to auto_coverages.py.')

        self.store_last_types()


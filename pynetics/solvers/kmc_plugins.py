import logging

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
    def __init__(self, possible_types):
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

        # write to file
        content = file_header + times_str + steps_str + cvgs_str
        with open('auto_coverages.py', 'w') as f:
            f.write(content)

        self.logger.info('write coverages info to auto_coverages.py.')

        return

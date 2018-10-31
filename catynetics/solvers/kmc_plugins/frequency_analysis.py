import logging
from operator import mul

try:
    from KMCLib import KMCAnalysisPlugin
except ImportError:
    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    print("!!!                                                    !!!")
    print("!!!         WARNING: KMCLibX is not installed          !!!")
    print("!!! Any kMC calculation using KMCLibX will be disabled !!!")
    print("!!!                                                    !!!")
    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

from ... import file_header
from ...compatutil import reduce
from ...mpicommons import mpi
from ...utilities.format_utilities import get_list_string, get_dict_string


class FrequencyAnalysis(KMCAnalysisPlugin):
    """
    KMC plugin to do On-The-Fly process occurence frequency analysis.
    """
    # {{{
    def __init__(self,
                 kmc_model,
                 filename="auto_frequency.py",
                 tof_start=0):
        """
        Constructor of FrequencyAnalysis object.

        Parameters:
        -----------
        kmc_model: KMC model object of catynetics.KineticModel.

        filename: The name of data file, str.

        """
        # LatticeModel object.
        self.__kmc_model = kmc_model

        # Recorder variables.
        self.__times = []
        self.__steps = []

        # Process pick statistics list.
        nprocess = len(kmc_model.processes)
        self.__process_occurencies = [0]*nprocess
        self.__steady_process_occurencies = [0]*nprocess

        # Set logger.
        if mpi.is_master:
            self.__logger = logging.getLogger("model.solvers.KMCSolver.FrequencyAnalysis")

        # TOF start step.
        self.__tof_start_time = 0.0
        self.__tof_end_time = 0.0

        # Name of data file.
        self.__filename = filename

    def setup(self, step, time, configuration, interactions):
        # Create statistic data file.
        variables_str = ("times = []\npicked_indices = []\n" +
                         "process_occurencies = []\n")
        with open(self.__filename, "w") as f:
            content = file_header + variables_str
            f.write(content)

    def registerStep(self, step, time, configuration, interactions):
        # Append picked index.
        picked_index = interactions.pickedIndex()

        # Add to collection list.
        self.__process_occurencies[picked_index] += 1

        # Collect steady frequency info.
        if step >= self.__kmc_model.tof_start:
            self.__steady_process_occurencies[picked_index] += 1

            if not self.__tof_start_time:
                self.__tof_start_time = time
                if mpi.is_master:
                    self.__logger.info("TOF analysis start at time = {:e}".format(time))

        self.__tof_end_time = time

    def finalize(self):
        """
        Write all data to files.
        """
        # Write process occurencies to file.
        occurencies_str = get_list_string("process_occurencies",
                                          self.__process_occurencies)

        # Calculate reaction occurencies.
        process_mapping = self.__kmc_model.process_mapping

        # Construct reaction occurencies dict.
        reaction_occurencies = {}
        steady_reaction_occurencies = {}

        for reaction in set(process_mapping):
            reaction_occurencies.setdefault(reaction, 0)
            steady_reaction_occurencies.setdefault(reaction, 0)

        # Fill the dict.
        for occurency, reaction in zip(self.__process_occurencies, process_mapping):
            reaction_occurencies[reaction] += occurency

        for occurency, reaction in zip(self.__steady_process_occurencies, process_mapping):
            steady_reaction_occurencies[reaction] += occurency

        reaction_occurencies_str = get_dict_string("reaction_occurencies",
                                                   reaction_occurencies)
        steady_reaction_occurencies_str = get_dict_string("steady_reaction_occurencies",
                                                          steady_reaction_occurencies)

        # Get number of active sites.
        repetitions = self.__kmc_model.repetitions
        basis_sites = self.__kmc_model.basis_sites
        nsites = reduce(mul, repetitions)*len(basis_sites)

        delta_t = self.__tof_end_time - self.__tof_start_time

        # Calculate rates.
        reaction_rates = {}
        for reaction, occurency in steady_reaction_occurencies.items():
            rate = occurency/(delta_t*nsites)
            reaction_rates.setdefault(reaction, rate)
        reaction_rates_str = get_dict_string("reaction_rates", reaction_rates)

        if mpi.is_master:
            with open(self.__filename, "a") as f:
                all_content = (occurencies_str + reaction_occurencies_str +
                               steady_reaction_occurencies_str + reaction_rates_str)
                f.write(all_content)

            # Info output.
            msg = "frequency informations are written to {}".format(self.__filename)
            self.__logger.info(msg)
    # }}}


from copy import deepcopy
import logging

import numpy as np
try:
    from KMCLib import KMCAnalysisPlugin
except ImportError:
    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    print "!!!                                                    !!!"
    print "!!!         WARNING: KMCLibX is not installed          !!!"
    print "!!! Any kMC calculation using KMCLibX will be disabled !!!"
    print "!!!                                                    !!!"
    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

try:
    from kynetix.solvers.plugin_backends.kmc_functions import *
except ImportError:
    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    print "!!!   WARNING: plugin backends extension not found.   !!!"
    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    from kynetix.solvers.kmc_functions import *

from kynetix import file_header
from kynetix.functions import get_list_string


class CoveragesAnalysis(KMCAnalysisPlugin):
    """
    KMC plugin to do On-The-Fly coverage analysis.
    """
    def __init__(self, kmc_model, filename="auto_coverages.py"):
        """
        Constructor of CoverageAnalysis object.

        Parameters:
        -----------
        kmc_model: KMC model object of Kynetix.KineticModel.
        """
        # LatticeModel object.
        self.__kmc_model = kmc_model

        # The ratio of a basis site occupied.
        self.__coverage_ratios = (1.0, 1./2, 1./2, 1./4)

        # Recorder variables.
        self.__times = []
        self.__steps = []
        self.__coverages = []

        self.__possible_types = kmc_model.possible_element_types()
        
        # Set logger.
        self.__logger = logging.getLogger("model.solvers.KMCSolver.CoveragesAnalysis")

        # Set data file name.
        self.__filename = filename

    def setup(self, step, time, configuration, interactions):
        # Append time and step.
        self.__times.append(time)
        self.__steps.append(step)

        # Remove empty type from possible_types.
        possible_types_copy = deepcopy(self.__possible_types)
        empty_type = self.__kmc_model.empty_type()
        empty_type_idx = possible_types_copy.index(empty_type)
        possible_types_copy.pop(empty_type_idx)

        # Collect species coverages.
        types = configuration.types()
        coverages = collect_coverages(types,
                                      possible_types_copy,
                                      self.__coverage_ratios)

        # Insert coverage of emtpy site.
        coverages = list(coverages)
        empty_coverage = 1.0 - sum(coverages)
        coverages.insert(empty_type_idx, empty_coverage)

        self.__coverages.append(coverages)

    def registerStep(self, step, time, configuration, interactions):
        # Do the same thing in setup().
        self.setup(step, time, configuration, interactions)

    def finalize(self):
        """
        Write all data to files.
        """
        # Get data strings.
        coverages = zip(*self.__coverages)
        coverages_str = get_list_string("coverages", coverages)
        times_str = get_list_string("times", self.__times)
        steps_str = get_list_string("steps", self.__steps)
        possible_types_str = get_list_string("possible_types", self.__possible_types)

        # Write to file.
        content = file_header + times_str + steps_str + coverages_str + possible_types_str
        with open(self.__filename, "w") as f:
            f.write(content)

        # Info output.
        msg = "coverages informations are written to {}".format(self.__filename)
        self.__logger.info(msg)

        return


class FrequencyAnalysis(KMCAnalysisPlugin):
    """
    KMC plugin to do On-The-Fly process occurence frequency analysis.
    """
    def __init__(self,
                 kmc_model,
                 filename="auto_frequency.py",
                 buffer_size=500):
        """
        Constructor of TOFAnalysis object.

        Parameters:
        -----------
        kmc_model: KMC model object of Kynetix.KineticModel.

        filename: The name of data file, str.

        buffer_size: The max length of recorder variables.
        """
        # LatticeModel object.
        self.__kmc_model = kmc_model

        # Recorder variables.
        self.__times = []
        self.__steps = []

        # Process indices.
        self.__picked_indices = []

        # Process pick statistics list.
        nprocess = len(kmc_model.processes())
        self.__process_occurencies = [0]*nprocess

        # Set logger.
        self.__logger = logging.getLogger("model.solvers.KMCSolver.FrequencyAnalysis")

        # Max length of recorder variables.
        self.__buffer_size = buffer_size

        # Name of data file.
        self.__filename = filename

    def setup(self, step, time, configuration, interactions):
        # Append time and step.
        self.__times.append(time)
        self.__steps.append(step)

        # Flush counter.
        self.__flush_counter = 0

        # Create statistic data file.
        variables_str = ("times = []\nsteps = []\npicked_indices = []\n" +
                         "process_occurencies = []\n")
        with open(self.__filename, "w") as f:
            content = file_header + variables_str
            f.write(content)

    def registerStep(self, step, time, configuration, interactions):
        # Append time and step.
        self.__times.append(time)
        self.__steps.append(step)

        # Append picked index.
        picked_index = interactions.pickedIndex()
        self.__picked_indices.append(picked_index)

        # Add to collection list.
        self.__process_occurencies[picked_index] += 1

        # Check and flush.
        if len(self.__picked_indices) > self.__buffer_size:
            self.__flush()

    def finalize(self):
        """
        Write all data to files.
        """
        # Flush data left to file.
        self.__flush()

        # Write process occurencies to file.
        occurencies_str = get_list_string("process_occurencies",
                                          self.__process_occurencies)
        with open(self.__filename, "a") as f:
            f.write(occurencies_str)

        # Info output.
        msg = "All frequency informations are written to {}".format(self.__filename)
        self.__logger.info(msg)

    def __flush(self):
        """
        Private function to flush data in buffer.
        """
        # Get picked index extension strings.
        var_name = "picked_indices_{}".format(self.__flush_counter)
        list_str = get_list_string(var_name, self.__picked_indices, ncols=40)
        extend_str = "picked_indices.extend({})\n\n".format(var_name)
        picked_indices_str = list_str + extend_str

        # Get times extension strings.
        var_name = "times_{}".format(self.__flush_counter)
        list_str = get_list_string(var_name, self.__times)
        extend_str = "times.extend({})\n\n".format(var_name)
        times_str = list_str + extend_str

        # Get steps extension strings.
        var_name = "steps_{}".format(self.__flush_counter)
        list_str = get_list_string(var_name, self.__steps, ncols=10)
        extend_str = "steps.extend({})\n\n".format(var_name)
        steps_str = list_str + extend_str

        # Get all content to be flushed.
        comment = ("\n# -------------------- flush {} " +
                   "---------------------\n").format(self.__flush_counter)
        content = comment + times_str + steps_str + picked_indices_str

        # Flush to file.
        with open(self.__filename, "a") as f:
            f.write(content)

        # Free buffers.
        self.__picked_indices = []
        self.__steps = []
        self.__times = []

        self.__flush_counter += 1


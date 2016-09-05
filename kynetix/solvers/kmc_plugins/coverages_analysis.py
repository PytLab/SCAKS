from copy import deepcopy
import logging

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
    from kynetix.solvers.kmc_plugins.plugin_backends.kmc_functions import *
except ImportError:
    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    print "!!!   WARNING: plugin backends extension not found.   !!!"
    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    from kynetix.solvers.kmc_plugins.kmc_functions import *

from kynetix import file_header
from kynetix import mpi_master
from kynetix.utilities.format_utilities import get_list_string


class CoveragesAnalysis(KMCAnalysisPlugin):
    """
    KMC plugin to do On-The-Fly coverage analysis.
    """
    # {{{
    def __init__(self, kmc_model,
                 filename="auto_coverages.py",
                 buffer_size=500):
        """
        Constructor of CoverageAnalysis object.

        Parameters:
        -----------
        kmc_model: KMC model object of Kynetix.KineticModel.

        filename: The name of data file, str.

        buffer_size: The max length of recorder variables.
        """
        # LatticeModel object.
        self.__kmc_model = kmc_model

        self.__coverage_ratios = kmc_model.coverage_ratios()

        # Recorder variables.
        self.__times = []
        self.__steps = []
        self.__coverages = []

        self.__possible_types = kmc_model.possible_element_types()

        # Set logger.
        if mpi_master:
            self.__logger = logging.getLogger("model.solvers.KMCSolver.CoveragesAnalysis")

        # Set data file name.
        self.__filename = filename

        # Data flush variables.
        self.__flush_counter = 0
        self.__buffer_size = buffer_size

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

        if mpi_master:
            # Create data file.
            times_str = "times = []\n"
            steps_str = "steps = []\n"
            coverages_str = "coverages = {}\n".format([[]]*len(self.__possible_types))
            possible_types_str = get_list_string("possible_types", self.__possible_types)
            with open(self.__filename, "w") as f:
                content = (file_header + times_str + steps_str +
                           coverages_str + "\n" + possible_types_str)
                f.write(content)

    def registerStep(self, step, time, configuration, interactions):
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

        buffer_full = len(self.__coverages) >= self.__buffer_size
        if mpi_master and buffer_full:
            self.__flush()

    def finalize(self):
        """
        Write all data to files.
        """
        if mpi_master:
            # Write to file.
            self.__flush()

            # Info output.
            msg = "coverages informations are written to {}".format(self.__filename)
            self.__logger.info(msg)

    def __flush(self):
        """
        Private helper function to flush data in buffer.
        """
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

        # Get coverages extension strings.
        coverages = zip(*self.__coverages)
        var_name = "coverages_{}".format(self.__flush_counter)
        list_str = get_list_string(var_name, coverages, ncols=1)

        extend_str = ""
        for idx in xrange(len(self.__possible_types)):
            sub_extend_str = "coverages[{}].extend(coverages_{}[{}])\n"
            sub_extend_str = sub_extend_str.format(idx, self.__flush_counter, idx)
            extend_str += sub_extend_str

        coverages_str = list_str + extend_str

        # Get all content to be flushed.
        comment = ("\n# -------------------- flush {} " +
                   "---------------------\n").format(self.__flush_counter)
        content = comment + times_str + steps_str + coverages_str

        # Flush to file.
        with open(self.__filename, "a") as f:
            f.write(content)

        # Free buffers.
        self.__coverages = []
        self.__steps = []
        self.__times = []

        self.__flush_counter += 1
    # }}}


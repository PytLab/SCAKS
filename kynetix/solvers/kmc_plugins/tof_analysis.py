import logging
from operator import mul

try:
    from KMCLib import KMCAnalysisPlugin
except ImportError:
    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    print "!!!                                                    !!!"
    print "!!!         WARNING: KMCLibX is not installed          !!!"
    print "!!! Any kMC calculation using KMCLibX will be disabled !!!"
    print "!!!                                                    !!!"
    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

from kynetix import file_header
from kynetix.mpicommons import mpi
from kynetix.utilities.format_utilities import get_list_string


class TOFAnalysis(KMCAnalysisPlugin):
    """
    KMC plugin to do On-The-Fly process instantaneous tof analysis.
    """
    # {{{
    def __init__(self,
                 kmc_model,
                 filename="auto_tofs.py",
                 buffer_size=500,
                 tof_start=0):
        """
        Constructor of TOFAnalysis object.

        Parameters:
        -----------
        kmc_model: KMC model object of Kynetix.KineticModel.

        filename: The name of data file, str.

        buffer_size: The max length of recorder variables.
        """
        # {{{
        # LatticeModel object.
        self.__kmc_model = kmc_model

        # Recorder variables.
        self.__times = []
        self.__steps = []

        # variables for instantaneous tof collection.
        nprocess = len(kmc_model.processes)
        self.__occurencies = [0]*nprocess

        self.__start_time = 0.0
        self.__tofs, self.__times = [], []

        self.__tof_interval = self.__kmc_model.tof_interval

        # Set logger.
        if mpi.is_master:
            self.__logger = logging.getLogger("model.solvers.KMCSolver.TOFAnalysis")

        # Max length of recorder variables.
        self.__flush_counter = 0
        self.__buffer_size = buffer_size

        # Name of data file.
        self.__filename = filename

        # Number of active sites.
        repetitions = self.__kmc_model.repetitions
        basis_sites = self.__kmc_model.basis_sites
        self.__nsites = reduce(mul, repetitions)*len(basis_sites)
        # }}}

    def setup(self, step, time, configuration, interactions):
        # Create statistic data file.
        variables_str = ("times = []\ntofs = []\n")
        process_mapping = self.__kmc_model.process_mapping
        process_mapping_str = get_list_string("processes", process_mapping, 1)

        with open(self.__filename, "w") as f:
            content = file_header + variables_str + process_mapping_str
            f.write(content)

    def registerStep(self, step, time, configuration, interactions):
        # Get picked process index.
        picked_index = interactions.pickedIndex()
        self.__occurencies[picked_index] += 1

        # Check time.
        delta_t = time - self.__start_time
        if delta_t > self.__tof_interval:
            # Calculate tof for all processes.
            tof = [occurency/(self.__nsites*delta_t)
                   for occurency in self.__occurencies]
            self.__tofs.append(tof)
            self.__times.append(self.__start_time)  # use start_time here.

            # Refresh.
            self.__start_time = time
            nprocess = len(self.__kmc_model.processes)
            self.__occurencies = [0]*nprocess

            # Flush buffer.
            buffer_full = (len(self.__times) >= self.__buffer_size)
            if buffer_full and mpi.is_master:
                self.__flush()

    def finalize(self):
        """
        Write all data to files.
        """
        if mpi.is_master and (len(self.__times) > 0):
            self.__flush()

            # Info output.
            msg = "Instantaneos TOF informations are written to {}".format(self.__filename)
            self.__logger.info(msg)

    def __flush(self):
        """
        Private function to flush data in buffer.
        """
        # Get times extension strings.
        var_name = "times_{}".format(self.__flush_counter)
        list_str = get_list_string(var_name, self.__times)
        extend_str = "times.extend({})\n\n".format(var_name)
        times_str = list_str + extend_str

        # Get tofs extension strings.
        var_name = "tofs_{}".format(self.__flush_counter)
        list_str = get_list_string(var_name, self.__tofs, ncols=1)
        extend_str = "tofs.extend({})\n\n".format(var_name)
        tofs_str = list_str + extend_str

        # Get all content to be flushed.
        comment = ("\n# -------------------- flush {} " +
                   "---------------------\n").format(self.__flush_counter)
        content = comment + times_str + tofs_str

        # Flush to file.
        with open(self.__filename, "a") as f:
            f.write(content)

        # Free buffers.
        self.__times = []
        self.__tofs = []

        self.__flush_counter += 1

    # }}}


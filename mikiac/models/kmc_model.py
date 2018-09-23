import logging
import os

from ..mpicommons import mpi
from .kinetic_model import KineticModel
from ..descriptors.descriptors import *


class KMCModel(KineticModel):

    # {{{
    # Basis vectors of unit cell.
    cell_vectors = SpaceVectors("cell_vectors",
                                   default=[[1.0, 0.0, 0.0],
                                            [0.0, 1.0, 0.0],
                                            [0.0, 0.0, 1.0]])

    # Basis sites coordinates.
    basis_sites = SpaceVectors("basis_sites",
                                  default=[[0.0, 0.0, 0.0]])

    # Supercell repetitions.
    repetitions = Sequence("repetitions",
                              default=(1, 1, 1),
                              entry_type=int)

    # POC.
    periodic = Sequence("periodic",
                           default=(True, True, True),
                           entry_type=bool)

    # kMC step number.
    nstep = Integer("nstep", default=1)

    # Random seed for kMC simulation.
    random_seed = Integer("random_seed", default=None)

    # Flag for if adding time in seed number.
    time_seed = Bool('time_seed', default=False)

    # Interval for trajectory dumping.
    trajectory_dump_interval = Integer("trajectory_dump_interval",
                                          default=1)

    # Random generator type.
    random_generator = String("random_generator",
                                 default="MT",
                                 candidates=["MT", "MINSTD", "RANLUX24", "RANLUNX48"])

    # kMC On-the-fly analysis type.
    analysis = Sequence("analysis",
                           default=[],
                           entry_type=str,
                           candidates=["CoveragesAnalysis",
                                       "FrequencyAnalysis",
                                       "TOFAnalysis",
                                       "EventAnalysis"])

    # Interval of doing on-the-fly analysis.
    analysis_interval = AnalysisInterval("analysis_interval", default=None)

    # All possible element types.
    possible_element_types = Sequence("possible_element_types",
                                         default=[],
                                         entry_type=str)

    # All possible site types.
    possible_site_types = Sequence("possible_site_types",
                                      default=[],
                                      entry_type=str)

    # Empty type.
    empty_type = String("empty_type", default="V")

    # Step from which TOF statistic begins.
    tof_start = Integer("tof_start", default=0)

    # Time limit.
    time_limit = Float("time_limit", default=float("inf"))

    # Coverage ratios.
    coverage_ratios = Sequence("coverage_ratios",
                                  default=[],
                                  entry_type=float)

    # Extra trajectory dump control range.
    extra_trajectories = Sequence("extra_trajectries",
                                     default=None,
                                     entry_type=int)

    # The time kMC simulation start.
    start_time = Float("start_time", default=0.0)

    # Interval for instantaneous TOF calculation.
    tof_interval = Float("tof_interval", default=10)

    # Flag for redistribution operation.
    do_redistribution = Bool("do_redistribution", default=False)

    # Interval for redistribution operation.
    redistribution_interval = Integer("redistribution_interval",
                                         default=1)

    # Default fast species.
    fast_species = Sequence("fast_species", default=None, entry_type=str)

    # Split number for constrained redistribution.
    nsplits = Sequence("nsplits", default=(1, 1, 1), entry_type=int)

    # Distributor type.
    distributor_type = String("distributor_type",
                                 default="RandomDistributor",
                                 candidates=["RandomDistributor", "ProcessRandomDistributor"])
    # }}}

    def __init__(self, **kwargs):
        """
        Parameters:
        -----------
        setup_file: kinetic model set up file, str.

        setup_dict: A dictionary contains essential setup parameters for kinetic model.
        
        logger_level: logging level, int.

        file_handler_level: logging level for file handler, int.

        console_handler_level: logging level for console handler, int.

        Example:
        --------
        >>> from .models.kinetic_model import MicroKineticModel
        >>> modelodel(setup_file="setup.mkm",
                      logger_level=logging.WARNING)
        """
        super(KMCModel, self).__init__(**kwargs)

    # Overwrite father's function.
    def _set_logger(self, filename="out.log"):
        super(KMCModel, self)._set_logger(filename)

    def run(self, scripting=True, trajectory_type="lattice"):
        """
        Function to do kinetic Monte Carlo simulation.

        Parameters:
        -----------
        scripting: generate lattice script or not, True by default, bool.

        trajectory_type: The type of trajectory to use, the default type is "lattice", str.
                         "xyz" | "lattice". 

        """
        # Get processes.
        self.__processes = self.solver.processes
        self.__process_mapping = self.solver.process_mapping

        # Run the lattice model.
        self.__solver.run(scripting=scripting, trajectory_type=trajectory_type)

    @Property
    def log_allowed(self):
        """
        Flag for if log output is allowed.
        """
        # Only master processor can output log.
        return True if mpi.is_master else False

    @Property
    def process_dicts(self):
        """
        Query function for process dicts list.
        """
        return self.__process_dicts

    @Property
    def processes(self):
        """
        Query function for processes list.
        """
        return self.__processes

    @Property
    def configuration(self):
        """
        Query function for KMCConfiguration of model.
        """
        return self.__configuration

    @Property
    def sitesmap(self):
        """
        Query function for KMCSitesMap of model.
        """
        return self.__sitesmap

    @Property
    def process_mapping(self):
        """
        Query function for process reaction type mapping.
        """
        return self.__process_mapping


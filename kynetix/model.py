import cPickle as cpkl
import copy
import inspect
import logging
import logging.config
import os
import sys

import kynetix.descriptors.descriptors as dc
import kynetix.descriptors.component_descriptors as cpdc
from kynetix import mpi_master, mpi_size, mpi_installed
from kynetix.database.thermo_data import kB_eV, h_eV
from kynetix.errors.error import *
from kynetix.functions import *
from kynetix.utilities.profiling_utitlities import do_cprofile


class KineticModel(object):
    """
    Main class for kinetic models.
    """

    # Attribute descriptors.
    # {{{
    setup_file = dc.String("setup_file", default="")
    setup_dict = dc.Dict("setup_dict", default={})
    verbosity = dc.Integer("verbosity",
                           default=logging.INFO,
                           candidates=range(0, 60, 10))

    parser = cpdc.Component("parser", default=None, candidates=["RelativeEnergyParser",
                                                                "CsvParser",
                                                                "KMCParser"])

    solver = cpdc.Component("solver", default=None, candidates=["KMCSolver",
                                                                "SteadyStateSolver"])

    corrector = cpdc.Component("corrector", default=None, candidates=["ThermodynamicCorrector"])

    table_maker = cpdc.Component("table_maker", default=None, candidates=["CsvMaker"])

    plotter = cpdc.Component("plotter", default=None, candidates=["EnergyProfilePlotter"])

    # Temperature.
    temperature = dc.Float("temperature", default=298.0)

    # Reaction expressions.
    rxn_expressions = dc.Sequence("rxn_expressions", default=[], entry_type=str)

    # Definition dict of species.
    species_definitions = dc.Dict("species_definitions", default={})

    # Model core components.
    components = dc.Sequence("components", default=["parser"], entry_type=str)

    # Data precision.
    decimal_precision = dc.Integer("decimal_precision", default=100)

    # Perturbation size for numerical jacobian matrix.
    perturbation_size = dc.Float("perturbation_size", default=0.01)

    # Direction of perturbation.
    perturbation_direction = dc.String("perturbation_direction", 
                                       default="right",
                                       candidates=["right", "left"])

    # Archived variables.
    archived_variables = dc.Sequence("archive_data",
                                     default=["steady_state_coverages"],
                                     entry_type=str)

    # Numerical representation.
    numerical_representation = dc.String("numerical_representation",
                                         default="mpmath",
                                         candidates=["mpmath", "gmpy", "sympy"])

    # Rootfinding iterator type.
    rootfinding = dc.String("rootfinding",
                            default="MDNewton",
                            candidates=["MDNewton", "ConstrainedNewton"])

    # Iteration tolerance.
    tolerance = dc.Float("tolerance", default=1e-8)

    # Max iteraction steps.
    max_rootfinding_iterations = dc.Integer("max_rootfinding_iterations",
                                            default=100)

    # Ode integration buffer size.
    ode_buffer_size = dc.Integer("ode_buffer_size", default=500)

    # Ode ouptut interval.
    ode_output_interval = dc.Integer("ode_output_interval", default=200)

    # File to store data.
    data_file = dc.String("data_file", default="data.pkl")

    # Species used for conversion from relative energy to absolute eneergy.
    ref_species = dc.Sequence("ref_species", default=[], entry_type=str)

    # Reference energies used to calculate formation energy.
    ref_energies = dc.Dict("ref_energies", default={})

    # Basis vectors of unit cell.
    cell_vectors = dc.SpaceVectors("cell_vectors",
                                   default=[[1.0, 0.0, 0.0],
                                            [0.0, 1.0, 0.0],
                                            [0.0, 0.0, 1.0]])

    # Basis sites coordinates.
    basis_sites = dc.SpaceVectors("basis_sites",
                                  default=[[0.0, 0.0, 0.0]])

    # Area of unit cell (m^2).
    unitcell_area = dc.Float("unitcell_area", default=0.0)

    # Ratio of active area.
    active_ratio = dc.Float("active_ratio", default=1.0)

    # Supercell repetitions.
    repetitions = dc.Sequence("repetitions",
                              default=(1, 1, 1),
                              entry_type=int)

    # POC.
    periodic = dc.Sequence("periodic",
                           default=(True, True, True),
                           entry_type=bool)

    # kMC step number.
    nstep = dc.Integer("nstep", default=1)

    # Random seed for kMC simulation.
    random_seed = dc.Integer("random_seed", default=None)

    # Interval for trajectory dumping.
    trajectory_dump_interval = dc.Integer("trajectory_dump_interval",
                                          default=1)

    # Random generator type.
    random_generator = dc.String("random_generator",
                                 default="MT",
                                 candidates=["MT", "MINSTD", "RANLUX24", "RANLUNX48"])

    # kMC On-the-fly analysis type.
    analysis = dc.Sequence("analysis",
                           default=[],
                           entry_type=str,
                           candidates=["CoveragesAnalysis",
                                       "FrequencyAnalysis",
                                       "TOFAnalysis"])

    # Interval of doing on-the-fly analysis.
    analysis_interval = dc.Sequence("analysis_interval",
                                    default=None,
                                    entry_type=int)

    # All possible element types.
    possible_element_types = dc.Sequence("possible_element_types",
                                         default=[],
                                         entry_type=str)

    # All possible site types.
    possible_site_types = dc.Sequence("possible_site_types",
                                      default=[],
                                      entry_type=str)

    # Empty type.
    empty_type = dc.String("empty_type", default="V")

    # Step from which TOF statistic begins.
    tof_start = dc.Integer("tof_start", default=0)

    # Time limit.
    time_limit = dc.Float("time_limit", default=float("inf"))

    # Coverage ratios.
    coverages_ratios = dc.Sequence("coverages_ratios",
                                   default=[],
                                   entry_type=float)

    # Extra trajectory dump control range.
    extra_trajectories = dc.Sequence("extra_trajectries",
                                     default=None,
                                     entry_type=int)

    # The time kMC simulation start.
    start_time = dc.Float("start_time", default=0.0)

    # Interval for instantaneous TOF calculation.
    tof_interval = dc.Float("tof_interval", default=10)

    # Flag for redistribution operation.
    do_redistribution = dc.Bool("do_redistribution", default=False)

    # Interval for redistribution operation.
    redistribution_interval = dc.Integer("redistribution_interval",
                                         default=1)

    # Default fast species.
    fast_species = dc.Sequence("fast_species", default=None, entry_type=str)

    # Split number for constrained redistribution.
    nsplits = dc.Sequence("nsplits", default=(1, 1, 1), entry_type=int)

    # Distributor type.
    distributor_type = dc.String("distributor_type",
                                 default="RandomDistributor",
                                 candidates=["RandomDistributor", "ProcessRandomDistributor"])
    # }}}

    def __init__(self, setup_file=None,
                       setup_dict=None,
                       verbosity=logging.INFO):
        """
        Parameters:
        -----------
        setup_file: kinetic model set up file, str.

        setup_dict: A dictionary contains essential setup parameters for kinetic model.
        
        verbosity: logging level, int.

        Example:
        --------
        >>> from kynetix.model import KineticModel
        >>> model = KineticModel(setup_file="setup.mkm",
                                 verbosity=logging.WARNING)
        """

        # {{{
        self.__class_name = self.__class__.__name__

        # Physical constants.
        self.__kB = kB_eV
        self.__h = h_eV

        # Get setup dict.
        if setup_file is None and setup_dict is None:
            msg = "setup_file or setup_dict must be supplied for kinetic model construction"
            raise ValueError(msg)

        # Setup dict has higher priority than setup file.
        if setup_dict is not None:
            self.setup_dict = setup_dict
        else:
            self.setup_file = setup_file
            globs, locs = {}, {}
            execfile(self.setup_file, globs, locs)
            self.setup_dict = locs

        self.verbosity = verbosity

        # Set logger.
        self.__set_logger()

        # Output MPI info.
        if mpi_installed and mpi_master:
            self.__logger.info("------------------------------------")
            self.__logger.info(" Model is runing in MPI Environment ")
            self.__logger.info(" Number of process: {}".format(mpi_size))
            self.__logger.info("------------------------------------")
            self.__logger.info(" ")

        # Energy flags.
        self.__has_absolute_energy = False
        self.__has_relative_energy = False
        self.__relative_energies = {}

        # Load setup file.
        self.__load(self.setup_dict)
        self.__logger.info('kinetic modeling...success!\n')
        # }}}

    def run_mkm(self, **kwargs):
        """
        Function to solve Micro-kinetic model using Steady State Approxmiation
        to get steady state coverages and turnover frequencies.

        Parameters:
        -----------
        init_cvgs: Initial guess for coverages, tuple of floats.

        correct_energy: add free energy corrections to energy data or not, bool

        solve_ode: solve ODE only or not, bool

        fsolve: use scipy.optimize.fsolve to get low-precision root or not, bool

        coarse_guess: use fsolve to do initial coverages preprocessing or not, bool

        XRC: calculate degree of rate control or nor, bool.

        product_name: Production name of the model, str. e.g. "CH3OH_g"

        data_file: The name of data file, str.

        """
        # {{{
        # Setup default parameters.
        init_cvgs = kwargs.pop("init_cvgs", None)
        relative = kwargs.pop("relative", False)
        correct_energy = kwargs.pop("correct_energy", False)
        solve_ode = kwargs.pop("solve_ode", False)
        fsolve = kwargs.pop("fsolve", False)
        coarse_guess = kwargs.pop("coarse_guess", True)
        XRC = kwargs.pop("XRC", False)
        product_name = kwargs.pop("product_name", None)
        data_file = kwargs.pop("data_file", "./rel_energy.py")

        if kwargs:
            for key in kwargs:
                msg = "Found redundant keyword argument: {}".format(key)
                self.__logger.warning(msg)

        if mpi_master:
            self.__logger.info('--- Solve Micro-kinetic model ---')

        # Get parser and solver.
        parser = self.__parser
        solver = self.__solver

        # Parse data.
        if mpi_master:
            self.__logger.info('reading data...')
        if relative:
            if mpi_master:
                self.__logger.info('use relative energy directly...')
        else:
            if mpi_master:
                self.__logger.info('convert relative to absolute energy...')
        parser.parse_data(filename=data_file, relative=relative)

        # -- solve steady state coverages --
        if mpi_master:
            self.__logger.info('passing data to solver...')
        solver.get_data()

        # solve ODE
        # !! do ODE integration AFTER passing data to solver !!
        if solve_ode:
            if mpi_master:
                self.__logger.info("initial coverages = %s", str(init_cvgs))
            solver.solve_ode(initial_cvgs=init_cvgs)
            return

        # set initial guess(initial coverage)
        # if there is converged coverage in current path,
        # use it as initial guess
        if init_cvgs:
            # Check init_cvgs type.
            if not isinstance(init_cvgs, (tuple, list)):
                msg = "init_cvgs must be a list or tuple, but {} received."
                msg = msg.format(type(init_cvgs))
                raise ParameterError(msg)

            # Check coverages length.
            if len(init_cvgs) != len(self.__adsorbate_names):
                msg = "init_cvgs must have {} elements, but {} is supplied"
                msg = msg.format(len(self.__adsorbate_names), len(init_cvgs))
                raise ParameterError(msg)

            if mpi_master:
                self.__logger.info('use user-defined coverages as initial guess...')

        elif os.path.exists("./data.pkl"):
            with open('data.pkl', 'rb') as f:
                data = cpkl.load(f)
            init_guess = 'steady_state_coverage'
            if init_guess in data:
                if mpi_master:
                    self.__logger.info('use coverages in data.pkl as initial guess...')
                init_cvgs = data[init_guess]
                coarse_guess = False
            else:
                if mpi_master:
                    self.__logger.info('use Boltzmann coverages as initial guess...')
                init_cvgs = solver.boltzmann_coverages()

        else:  # use Boltzmann coverage
            if mpi_master:
                self.__logger.info('use Boltzmann coverages as initial guess...')
            init_cvgs = solver.boltzmann_coverages()

        # Solve steady state coverages.
        # Use scipy.optimize.fsolve or not (fast but low-precision).
        if fsolve:
            if mpi_master:
                self.__logger.info('using fsolve to get steady state coverages...')
            ss_cvgs = solver.fsolve_steady_state_cvgs(init_cvgs)
        else:
            if coarse_guess:
                if mpi_master:
                    self.__logger.info('getting coarse steady state coverages...')
                init_cvgs = solver.coarse_steady_state_cvgs(init_cvgs)  # coarse root
            if mpi_master:
                self.__logger.info('getting precise steady state coverages...')
            ss_cvgs = solver.get_steady_state_cvgs(init_cvgs)

        # Get TOFs for gases.
        tofs = solver.get_tof(ss_cvgs)

        # Get reversibilities.
        rf, rr = solver.get_rates(ss_cvgs)
        reversibilities = solver.get_reversibilities(rf, rr)

        # Calculate XRC.
        if XRC:
            if product_name is None:
                raise ParameterError("production name must be provided to get XRC.")
            solver.get_single_XRC(product_name, epsilon=1e-5)

        return
        # }}}

    def run_kmc(self,
                processes_file=None,
                configuration_file=None,
                sitesmap_file=None,
                scripting=True,
                trajectory_type="lattice"):
        """
        Function to do kinetic Monte Carlo simulation.

        Parameters:
        -----------
        processes_file: The name of processes definition file, str.
                        the default name is "kmc_processes.py".

        configuration_file: The name of configuration definition file, str.
                            the default name is "kmc_processes.py".

        sitesmap_file: The name of sitesmap definition file, str.
                       the default name is "kmc_processes.py".

        scripting: generate lattice script or not, True by default, bool.

        trajectory_type: The type of trajectory to use, the default type is "lattice", str.
                         "xyz" | "lattice". 

        """
        parser = self.__parser

        # Parse processes, configuration, sitesmap.
        self.__processes = parser.parse_processes(filename=processes_file)
        self.__configuration = parser.parse_configuration(filename=configuration_file)
        self.__sitesmap = parser.construct_sitesmap(filename=sitesmap_file)

        # Set process reaction mapping.
        self.__process_mapping = parser.process_mapping()

        # Run the lattice model.
        self.__solver.run(scripting=scripting,
                          trajectory_type=trajectory_type)

    def __set_logger(self):
        """
        Private function to get logging.logger instance as logger of kinetic model.
        """
        # {{{
        logger = logging.getLogger('model')
        if os.path.exists('./logging.conf'):
            logging.config.fileConfig('./logging.conf')
        else:
            # Set logging level.
            if hasattr(self, "_" + self.__class_name + "__verbosity"):
                logger.setLevel(self.__verbosity)
            else:
                logger.setLevel(logging.INFO)

            # Create handlers.
            std_hdlr = logging.FileHandler('out.log')
            std_hdlr.setLevel(logging.DEBUG)
            console_hdlr = logging.StreamHandler()
            console_hdlr.setLevel(logging.INFO)

            # Create formatter and add it to the handlers.
            formatter = logging.Formatter('%(name)s   %(levelname)-8s %(message)s')
            std_hdlr.setFormatter(formatter)
            console_hdlr.setFormatter(formatter)

            # Add the handlers to logger.
            logger.addHandler(std_hdlr)
            logger.addHandler(console_hdlr)

        self.__logger = logger
        # }}}

    def set_logger_level(self, handler_type, level):
        """
        Set the logging level of logger handler.

        Parameters:
        -----------
        handler_type: logger handler name, str.

        level: logging level, int.

        Returns
        -------
        Old logging level of the handler, int.
        """
        # Locate handler.
        handler = None
        for h in self.__logger.handlers:
            if h.__class__.__name__ == handler_type:
                handler = h
                break

        if handler is None:
            raise ValueError("Unknown handler type '{}'".format(handler_type))

        # Reset logging level.
        old_level = handler.level
        handler.setLevel(level)

        return old_level

    def __load(self, setup_dict):
        """
        Load 'setup_file' into kinetic model by exec setup file
        and assigning all local variables as attrs of model.
        For tools, create the instances of tool classes and
        assign them as the attrs of model.
        """
        # {{{
        if mpi_master:
            self.__logger.info('Loading Kinetic Model...\n')
            self.__logger.info('read in parameters...')

        setup_dict_copy = copy.deepcopy(setup_dict)

        # Set model attributes in setup file.
        for key, value in setup_dict.iteritems():
            # Parser & solver will be set later.
            if key in ["parser", "solver"]:
                continue

            # Set parameters in setup dict as attiributes of model.
            setattr(self, key, value)

            # Output info.
            specials = ("rxn_expressions", "species_definitions")
            if key not in specials:
                if mpi_master:
                    self.__logger.info('{} = {}'.format(key, str(value)))

            # If it is a iterable, loop to output.
            else:
                if mpi_master:
                    self.__logger.info("{} =".format(key))
                if type(setup_dict[key]) is dict:
                    for k, v in value.iteritems():
                        if mpi_master:
                            self.__logger.info("        {}: {}".format(k, v))
                else:
                    for item in value:
                        if mpi_master:
                            self.__logger.info("        {}".format(item))

            # Delete.
            del setup_dict_copy[key]

        # Instantialize parser.
        self.parser = setup_dict["parser"]
        del setup_dict_copy["parser"]

        # use parser parse essential attrs for other tools
        # Parse elementary rxns
        if mpi_master:
            self.__logger.info('Parsing elementary rxns...')
        if self.__rxn_expressions:
            (self.__adsorbate_names,
             self.__gas_names,
             self.__liquid_names,
             self.__site_names,
             self.__transition_state_names,
             self.__elementary_rxns_list) = \
                self.__parser.parse_elementary_rxns(self.__rxn_expressions)

        # Instantialize solver.
        if "solver" in setup_dict:
            self.solver = setup_dict["solver"]
            del setup_dict_copy["solver"]

        # Check if there is redundant parameters.
        if setup_dict_copy:
            for key in setup_dict_copy:
                msg = "Found redundant parameter '{}'".format(key)
                self.__logger.warning(msg)
        # }}}

    @dc.Property
    def kB(self):
        return self.__kB

    @dc.Property
    def h(self):
        return self.__h

    @dc.Property
    def logger(self):
        """
        Query function for model logger.
        """
        return self.__logger

    @dc.Property
    def elementary_rxns_list(self):
        """
        Query function for elementary reactions list.
        """
        return self.__elementary_rxns_list

    @dc.Property
    def site_names(self):
        """
        Query function for site names in model.
        """
        return self.__site_names

    @dc.Property
    def adsorbate_names(self):
        """
        Query function for adsorbate names in model.
        """
        return self.__adsorbate_names

    @dc.Property
    def gas_names(self):
        """
        Query function for gas names in model.
        """
        return self.__gas_names

    @dc.Property
    def liquid_names(self):
        """
        Query function for liquid names in model.
        """
        return self.__liquid_names

    @dc.Property
    def transition_state_names(self):
        """
        Query function for transition state species names in model.
        """
        return self.__transition_state_names

    @dc.Property
    def has_relative_energy(self):
        """
        Query function for relative energy flag.
        """
        return self.__has_relative_energy

    @dc.Property
    def has_absolute_energy(self):
        """
        Query function for absolute energy flag.
        """
        return self.__has_absolute_energy

    @dc.Property
    @return_deepcopy
    def relative_energies(self):
        """
        Query function for relative energy in data file.
        """
        return self.__relative_energies

    # ------------------------------------
    # KMC Parameters query functions.

    @dc.Property
    def processes(self):
        """
        Query function for processes list.
        """
        return self.__processes

    @dc.Property
    def configuration(self):
        """
        Query function for KMCConfiguration of model.
        """
        return self.__configuration

    @dc.Property
    def sitesmap(self):
        """
        Query function for KMCSitesMap of model.
        """
        return self.__sitesmap

    @dc.Property
    def process_mapping(self):
        """
        Query function for process reaction type mapping.
        """
        return self.__process_mapping


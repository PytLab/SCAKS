import cPickle as cpkl
import inspect
import logging
import logging.config
import os
import sys

from kynetix import mpi_master, mpi_size, mpi_installed
from kynetix.database.thermo_data import kB_eV, h_eV
from kynetix.errors.error import *
from kynetix.functions import *
from kynetix.utilities.check_utilities import *


class KineticModel(object):
    """
    Main class for kinetic models.
    """
    def __init__(self, **kwargs):
        """
        Parameters:
        -----------
        setup_file: kinetic model set up file, str.
        
        verbosity: logging level, int.

        Example:
        --------
        >>> from kynetix.model import KineticModel
        >>> model = KineticModel(setup_file="setup.mkm",
                                 verbosity=logging.WARNING)
        """

        # Get class name.
        self.__class_name = self.__class__.__name__

        # Set physical constants.
        self.__kB = kB_eV  # Boltzmann constant from NIST, eV/K
        self.__h = h_eV    # Planck constant from NIST, eV s

        # Parse in keyword args.
        for key in kwargs:
            setattr(self, "_" + self.__class_name + "__" + key, kwargs[key])

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
        if hasattr(self, '_' + self.__class_name + '__setup_file'):
            if mpi_master:
                self.__logger.info('setup file [ {} ] is found'.format(self.__setup_file))
            model_name = self.__setup_file.rsplit('.', 1)[0]
            self.__model_name = model_name
            self.__load(self.__setup_file)
            if mpi_master:
                self.__logger.info('kinetic modeling...success!\n')
        else:
            if mpi_master:
                self.__logger.warning('setup file not read...')

    def __check_inputs(self, inputs_dict):
        """
        Private helper function to check all parameters in setup file.

        Parameters:
        -----------
        inputs_dict: A dict store all setup information.

        Returns:
        --------
        inputs_dict: The valid input dict (REFERENCE of input).
        """
        # {{{
        invalid_parameters = []
        for key, value in inputs_dict.iteritems():
            # Check parameter validity.
            if key not in type_rules and mpi_master:
                msg = (("Parameter [{}] is not a valid setup parameter, " +
                        "it will be ignored.")).format(key)
                self.__logger.warning(msg)

                # Collect the invalid parameter.
                invalid_parameters.append(key)
                continue

            rule = type_rules[key]
            if len(rule) == 1:
                # If it is check function.
                if hasattr(rule[0], "__call__"):
                    check_func = rule[0]
                    check_func(value)
                # If it is a type.
                elif type(rule[0]) is type:
                    if not isinstance(value, rule[0]):
                        msg = "{} should be a {}".format(key, rule[0])
                        raise SetupError(msg)
            else:  # Call corresponding check function.
                check_func, arg = rule
                check_func(value, arg, key)

        # Clean input dict.
        for invalid_param in invalid_parameters:
            del inputs_dict[invalid_param]

        return inputs_dict
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
        init_cvgs = setdefault_args("init_cvgs", kwargs, None)
        relative = setdefault_args("relative", kwargs, False)
        correct_energy = setdefault_args("correct_energy", kwargs, False)
        solve_ode = setdefault_args("solve_ode", kwargs, False)
        fsolve = setdefault_args("fsolve", kwargs, False)
        coarse_guess = setdefault_args("coarse_guess", kwargs, True)
        XRC = setdefault_args("XRC", kwargs, False)
        product_name = setdefault_args("product_name", kwargs, None)
        data_file = setdefault_args("data_file", kwargs, "./rel_energy.py")

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

    def __set_parser(self, parser_name):
        """
        Private function to import parser and
        set the instance of it as attr of model
        """
        # https://docs.python.org/2/library/functions.html#__import__
        basepath = os.path.dirname(
            inspect.getfile(inspect.currentframe()))
        if basepath not in sys.path:
            sys.path.append(basepath)
        # from  loggers import logger
        _module = __import__('parsers', globals(), locals())
        parser_instance = getattr(_module, parser_name)(owner=self)
        setattr(self, "_" + self.__class_name + '__parser', parser_instance)
        if mpi_master:
            self.__logger.info('parser is set.')

    def __set_logger(self):
        """
        Private function to get logging.logger instance as logger of kinetic model.
        """
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

    def __load(self, setup_file):
        """
        Load 'setup_file' into kinetic model by exec setup file
        and assigning all local variables as attrs of model.
        For tools, create the instances of tool classes and
        assign them as the attrs of model.
        """
        if mpi_master:
            self.__logger.info('Loading Kinetic Model...\n')

        defaults = dict(
            data_file='data.pkl',
            grid_type='square',
            decimal_precision=100,
            tools=['parser'],
            parser='RelativeEnergyParser',
            table_maker='CsvMaker',
            solver='SteadyStateSolver',
            corrector='ThermodynamicCorrector',
            plotter='EnergyProfilePlotter',
        )

        if mpi_master:
            self.__logger.info('read in parameters...')

        # Exec setup file set local variables as attrs of model
        globs = {}
        locs = defaults
        execfile(setup_file, globs, locs)

        # Check parameters validity.
        if mpi_master:
            self.__logger.info("Check setup file validity...")
        self.__check_inputs(locs)

        # Customize model tools
        if 'tools' in locs:
            if 'parser' not in locs['tools']:
                raise ParameterError('[ parser ] must be in tools.')
            self.__tools = locs['tools']
            del locs['tools']

            if mpi_master:
                self.__logger.info('tools = {}'.format(str(self.__tools)))

        # Set model attributes in setup file.
        for key in locs.keys():
            # ignore tools which will be loaded later
            if key in self.__tools:
                continue
            setattr(self, "_" + self.__class_name + "__" + key, locs[key])

            # Output info.
            specials = ("rxn_expressions", "species_definitions")
            if key not in specials:
                if mpi_master:
                    self.__logger.info('{} = {}'.format(key, str(locs[key])))

            # If it is a iterable, loop to output.
            else:
                if mpi_master:
                    self.__logger.info("{} =".format(key))
                if type(locs[key]) is dict:
                    for k, v in locs[key].iteritems():
                        if mpi_master:
                            self.__logger.info("        {}: {}".format(k, v))
                else:
                    for item in locs[key]:
                        if mpi_master:
                            self.__logger.info("        {}".format(item))

        # assign parser ahead to provide essential attrs for other tools
        if mpi_master:
            self.__logger.info('instantiate {}'.format(str(locs['parser'])))
        self.__set_parser(locs['parser'])

        # if parser is kmc_parser use kmc_solver correspondingly
        if locs['parser'] == 'KMCParser':
            locs['solver'] = 'KMCSolver'
            if mpi_master:
                self.__logger.info('set solver [ KMCSolver ].')

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

        # load tools of model
        if mpi_master:
            self.__logger.info('instantiate model tools...')
        for key in self.__tools:
            # Auto-import classes.
            if key == 'parser':  # ignore parser which is loaded before
                continue
            if locs[key]:
                if not key.endswith('s'):
                    pyfile = key + 's'
                else:
                    pyfile = key
                basepath = os.path.dirname(inspect.getfile(inspect.currentframe()))
                if basepath not in sys.path:
                    sys.path.append(basepath)
                sublocs = {}
                _temp = __import__(pyfile, globals(), sublocs, [locs[key]])
                tool_instance = getattr(_temp, locs[key])(owner=self)
                setattr(self, "_" + self.__class_name + "__" + key, tool_instance)
                if mpi_master:
                    self.__logger.info('{} = {}'.format(key, locs[key]))
            else:
                setattr(self, "_" + self.__class_name + "__" + key, None)
                if mpi_master:
                    self.__logger.warning('{} is set to None.'.format(key))

        # Set kMC parameters.
        if (not isinstance(self.__solver, str)) and (locs["solver"] == "KMCSolver"):
            self.__set_kmc_parameters(locs)

    def __set_kmc_parameters(self, locs):
        """
        Private helper function to set KMC related parameters.
        """
        # kMC parameters check
        solver_type = repr(self.__solver).split('.')[2]
        # pseudo random generator
        if solver_type == 'kmc_solver' and 'random_generator' not in locs:
            if mpi_master:
                self.__logger.info('pseudo random generator type was not set.')
                self.__logger.info('use Mersenne-Twister by default.')

        # kMC control parameters
        if 'kmc_continue' in locs and locs['kmc_continue']:
            # read in step and time info
            if mpi_master:
                self.__logger.info('reading in kMC control parameters...')

            # read from auto_coverages.py
            if mpi_master:
                self.__logger.info('reading auto_coverages.py...')
            has_auto_coverages = os.path.exists('./auto_coverages.py')
            if has_auto_coverages:
                #raise FilesError('No auto_coverages.py in current path.')
                cvg_globs, cvg_locs = {}, {}
                execfile('./auto_coverages.py', cvg_globs, cvg_locs)
                cvg_step, cvg_time = cvg_locs['steps'][-1], cvg_locs['times'][-1]
                if mpi_master:
                    self.__logger.info('auto_coverages.py was read.')
            else:
                if mpi_master:
                    self.__logger.info('auto_coverages.py not read.')

            # read from auto_TOFs.py
            if mpi_master:
                self.__logger.info('reading auto_TOFs.py...')
            if not os.path.exists('./auto_TOFs.py'):
                raise FilesError('No auto_TOFs.py in current path.')
            tof_globs, tof_locs = {}, {}
            execfile('./auto_TOFs.py', tof_globs, tof_locs)
            tof_step, tof_time = tof_locs['steps'][-1], tof_locs['times'][-1]
            # get total_rates calculated from last loop
            total_rates = tof_locs['total_rates']
            # get time when the first TOF analysis start
            tof_start_time = tof_locs['start_time']
            if mpi_master:
                self.__logger.info('auto_TOFs.py was read.')

            # read from auto_last_types.py
            if mpi_master:
                self.__logger.info('reading auto_last_types.py...')
            if not os.path.exists('auto_last_types.py'):
                raise FilesError('No auto_last_types.py in current path.')
            types_globs, types_locs = {}, {}
            execfile('auto_last_types.py', types_globs, types_locs)
            last_types = types_locs['types']
            if mpi_master:
                self.__logger.info('auto_last_types.py was read.')

            # check data consistency
            if has_auto_coverages:
                if not cvg_step == tof_step:
                    raise FilesError('Steps(%d, %d) are not consistent.' %
                                     (cvg_step, tof_step))
                if not (abs(cvg_time - tof_time) < 1e-5):
                    raise FilesError('Times(%s, %s) are not consistent.' %
                                     (cvg_time, tof_time))

            # set model attributes
            # -----------------------------------------------------------------
            # **start_time** is the time when all analysis object end
            # in last kmc loop which is the start point in current kmc loop.
            #
            # **tof_start_time** is the time when true TOF analysis
            # begins(no TOF calculation at this point) in last kmc loop
            # which is also used in continous kmc loop.
            # -----------------------------------------------------------------
            self.__start_step, self.__start_time = tof_step, tof_time
            self.__tof_start_time = tof_start_time
            self.__start_types = last_types
            self.__total_rates = total_rates
            if mpi_master:
                self.__logger.info('kMC analysis starts from step = %d, time = %e s',
                                   self.__start_step, self.__start_time)

    def setup_file(self):
        """
        Query function for setup file.
        """
        return self.__setup_file

    def name(self):
        """
        Query function for model name.
        """
        return self.__model_name

    def kB(self):
        """
        Query function for Boltzmann constant.
        """
        return self.__kB

    def h(self):
        """
        Query function for Plank constant.
        """
        return self.__h

    def logger(self):
        """
        Query function for model logger.
        """
        return self.__logger

    @return_deepcopy
    def tools(self):
        """
        Query function for model tools.
        """
        return self.__tools

    def parser(self):
        """
        Query function for model parser object.
        """
        return self.__parser

    def solver(self):
        """
        Query function for model solver object.
        """
        return self.__solver

    def corrector(self):
        """
        Query function for model corrector.
        """
        return self.__corrector

    def plotter(self):
        """
        Query function for model plotter.
        """
        return self.__plotter

    @return_deepcopy
    def rxn_expressions(self):
        """
        Query function for reaction expressions in model.
        """
        return self.__rxn_expressions

    @return_deepcopy
    def species_definitions(self):
        """
        Query function for species definitions in model.
        """
        return self.__species_definitions

    @return_deepcopy
    def decimal_precision(self):
        """
        Query function for data precision.
        """
        return self.__decimal_precision

    def elementary_rxns_list(self):
        """
        Query function for elementary reactions list.
        """
        return self.__elementary_rxns_list

    def temperature(self):
        """
        Query function for system temperature.
        """
        return self.__temperature

    @return_deepcopy
    def site_names(self):
        """
        Query function for site names in model.
        """
        return self.__site_names

    @return_deepcopy
    def adsorbate_names(self):
        """
        Query function for adsorbate names in model.
        """
        return self.__adsorbate_names

    @return_deepcopy
    def gas_names(self):
        """
        Query function for gas names in model.
        """
        return self.__gas_names

    @return_deepcopy
    def liquid_names(self):
        """
        Query function for liquid names in model.
        """
        return self.__liquid_names

    @return_deepcopy
    def transition_state_names(self):
        """
        Query function for transition state species names in model.
        """
        return self.__transition_state_names

    def gas_thermo_mode(self):
        """
        Query function for gas mode in model.
        """
        return self.__gas_thermo_mode

    @return_deepcopy
    def ref_species(self):
        """
        Query function for reference species name.
        """
        return self.__ref_species

    def surface_name(self):
        """
        Query function for surface name.
        """
        return self.__surface_name

    def verbosity(self):
        """
        Query function for logging verbosity.
        """
        return self.__verbosity

    def has_relative_energy(self):
        """
        Query function for relative energy flag.
        """
        return self.__has_relative_energy

    def has_absolute_energy(self):
        """
        Query function for absolute energy flag.
        """
        return self.__has_absolute_energy

    @return_deepcopy
    def relative_energies(self):
        """
        Query function for relative energy in data file.
        """
        return self.__relative_energies

    def data_file(self):
        """
        Query function for data archive file name.
        """
        return self.__data_file

    def table_maker(self):
        """
        Query function for table_maker object.
        """
        return self.__table_maker

    @return_deepcopy
    def ref_energies(self):
        """
        Query function for reference energy dict.
        """
        return self.__ref_energies

    # ------------------------------------
    # KMC Parameters query functions.

    @return_deepcopy
    def cell_vectors(self):
        """
        Query function for cell base vectors.
        """
        return self.__cell_vectors

    @return_deepcopy
    def basis_sites(self):
        """
        Query function for basis sites.
        """
        return self.__basis_sites

    def unitcell_area(self):
        """
        Query function for area of unitcell.
        """
        return self.__unitcell_area

    def active_ratio(self):
        """
        Query function for active ratio(Ast/Auc).
        """
        return self.__active_ratio

    def repetitions(self):
        """
        Query function for lattice repetitions.
        """
        return self.__repetitions

    def periodic(self):
        """
        Query function for lattice periodic.
        """
        return self.__periodic

    def nstep(self):
        """
        Query function for number of kmc step.
        """
        return self.__nstep

    def random_seed(self):
        """
        Query function for random seed.
        """
        try:
            return self.__seed
        except AttributeError:
            return None

    def trajectory_dump_interval(self):
        """
        Query function for trajectory dump interval.
        """
        try:
            return self.__trajectory_dump_interval
        except AttributeError:
            return None

    def random_generator(self):
        """
        Query function for random generator name.
        """
        try:
            return self.__random_generator
        except AttributeError:
            return None

    def analysis(self):
        """
        Query function for analysis names.
        """
        return self.__analysis

    def analysis_interval(self):
        """
        Query function for analysis interval.
        """
        try:
            return self.__analysis_interval
        except AttributeError:
            return None

    def analysis_dump_interval(self):
        """
        Query function for analysis dump interval.
        """
        return self.__analysis_dump_interval

    @return_deepcopy
    def color_dict(self):
        """
        Query function for color dict for elements on surface.
        """
        return self.__color_dict

    @return_deepcopy
    def circle_attrs(self):
        """
        Query function for circle attributes for circle plotting.
        """
        return self.__circle_attrs

    def processes(self):
        """
        Query function for processes list.
        """
        return self.__processes

    def configuration(self):
        """
        Query function for KMCConfiguration of model.
        """
        return self.__configuration

    def sitesmap(self):
        """
        Query function for KMCSitesMap of model.
        """
        return self.__sitesmap

    def possible_element_types(self):
        """
        Query function for possible element types.
        """
        return self.__possible_element_types

    def empty_type(self):
        """
        Query function for empty element type.
        """
        return self.__empty_type

    def possible_site_types(self):
        """
        Query function for possible site types.
        """
        return self.__possible_site_types

    def process_mapping(self):
        """
        Query function for process reaction type mapping.
        """
        return self.__process_mapping

    def tof_start(self):
        """
        Query function for TOF collection starting step.
        """
        try:
            return self.__tof_start
        except AttributeError:
            return 0

    def time_limit(self):
        """
        Query function for KMC loop time upper bound.
        """
        try:
            return self.__time_limit
        except AttributeError:
            return float("inf")

    def coverage_ratios(self):
        """
        Query function for coverage ratios for all basis sites.
        """
        try:
            return self.__coverage_ratios
        except AttributeError:
            return [1.0]*len(self.__basis_sites)


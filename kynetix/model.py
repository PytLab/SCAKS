import os
import sys
import inspect
import logging
import logging.config
import cPickle as cpkl

from kynetix.functions import *
from kynetix.errors.error import *
from kynetix.database.thermo_data import kB_eV, h_eV


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

        # set physical constants
        self.__kB = kB_eV  # Boltzmann constant from NIST, eV/K
        self.__h = h_eV    # Planck constant from NIST, eV s

        # parse in keyword args
        for key in kwargs:
            setattr(self, "_" + self.__class_name + "__" + key, kwargs[key])

        # set logger
        self.__set_logger()

        self.has_absolute_energy = False
        self.has_relative_energy = False

        # load setup file
        if hasattr(self, '_' + self.__class_name + '__setup_file'):
            self.__logger.info('setup file [ {} ] is found'.format(self.__setup_file))
            model_name = self.__setup_file.rsplit('.', 1)[0]
            self.__model_name = model_name
            self.__load(self.__setup_file)
            self.__logger.info('kinetic modeling...success!\n')
        else:
            self.__logger.warning('setup file not read...')

    def run_mkm(self, init_cvgs=None, relative=False, correct_energy=False,
                solve_ode=False, fsolve=False, coarse_guess=True):
        '''
        Function to solve Micro-kinetic model using Steady State Approxmiation
        to get steady state coverages and turnover frequencies.

        Parameters:
        -----------
        correct_energy: add free energy corrections to energy data or not, bool

        solve_ode: solve ODE only or not, bool

        fsolve: use scipy.optimize.fsolve to get low-precision root or not, bool

        coarse_guess: use fsolve to do initial coverages preprocessing or not, bool

        '''
        self.__logger.info('--- Solve Micro-kinetic model ---')

        # parse data
        self.__logger.info('reading data...')
        self.parser.chk_data_validity()
        if relative:
            self.__logger.info('use relative energy directly...')
        else:
            self.__logger.info('convert relative to absolute energy...')
        self.parser.parse_data(relative=relative)

        # -- solve steady state coverages --
        self.__logger.info('passing data to solver...')
        self.solver.get_data()
        self.__logger.info('getting rate constants...')
        self.solver.get_rate_constants()
        self.__logger.info('getting rate expressions...')
        self.solver.get_rate_expressions(self.solver.rxns_list)

        # energy correction
        if correct_energy:
            self.__logger.info('free energy correction...')
            self.solver.correct_energies()

        # solve ODE
        # !! do ODE integration AFTER passing data to solver !!
        if solve_ode:
            self.__logger.info("initial coverages = %s", str(init_cvgs))
            self.solver.solve_ode(initial_cvgs=init_cvgs)
            return

        # set initial guess(initial coverage)
        # if there is converged coverage in current path,
        # use it as initial guess
        if init_cvgs:
            # check init_cvgs validity
            if not isinstance(init_cvgs, (tuple, list)):
                msg = ('init_cvgs must be a list or tuple, but %s supplied.' %
                       str(type(init_cvgs)))
                raise ParameterError(msg)
            if len(init_cvgs) != len(self.adsorbate_names):
                msg = ('init_cvgs must have %d elements, but %d is supplied' %
                       (len(self.adsorbate_names), len(init_cvgs)))
                raise ParameterError(msg)
            self.__logger.info('use user-defined coverages as initial guess...')

        elif os.path.exists("./data.pkl"):
            with open('data.pkl', 'rb') as f:
                data = cpkl.load(f)
            init_guess = 'steady_state_coverage'
            if init_guess in data:
                self.__logger.info('use coverages in data.pkl as initial guess...')
                init_cvgs = data[init_guess]
                coarse_guess = False
            else:
                self.__logger.info('use Boltzmann coverages as initial guess...')
                init_cvgs = self.solver.boltzmann_coverages()

        else:  # use Boltzmann coverage
            self.__logger.info('use Boltzmann coverages as initial guess...')
            init_cvgs = self.solver.boltzmann_coverages()

        # solve steady state coverages
        # use scipy.optimize.fsolve or not (fast but low-precision)
        if fsolve:
            self.__logger.info('using fsolve to get steady state coverages...')
            ss_cvgs = self.solver.fsolve_steady_state_cvgs(init_cvgs)
        else:
            if coarse_guess:
                self.__logger.info('getting coarse steady state coverages...')
                init_cvgs = self.solver.coarse_steady_state_cvgs(init_cvgs)  # coarse root
            self.__logger.info('getting precise steady state coverages...')
            ss_cvgs = self.solver.get_steady_state_cvgs(init_cvgs)

        # get TOFs for gases
        tofs = self.solver.get_cvg_tof(ss_cvgs)

        return

    def run_kmc(self):
        '''
        Function to do kinetic Monte Carlo simulation to
        get steady state coverages and turnover frequencies.
        '''
        self.solver.run()

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
            # create handlers
            std_hdlr = logging.FileHandler('out.log')
            std_hdlr.setLevel(logging.DEBUG)
            console_hdlr = logging.StreamHandler()
            console_hdlr.setLevel(logging.INFO)
            # create formatter and add it to the handlers
            formatter = logging.Formatter('%(name)s   %(levelname)-8s %(message)s')
            std_hdlr.setFormatter(formatter)
            console_hdlr.setFormatter(formatter)
            # add the handlers to logger
            logger.addHandler(std_hdlr)
            logger.addHandler(console_hdlr)

        self.__logger = logger

    def __load(self, setup_file):
        """
        Load 'setup_file' into kinetic model by exec setup file
        and assigning all local variables as attrs of model.
        For tools, create the instances of tool classes and
        assign them as the attrs of model.
        """
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

        self.__logger.info('read in parameters...')

        #exec setup file set local variables as attrs of model
        globs = {}
        locs = defaults
        execfile(setup_file, globs, locs)

        # customize model tools
        if 'tools' in locs:
            if 'parser' not in locs['tools']:
                raise ParameterError('[ parser ] must be in tools.')
            self.__tools = locs['tools']
            del locs['tools']

            self.__logger.info('tools = {}'.format(str(self.__tools)))

        # assign parser ahead to provide essential attrs for other tools
        self.__logger.info('instantiate {}'.format(str(locs['parser'])))
        self.__set_parser(locs['parser'])

        # if parser is kmc_parser use kmc_solver correspondingly
        if locs['parser'] == 'KMCParser':
            locs['solver'] = 'KMCLibSolver'
            self.__logger.info('set solver [ KMCSolver ].')

        # assign other tools
        for key in locs.keys():
            # ignore tools which will be loaded later
            if key in self.__tools:
                continue
            # check type of variables
            # Add later ...
            setattr(self, "_" + self.__class_name + "__" + key, locs[key])
            self.__logger.info('{} = {}'.format(key, str(locs[key])))

        #use parser parse essential attrs for other tools
        #parse elementary rxns
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
        self.__logger.info('instantiate model tools...')
        for key in self.__tools:
            # Auto-import classes.
            if key == 'parser':  # ignore parser which is loaded before
                continue
            try:
                if locs[key]:
                    if not key.endswith('s'):
                        pyfile = key + 's'
                    else:
                        pyfile = key
                    basepath = os.path.dirname(
                        inspect.getfile(inspect.currentframe()))
                    if basepath not in sys.path:
                        sys.path.append(basepath)
                    sublocs = {}
                    _temp = \
                        __import__(pyfile, globals(), sublocs, [locs[key]])
                    tool_instance = getattr(_temp, locs[key])(owner=self)
                    setattr(self, "_" + self.__class_name + "__" + key, tool_instance)
                    self.__logger.info('{} = {}'.format(key, locs[key]))
                else:
                    setattr(self, "_" + self.__class_name + "__" + key, None)
                    self.__logger.warning('{} is set to None.'.format(key))
            except ImportError:
                raise ToolsImportError(key.capitalize()+' '+locs[key] +
                                       ' could not be imported. ' +
                                       'Ensure that the class ' +
                                       'exists and is spelled properly.')
        # Set kMC parameters.
        if not isinstance(self.__solver, str):
            self.__set_kmc_parameters()

    def __set_kmc_parameters(self):
        """
        Private helper function to set KMC related parameters.
        """
        # kMC parameters check
        solver_type = repr(self.__solver).split('.')[2]
        # pseudo random generator
        if solver_type == 'kmc_solver' and 'random_generator' not in locs:
            self.__logger.info('pseudo random generator type was not set.')
            self.__logger.info('use Mersenne-Twister by default.')

        # kMC control parameters
        if 'kmc_continue' in locs and locs['kmc_continue']:
            # read in step and time info
            self.__logger.info('reading in kMC control parameters...')

            # read from auto_coverages.py
            self.__logger.info('reading auto_coverages.py...')
            has_auto_coverages = os.path.exists('./auto_coverages.py')
            if has_auto_coverages:
                #raise FilesError('No auto_coverages.py in current path.')
                cvg_globs, cvg_locs = {}, {}
                execfile('./auto_coverages.py', cvg_globs, cvg_locs)
                cvg_step, cvg_time = cvg_locs['steps'][-1], cvg_locs['times'][-1]
                self.__logger.info('auto_coverages.py was read.')
            else:
                self.__logger.info('auto_coverages.py not read.')

            # read from auto_TOFs.py
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
            self.__logger.info('auto_TOFs.py was read.')

            # read from auto_last_types.py
            self.__logger.info('reading auto_last_types.py...')
            if not os.path.exists('auto_last_types.py'):
                raise FilesError('No auto_last_types.py in current path.')
            types_globs, types_locs = {}, {}
            execfile('auto_last_types.py', types_globs, types_locs)
            last_types = types_locs['types']
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

    def rxn_expressions(self):
        """
        Query function for reaction expressions in model.
        """
        return self.__rxn_expressions

    def species_definitions(self):
        """
        Query function for species definitions in model.
        """
        return self.__species_definitions

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

    def site_names(self):
        """
        Query function for site names in model.
        """
        return self.__site_names

    def adsorbate_names(self):
        """
        Query function for adsorbate names in model.
        """
        return self.__adsorbate_names

    def gas_names(self):
        """
        Query function for gas names in model.
        """
        return self.__gas_names

    def liquid_names(self):
        """
        Query function for liquid names in model.
        """
        return self.__liquid_names

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

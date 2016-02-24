import os
import sys
import inspect
import re
import logging
import logging.config
import cPickle as cpkl

from functions import *
from errors.error import *
from .database.thermo_data import kB_eV, h_eV


class KineticModel(object):
    """
    Main class for kinetic models.
    """
    def __init__(self, **kwargs):

        # set physical constants
        self._kB = kB_eV  # Boltzmann constant from NIST, eV/K
        self._h = h_eV    # Planck constant from NIST, eV s

        #set kinetic attrs
        self._attr_type_dict = {
            'adsorbate_names': tuple,
            'transition_state_names': tuple,
            'gas_names': tuple,
            'descriptor_names': tuple,
            'surface_names': tuple,
            'species_definitions': dict,
            'elementary_rxns': tuple,
            'setup_file': str,
            'rxn_expression': list
        }

        # set logger
        self.set_logger()

        # parse in keyword args
        for key in kwargs:
            if key in self._attr_type_dict:
                try:
                    val = self._attr_type_dict[key](kwargs[key])
                except:
                    raise AttributeError('Argument \''+key+'\' '
                                         'is in wrong type.')
                setattr(self, key, val)
            else:  # if key is not essential attr for model
                self.logger.warning("redundant keyword - [ %s ]", str(key))

        #set elementary parse regex(compiled)
        self.regex_dict = {}

        states_regex = re.compile(r'([^\<\>]*)(?:\<?\-\>)' +
                                  r'(?:([^\<\>]*)(?:\<?\-\>))?([^\<\>]*)')
        self.regex_dict['IS_TS_FS'] = [states_regex, ['IS', 'TS', 'FS']]

        species_regex = re.compile(r'(\d*)([^\_\+\*\<\>]+)_(\d*)(\w+)')
        self.regex_dict['species'] = \
            [species_regex, ['stoichiometry', 'name', 'site_number', 'site']]

        site_regex = re.compile(r'(\d*)(?:\*\_)(\w+)')
        self.regex_dict['empty_site'] = \
            [site_regex, ['stoichiometry', 'site']]

#        self.regex_dict['species_separator'] = \
#        			[r'(?:\A|\s*\+\s*|\s+)', []]

        self.has_absolute_energy = False
        self.has_relative_energy = False

        # load setup file
        if hasattr(self, 'setup_file'):
            self.logger.info('setup file [ %s ] is found', self.setup_file)
            model_name = self.setup_file.rsplit('.', 1)[0]
            setattr(self, 'model_name', model_name)
            self.load(self.setup_file)
            self.logger.info('kinetic modeling...success!\n')
        else:
            self.logger.warning('setup file not read...')

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
        self.logger.info('--- Solve Micro-kinetic model ---')

        # parse data
        self.logger.info('reading data...')
        self.parser.chk_data_validity()
        if relative:
            self.logger.info('use relative energy directly...')
        else:
            self.logger.info('convert relative to absolute energy...')
        self.parser.parse_data(relative=relative)

        # -- solve steady state coverages --
        self.logger.info('passing data to solver...')
        self.solver.get_data()
        self.logger.info('getting rate constants...')
        self.solver.get_rate_constants()
        self.logger.info('getting rate expressions...')
        self.solver.get_rate_expressions(self.solver.rxns_list)

        # energy correction
        if correct_energy:
            self.logger.info('free energy correction...')
            self.solver.correct_energies()

        # solve ODE
        # !! do ODE integration AFTER passing data to solver !!
        if solve_ode:
            self.logger.info("initial coverages = %s", str(init_cvgs))
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
            self.logger.info('use user-defined coverages as initial guess...')

        elif os.path.exists("./data.pkl"):
            with open('data.pkl', 'rb') as f:
                data = cpkl.load(f)
            init_guess = 'steady_state_coverage'
            if init_guess in data:
                self.logger.info('use coverages in data.pkl as initial guess...')
                init_cvgs = data[init_guess]
                coarse_guess = False
            else:
                self.logger.info('use Boltzmann coverages as initial guess...')
                init_cvgs = self.solver.boltzmann_coverages()

        else:  # use Boltzmann coverage
            self.logger.info('use Boltzmann coverages as initial guess...')
            init_cvgs = self.solver.boltzmann_coverages()

        # solve steady state coverages
        # use scipy.optimize.fsolve or not (fast but low-precision)
        if fsolve:
            self.logger.info('using fsolve to get steady state coverages...')
            ss_cvgs = self.solver.fsolve_steady_state_cvgs(init_cvgs)
        else:
            if coarse_guess:
                self.logger.info('getting coarse steady state coverages...')
                init_cvgs = self.solver.coarse_steady_state_cvgs(init_cvgs)  # coarse root
            self.logger.info('getting precise steady state coverages...')
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

    def set_parser(self, parser_name):
        """
        Import parser and set the instance of it as attr of model
        """
        # hacked from CatMap (catmap/model.py)
        # https://docs.python.org/2/library/functions.html#__import__
        basepath = os.path.dirname(
            inspect.getfile(inspect.currentframe()))
        if basepath not in sys.path:
            sys.path.append(basepath)
        # from  loggers import logger
        _module = __import__('parsers', globals(), locals())
        parser_instance = getattr(_module, parser_name)(owner=self)
        setattr(self, 'parser', parser_instance)
        self.logger.info('parser is set.')

    def set_logger(self):
        """
        Get logging.logger instance as logger of kinetic model.
        """
        logger = logging.getLogger('model')
        if os.path.exists('./logging.conf'):
            logging.config.fileConfig('./logging.conf')
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

        self.logger = logger

    def load(self, setup_file):
        """
        Load 'setup_file' into kinetic model by exec setup file
        and assigning all local variables as attrs of model.
        For tools, create the instances of tool classes and
        assign them as the attrs of model.
        """
        self.logger.info('Loading Kinetic Model...\n')

        defaults = dict(
            data_file='data.pkl',
            grid_type='square',
            decimal_precision=100,
            tools=['parser', 'plotter'],
            parser='RelativeEnergyParser',
            table_maker='CsvMaker',
            solver='SteadyStateSolver',
            corrector='ThermodynamicCorrector',
            plotter='EnergyProfilePlotter',
        )

        self.logger.info('read in parameters...')

        #exec setup file set local variables as attrs of model
        globs = {}
        locs = defaults
        execfile(setup_file, globs, locs)

        # customize model tools
        if 'tools' in locs:
            if 'parser' not in locs['tools']:
                raise ParameterError('[ parser ] must be in tools.')
            self.tools = locs['tools']
            del locs['tools']

            self.logger.info('tools = %s', str(self.tools))

        # assign parser ahead to provide essential attrs for other tools
        self.logger.info('instantiate %s', str(locs['parser']))
        self.set_parser(locs['parser'])

        # if parser is kmc_parser use kmc_solver correspondingly
        if locs['parser'] == 'KMCParser':
            locs['solver'] = 'KMCLibSolver'
            self.logger.info('set solver [ KMCSolver ].')

        # assign other tools
        for key in locs.keys():
            # ignore tools which will be loaded later
            if key in self.tools:
                continue
            # check type of variables
            if key in self._attr_type_dict:
                #chech attr type
                if type(locs[key]) != self._attr_type_dict[key]:
                    try:
                        locs[key] = self._attr_type_dict[key](locs[key])
                        setattr(self, key, locs[key])
                    except:
                        raise ValueError('\''+key+'\' is in wrong type. ' +
                                         str(self._attr_type_dict[key]) +
                                         ' object is expected.')
            setattr(self, key, locs[key])
            self.logger.info('%s = %s', key, str(locs[key]))

        #use parser parse essential attrs for other tools
        #parse elementary rxns
        self.logger.info('Parsing elementary rxns...')
        if self.rxn_expressions:
            self.parser.parse_elementary_rxns(self.rxn_expressions)

        # load tools of model
        self.logger.info('instantiate model tools...')
        for key in self.tools:
            # auto-import classes (HACKED from CatMap)
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
                    setattr(self, key, tool_instance)
                    self.logger.info('%s = %s', key, locs[key])
                else:
                    setattr(self, key, None)
                    self.logger.warning('%s is set to None.')
            except ImportError:
                raise ToolsImportError(key.capitalize()+' '+locs[key] +
                                       ' could not be imported. ' +
                                       'Ensure that the class ' +
                                       'exists and is spelled properly.')
        #HACK END

        # kMC parameters check
        solver_type = repr(self.solver).split('.')[2]
        # pseudo random generator
        if solver_type == 'kmc_solver' and 'random_generator' not in locs:
            self.logger.info('pseudo random generator type was not set.')
            self.logger.info('use Mersenne-Twister by default.')

        # kMC control parameters
        if 'kmc_continue' in locs and locs['kmc_continue']:
            # read in step and time info
            self.logger.info('reading in kMC control parameters...')

            # read from auto_coverages.py
            self.logger.info('reading auto_coverages.py...')
            has_auto_coverages = os.path.exists('./auto_coverages.py')
            if has_auto_coverages:
                #raise FilesError('No auto_coverages.py in current path.')
                cvg_globs, cvg_locs = {}, {}
                execfile('./auto_coverages.py', cvg_globs, cvg_locs)
                cvg_step, cvg_time = cvg_locs['steps'][-1], cvg_locs['times'][-1]
                self.logger.info('auto_coverages.py was read.')
            else:
                self.logger.info('auto_coverages.py not read.')

            # read from auto_TOFs.py
            self.logger.info('reading auto_TOFs.py...')
            if not os.path.exists('./auto_TOFs.py'):
                raise FilesError('No auto_TOFs.py in current path.')
            tof_globs, tof_locs = {}, {}
            execfile('./auto_TOFs.py', tof_globs, tof_locs)
            tof_step, tof_time = tof_locs['steps'][-1], tof_locs['times'][-1]
            # get total_rates calculated from last loop
            total_rates = tof_locs['total_rates']
            # get time when the first TOF analysis start
            tof_start_time = tof_locs['start_time']
            self.logger.info('auto_TOFs.py was read.')

            # read from auto_last_types.py
            self.logger.info('reading auto_last_types.py...')
            if not os.path.exists('auto_last_types.py'):
                raise FilesError('No auto_last_types.py in current path.')
            types_globs, types_locs = {}, {}
            execfile('auto_last_types.py', types_globs, types_locs)
            last_types = types_locs['types']
            self.logger.info('auto_last_types.py was read.')

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
            self.start_step, self.start_time = tof_step, tof_time
            self.tof_start_time = tof_start_time
            self.start_types = last_types
            self.total_rates = total_rates
            self.logger.info('kMC analysis starts from step = %d, time = %e s',
                             self.start_step, self.start_time)

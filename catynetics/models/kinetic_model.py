import copy
import logging
import os
from operator import add


from ..compatutil import reduce
from ..descriptors.descriptors import *
from ..descriptors.component_descriptors import Component
from ..mpicommons import mpi
from ..database.thermo_data import kB_eV, h_eV
from ..errors.error import *
from ..functions import *
from ..utilities.profiling_utitlities import do_cprofile


class KineticModel(object):
    ''' The base class for kinetic models.

    :param setup_file: kinetic model set up file
    :type setup_file: str

    :param setup_dict: A dictionary contains essential setup parameters for kinetic model
    :type setup_dict: dict
    
    :param logger_level: logging level for logger
    :type logger_level: int

    :param file_handler_level: logging level for file handler
    :type file_handler_level: int

    :param console_handler_level: logging level for console handler
    :type console_handler_level: int

    Example::

        >>> from catynetics.models import KineticModel
        >>> model = KineticModel(setup_file="setup.mkm", logger_level=30)

    Attributes:
        setup_file(:obj:`str`): kinetic model setup file

        setup_dict(:obj:`dict`): A dictionary contains essnetial setup parameters for kinetic model

        logger_level(:obj:`int`): The logging level for built-in logger

        file_handler_level(:obj:`int`): The logging level for file handler

        console_handler_level(:obj:`int`): The logging level for console handler

        parser(:obj:`object`): Kinetic Model parser,
        possible type could be :obj:`catynetics.parsers.RelativeEnergyParser`,
        :obj:`catynetics.parsers.AbsoluteEnergyParser`, :obj:`catynetics.parsers.KMCParser`

        solver(:obj:`object`): Kinetic Model solver,
        possible type could be :obj:`catynetics.solvers.SteadyStateSolver`,
        :obj:`catynetics.solvers.KMCSolver`

        corrector(:obj:`object`): Kinetic Model data corrector,
        possible type could be :obj:`ThermodynamicCorrector`

        temperature(:obj:`float`): Environment temperature in K, default is 298.0K

        rxn_expressions(:obj:`list` of :obj:`str`): Reaction expressions for elementary reactions,
        default is []

        species_definitions(:obj:`dict`): All species definition in kinetic model after parsing

        rate_algo(:obj:`str`): Algorithm for rate calculation, could be 'TST' for
        transition state theory or 'CT' for collision theory, default value is 'TST'

    '''

    # Attribute descriptors.
    # {{{
    setup_file = String("setup_file", default="")
    setup_dict = Dict("setup_dict", default={})

    logger_level = Integer("logger_level",
                              default=logging.INFO,
                              candidates=range(0, 60, 10))

    # Logging level for file handler.
    file_handler_level = Integer("file_handler_level",
                                    default=logging.DEBUG,
                                    candidates=range(0, 60, 10))

    # Logging level for console handler.
    console_handler_level = Integer("console_handler_level",
                                       default=logging.INFO,
                                       candidates=range(0, 60, 10))

    parser = Component("parser", default=None, candidates=["RelativeEnergyParser",
                                                           "AbsoluteEnergyParser",
                                                           "KMCParser"])

    solver = Component("solver", default=None, candidates=["KMCSolver", "SteadyStateSolver"])

    corrector = Component("corrector",
                          default=None,
                          candidates=["ThermodynamicCorrector"])

    plotter = Component("plotter",
                        default=None,
                        candidates=["EnergyProfilePlotter"])

    # Temperature.
    temperature = Float("temperature", default=298.0)

    # Reaction expressions.
    rxn_expressions = Sequence("rxn_expressions", default=[], entry_type=str)

    # Definition dict of species.
    species_definitions = SpeciesDefinitions("species_definitions", default={})

    # Algorithm for rate calculation.
    rate_algo = String("rate_algo", default="TST")

    # Model core components.
    components = Sequence("components", default=["parser"], entry_type=str)

    # Area of unit cell (m^2).
    unitcell_area = Float("unitcell_area", default=0.0)

    # Ratio of active area.
    active_ratio = Float("active_ratio", default=1.0)

    # }}}

    def __init__(self, **kwargs):
        # {{{

        # Get all kwargs.
        setup_file = kwargs.pop("setup_file", None)
        setup_dict = kwargs.pop("setup_dict", None)

        self.logger_level = kwargs.pop("logger_level", logging.INFO)
        self.file_handler_level = kwargs.pop("file_handler_level", logging.DEBUG)
        self.console_handler_level = kwargs.pop("console_handler_level", logging.INFO)

        # Physical constants.
        self._kB = kB_eV
        self._h = h_eV

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
            exec(open(self.setup_file, "rb").read(), globs, locs)
            self.setup_dict = locs

        # Set logger.
        self._set_logger()

        # Output MPI info.
        if self.log_allowed:
            self._logger.info("------------------------------------")
            self._logger.info(" Model is runing in MPI Environment ")
            self._logger.info(" Number of process: {}".format(mpi.size))
            self._logger.info("------------------------------------")
            self._logger.info(" ")

        # Energy flags.
        self._has_absolute_energy = False
        self._has_relative_energy = False
        self._relative_energies = {}

        # Load setup file.
        self._load(self.setup_dict)
        if self.log_allowed:
            self._logger.info('{} creating...success!\n'.format(self.__class__.__name__))
        # }}}

    def _set_logger(self, filename=None):
        """ Private function to get logging.logger instance as logger of kinetic model.
        """
        # {{{
        # Create logger.
        logger = logging.getLogger('model')
        logger.setLevel(self.logger_level)

        # Set log file name.
        if filename is None:
            if 1 == mpi.size:
                filename = "out.log"
            else:
                mpi.barrier()
                if not os.path.exists("./log") and mpi.is_master:
                    os.mkdir("./log")
                filename = "./log/out_{}.log".format(mpi.rank)

        # Create handlers.
        std_hdlr = logging.FileHandler(filename)
        std_hdlr.setLevel(self.file_handler_level)
        console_hdlr = logging.StreamHandler()
        console_hdlr.setLevel(self.console_handler_level)

        # Create formatter and add it to the handlers.
        formatter = logging.Formatter('%(name)s   %(levelname)-8s %(message)s')
        std_hdlr.setFormatter(formatter)
        console_hdlr.setFormatter(formatter)

        # Add the handlers to logger.
        logger.addHandler(std_hdlr)
        logger.addHandler(console_hdlr)

        self._logger = logger
        # }}}

    def set_logger_level(self, handler_type, level):
        """ Set the logging level of logger handler.

        :param handler_type: logger handler name
        :type handler_type: str

        :param level: logging level
        :type level: int

        :return: Old logging level of the handler
        :rtype: int
        """
        # Locate handler.
        handler = None
        for h in self._logger.handlers:
            if h.__class__.__name__ == handler_type:
                handler = h
                break

        if handler is None:
            raise ValueError("Unknown handler type '{}'".format(handler_type))

        # Reset logging level.
        old_level = handler.level
        handler.setLevel(level)

        return old_level

    def clear_handlers(self):
        """ Clear all handlers in logger.
        """
        self._logger.handlers = []

    def __mro_class_attrs(self):
        """ Private helper function to get all class attribute names(include father classes) .
        """
        d = {}

        for cls in self.__class__.__mro__:
            d.update(cls.__dict__)

        return d.keys()

    def _load(self, setup_dict):
        """ Load 'setup_file' into kinetic model by exec setup file
        and assigning all local variables as attrs of model.
        For tools, create the instances of tool classes and
        assign them as the attrs of model.
        """
        # {{{
        if self.log_allowed:
            self._logger.info('Loading Kinetic Model...\n')
            self._logger.info('read in parameters...')

        # Set model attributes in setup file.
        class_attrs = self.__mro_class_attrs()
        for key, value in setup_dict.items():
            # Check redundant parameter.
            if key not in class_attrs:
                if self.log_allowed:
                    msg = "Found redundant parameter '{}'".format(key)
                    self._logger.warning(msg)
                continue

            # Parser & solver will be set later.
            if key in ["parser", "solver"]:
                continue

            # Set parameters in setup dict as attiributes of model.
            setattr(self, key, value)

            # Output info.
            specials = ("rxn_expressions", "species_definitions")
            if key not in specials:
                if self.log_allowed:
                    self._logger.info('{} = {}'.format(key, str(value)))

            # If it is a iterable, loop to output.
            else:
                if self.log_allowed:
                    self._logger.info("{} =".format(key))
                if type(setup_dict[key]) is dict:
                    for k, v in value.items():
                        if self.log_allowed:
                            self._logger.info("        {}: {}".format(k, v))
                else:
                    for item in value:
                        if self.log_allowed:
                            self._logger.info("        {}".format(item))

        # Instantialize parser.
        if "parser" not in setup_dict:
            raise ParameterError("Parser is not set in kinetic model.")
        self.parser = setup_dict["parser"]

        # use parser parse essential attrs for other tools
        # Parse elementary rxns
        if self.log_allowed:
            self._logger.info('Parsing elementary rxns...')
        if self.rxn_expressions:
            (self.__adsorbate_names,
             self.__gas_names,
             self.__liquid_names,
             self.__site_names,
             self.__transition_state_names,
             self.__elementary_rxns_list) = \
                self.parser.parse_elementary_rxns(self.rxn_expressions)

        # Instantialize solver.
        if "solver" in setup_dict:
            self.solver = setup_dict["solver"]
        # }}}

    def generate_relative_energies_file(self, filename="rel_energy.py"):
        """ Generate a energy input file containing relative energies for 
        all elementary reactions.

        :param filename: The name of relative input file, default value is 'rel_energy.py'.
        :type filename: str

        """
        content = ("# Relative Energies for all elementary reactions.\n" +
                   "Ga, dG = [], []\n\n")

        for rxn_expression in self.rxn_expressions:
            rxn_content = "# {}\nGa.append()\ndG.append()\n\n".format(rxn_expression)
            content += rxn_content

        with open(filename, "w") as f:
            f.write(content)

    def generate_absolute_energies_file(self, filename="abs_energy.py"):
        """ Generate a energy input file containing absolute energies
        for all species(including sites).

        :param filename: The name of absolute energy input file, default value is 'abs_energy.py'.
        :type filename: str
        """
        content = ("# Absolute energies for all species.\n" +
                   "absolute_energies = {\n\n")

        all_species = reduce(add, [self.gas_names,
                                   self.liquid_names,
                                   self.adsorbate_names,
                                   self.transition_state_names,
                                   self.site_names])

        for sp in all_species:
            content += "    '{}': 0.0, # eV\n\n".format(sp)

        content += "}\n\n"

        with open(filename, "w") as f:
            f.write(content)

    def run(self, *kwargs):
        """ Need implementation.
        """
        raise NotImplementedError

    def set_solver(self, solver):
        ''' Set solver for current model

        :param solver: Valid solver object
        :type solver: any
        '''
        setattr(self, mangled_name(self, 'solver'), solver)

    def set_plotter(self, plotter):
        ''' Set plotter for current model

        :param plotter: Valid plotter object
        :type plotter: any
        '''
        setattr(self, mangled_name(self, 'plotter'), solver)

    def set_corrector(self, corrector):
        ''' Set corrector for current model

        :param corrector: Valid corrector object
        :type corrector: any
        '''
        setattr(self, mangled_name(self, 'corrector'), corrector)

    @Property
    def log_allowed(self):
        """
        Flag for if log output is allowed.
        """
        # This function should be overridden in sub-class.
        return False

    @Property
    def kB(self):
        return self._kB

    @Property
    def h(self):
        return self._h

    @Property
    def logger(self):
        """
        Query function for model logger.
        """
        return self._logger

    @Property
    def elementary_rxns_list(self):
        """
        Query function for elementary reactions list.
        """
        return self.__elementary_rxns_list

    @Property
    def site_names(self):
        """
        Query function for site names in model.
        """
        return self.__site_names

    @Property
    def adsorbate_names(self):
        """
        Query function for adsorbate names in model.
        """
        return self.__adsorbate_names

    @Property
    def gas_names(self):
        """
        Query function for gas names in model.
        """
        return self.__gas_names

    @Property
    def liquid_names(self):
        """
        Query function for liquid names in model.
        """
        return self.__liquid_names

    @Property
    def transition_state_names(self):
        """
        Query function for transition state species names in model.
        """
        return self.__transition_state_names

    @Property
    def has_relative_energy(self):
        """
        Query function for relative energy flag.
        """
        return self._has_relative_energy

    @Property
    def has_absolute_energy(self):
        """
        Query function for absolute energy flag.
        """
        return self._has_absolute_energy

    @Property
    def relative_energies(self):
        """
        Query function for relative energy in data file.
        """
        return self._relative_energies

    @Property
    def absolute_energies(self):
        """
        Query function for absolute energy in data file.
        """
        return self._absolute_energies


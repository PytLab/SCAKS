import cPickle as cpkl
import copy
import logging
import logging.config
import os
import sys

import kynetix.descriptors.descriptors as dc
import kynetix.descriptors.component_descriptors as cpdc
from kynetix import mpi_master, mpi_size, mpi_installed, mpi_rank
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
    species_definitions = dc.Dict("species_definitions", default={}, deepcopy=True)

    # Model core components.
    components = dc.Sequence("components", default=["parser"], entry_type=str)

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
            execfile(self.setup_file, globs, locs)
            self.setup_dict = locs

        self.verbosity = verbosity

        # Set logger.
        self._set_logger()

        # Output MPI info.
        if self.log_allowed:
            self._logger.info("------------------------------------")
            self._logger.info(" Model is runing in MPI Environment ")
            self._logger.info(" Number of process: {}".format(mpi_size))
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
        """
        Private function to get logging.logger instance as logger of kinetic model.
        """
        # {{{
        # Create logger.
        logger = logging.getLogger('model')
        logger.setLevel(self.verbosity)

        # Set log file name.
        if filename is None:
            if not os.path.exists("./log"):
                os.mkdir("./log")
            if not mpi_installed or 1 == mpi_size:
                filename = "./log/out.log"
            else:
                filename = "./log/out_{}.log".format(mpi_rank)

        # Create handlers.
        std_hdlr = logging.FileHandler(filename)
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

        self._logger = logger
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

    def __mro_class_attrs(self):
        """
        Private helper function to get all class attribute names(include father classes) .
        """
        d = {}

        for cls in self.__class__.__mro__:
            d.update(cls.__dict__)

        return d.keys()

    def _load(self, setup_dict):
        """
        Load 'setup_file' into kinetic model by exec setup file
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
        for key, value in setup_dict.iteritems():
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
                    for k, v in value.iteritems():
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

    def run(self, *kwargs):
        pass

    @dc.Property
    def log_allowed(self):
        """
        Flag for if log output is allowed.
        """
        # This function should be overridden in sub-class.
        return False

    @dc.Property
    def kB(self):
        return self._kB

    @dc.Property
    def h(self):
        return self._h

    @dc.Property
    def logger(self):
        """
        Query function for model logger.
        """
        return self._logger

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
        return self._has_relative_energy

    @dc.Property
    def has_absolute_energy(self):
        """
        Query function for absolute energy flag.
        """
        return self._has_absolute_energy

    @dc.Property
    @return_deepcopy
    def relative_energies(self):
        """
        Query function for relative energy in data file.
        """
        return self._relative_energies


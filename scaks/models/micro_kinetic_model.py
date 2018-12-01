'''
Module for micro-kinetic model class definition.
'''

import logging
import os
from functools import wraps

from .kinetic_model import KineticModel
from ..errors.error import ParameterError
from ..mpicommons import mpi
from ..descriptors.descriptors import *
from ..utilities.profiling_utitlities import do_cprofile
from ..compatutil import pickle
from ..plugins.analysis import OnTheFlyAnalysis


class MicroKineticModel(KineticModel):
    ''' Class for micro-kinetic model.
    
    :param setup_file: kinetic model set up file
    :type setup_file: int


    :param setup_dict: A dictionary contains essential setup parameters for kinetic model.
    :type setup_dict: dict

    :param logger_level: logging level
    :type logger_level: int

    :param file_handler_level: logging level for file handler
    :type file_handler_level: int

    :param console_handler_level: logging level for console handler
    :type console_handler_level: int

    :returns: :obj:`None`

    Example::

        >>> from scaks.models.kinetic_model import MicroKineticModel
        >>> model = MicroKineticModel(setup_file="setup.mkm", logger_level=30)

    Attributes:
        decimal_precision (:obj:`int`): Data precision, default value is 100

        perturbation_size (:obj:`float`): Perturbation size for numerical jacobian matrix,
        default value is 0.01

        perturbation_direction (:obj:`str`): Direction of perturbation, :obj:`left` or :obj:`right`,
        default is :obj:`'right'`

        archived_variables (:obj:`str`): Archived variables. Possible values:
        'initial_guess', 'steady_state_coverages', 'steady_state_error', 'rates',
        'net_rates', 'reversibilities', 'tofs'

        numerical_representation (:obj:`str`): Numerical representation method,
        value could be 'mpmath' or 'sympy'.

        rootfinding(:obj:`str`): Rootfinding iterator type, default value is 'MDNewton',
        possible value can be 'MDNewton' or 'ConstrainedNewton'

        tolerance (:obj:`float`): Iteration tolerance, default is 1e-8

        max_rootfinding_iterations (:obj:`int`): Max iteraction steps, default is 100

        ode_buffer_size (:obj:`int`): Ode integration buffer size, default is 500

        ode_output_interval (:obj:`int`): Ode ouptut interval, default is 200

        data_file (:obj:`str`): Filename for data store, default is 'data.pkl'

        log_allowed (:obj:`bool`): If log output is allowed, default is `True`

        TOFs (:obj:`list` of :obj:`float`): Turnover frequencies calculated.

        reversibilities (:obj:`list` of :obj:`float`): Reversibilities calculated.

        error (:obj:`list` of :obj:`float`): Residual errors after iterations.

        analysis (:obj:`list` of :obj:`scaks.plugins.analysis.OnTheFlyAnalysis`: On-the-fly analysis plugins
    '''

    # {{{
    decimal_precision = Integer("decimal_precision", default=100)

    perturbation_size = Float("perturbation_size", default=0.01)

    perturbation_direction = String("perturbation_direction",
                                    default="right",
                                    candidates=["right", "left"])

    archived_variables = Sequence("archive_data",
                                  default=["steady_state_coverages"],
                                  entry_type=str)

    numerical_representation = String("numerical_representation",
                                      default="mpmath",
                                      candidates=["mpmath", "sympy"])

    rootfinding = String("rootfinding",
                         default="MDNewton",
                         candidates=["MDNewton", "ConstrainedNewton"])

    tolerance = Float("tolerance", default=1e-8)

    max_rootfinding_iterations = Integer("max_rootfinding_iterations",
                                         default=100)

    ode_buffer_size = Integer("ode_buffer_size", default=500)

    ode_output_interval = Integer("ode_output_interval", default=200)

    # Reference energies used to calculate formation energy.
    ref_energies = RefEnergies("ref_energies", default={})

    # On-the-fly analysis plugins
    analysis = Sequence('analysis', default=[], entry_type=OnTheFlyAnalysis)
    # }}}

    def __init__(self, **kwargs):
        super(MicroKineticModel, self).__init__(**kwargs)

        # Create data directory if need.
        if mpi.size != 1 and not os.path.exists("./data"):
            mpi.barrier()
            if mpi.is_master:
                os.mkdir("./data")

        # Model attributes definitions.
        self.__ss_cvgs = None          # steady-state coverages
        self.__tofs = None              # turn-over frequencies
        self.__reversibilities = None  # reversibilities

    def _set_logger(self, filename=None):
        super(MicroKineticModel, self)._set_logger(filename)
        # if not master processor, no INFO to console.
        if not mpi.is_master:
            self.set_logger_level("StreamHandler", logging.WARNING)

    def run(self, **kwargs):
        """
        Run function to solve Micro-kinetic model using Steady State Approxmiation
        to get steady state coverages and turnover frequencies (TOF).

        :param init_cvgs: Initial guess for coverages
        :type init_cvgs: :obj:`list` of :obj:`float`

        :param relative_energies: Relative energies for all elementary reactions.
        :type relative_energies: dict

        .. note::
            keys :obj:`Gaf` and :obj:`Gar` must be in relative energies dict.
            e.g. :obj:`{"Gaf": [...], "Gar": [...], "dG": [...]}`

        :param fsolve: use scipy.optimize.fsolve to get low-precision root or not
        :type fsolve: bool

        :param coarse_guess: use fsolve to do initial coverages preprocessing or not
        :type coarse_guess: bool

        :param XRC: calculate degree of rate control or nor
        :type XRC: bool

        :param epsilon: the change of energy for XRC calculation, default is 10-5.
        :type epsilon: float

        :param product_name: Production name of the model. e.g. :obj:`'CH3OH_g'`
        :type product_name: str

        :returns: :obj:`None`
        """
        # {{{
        # Setup default parameters.
        init_cvgs = kwargs.pop("init_cvgs", None)
        relative_energies = kwargs.pop("relative_energies", None)
        fsolve = kwargs.pop("fsolve", False)
        coarse_guess = kwargs.pop("coarse_guess", True)
        XRC = kwargs.pop("XRC", False)
        epsilon = kwargs.pop("epsilon", 1e-5)
        product_name = kwargs.pop("product_name", None)

        if kwargs:
            for key in kwargs:
                msg = "Found redundant keyword argument: {}".format(key)
                self._logger.warning(msg)

        if self.log_allowed:
            self._logger.info('--- Solve Micro-kinetic model ---')

        solver = self.solver

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
            if len(init_cvgs) != len(self.adsorbate_names):
                msg = "init_cvgs must have {} elements, but {} is supplied"
                msg = msg.format(len(self.adsorbate_names), len(init_cvgs))
                raise ParameterError(msg)

            if self.log_allowed:
                self._logger.info('use user-defined coverages as initial guess...')

        elif os.path.exists(self.data_file):
            with open(self.data_file, 'rb') as f:
                data = pickle.load(f)
            init_guess = 'steady_state_coverage'
            if init_guess in data:
                if self.log_allowed:
                    msg = 'use coverages in {} as initial guess...'.format(self.data_file)
                    self._logger.info(msg)
                init_cvgs = data[init_guess]
                coarse_guess = False
            else:
                if self.log_allowed:
                    self._logger.info('Do ODE integration to get initial guess...')
                ode_traj = solver.solve_ode()
                init_cvgs = ode_traj[-1]

        else:
            if self.log_allowed:
                self._logger.info('Do ODE integration to get initial guess...')
            ode_traj = solver.solve_ode()
            init_cvgs = ode_traj[-1]

        # Solve steady state coverages.
        # Use scipy.optimize.fsolve or not (fast but low-precision).
        if fsolve:
            if self.log_allowed:
                self._logger.info('using fsolve to get steady state coverages...')
            self.__ss_cvgs = solver.fsolve_steady_state_cvgs(c0=init_cvgs,
                                                             relative_energies=relative_energies)
        else:
            if coarse_guess:
                if self.log_allowed:
                    self._logger.info('getting coarse steady state coverages...')
                init_cvgs = solver.coarse_steady_state_cvgs(c0=init_cvgs,
                                                            relative_energies=relative_energies)
            if self.log_allowed:
                self._logger.info('getting precise steady state coverages...')
            self.__ss_cvgs = solver.get_steady_state_cvgs(c0=init_cvgs,
                                                          relative_energies=relative_energies)

        # Output rate constants for all elementary reactions.
        solver.get_rate_constants(relative_energies=relative_energies, log=True)

        # Get steady state rates for all elementary reactions.
        rf, rr = solver.get_rates(cvgs_tuple=self.__ss_cvgs,
                                  relative_energies=relative_energies,
                                  log=True)

        # Get TOFs for gases.
        self.__tofs = solver.get_tof(cvgs=self.__ss_cvgs,
                                     relative_energies=relative_energies)

        # Get reversibilities.
        self.__reversibilities = solver.get_reversibilities(rf, rr)

        # Get residual error.
        self.__error = solver.error

        # Calculate XRC.
        if XRC:
            if product_name is None:
                raise ParameterError("production name must be provided to get XRC.")

            solver.get_single_XRC(product_name,
                                  epsilon=epsilon,
                                  relative_energies=relative_energies)
        # }}}

    def hybrid_method_register(self, fn):
        ''' A decorator for hybrid method function register to current model

        :param fn: A function cooperate with Newton's method to provide a hybrid
                   iteration method.
        :type fn: function
        '''
        @wraps(fn)
        def _hybrid_method(model):
            if not isinstance(model, MicroKineticModel):
                raise TypeError('model must be a MicroKineticModel object')
            return fn(model)

        self.hybrid_method = _hybrid_method

    def analysis_register(self, analysis_cls):
        ''' A decorator for analysis regsiter.

        :param analysis_cls: The analysis to be registered
        :type analysis_cls: :obj:`scaks.plugin_interfaces.OnTheFlyAnalysis`
        '''
        if not issubclass(analysis_cls, OnTheFlyAnalysis):
            raise TypeError('analysis class must be subclass of OnTheFlyAnalysis')
        
        # Add analysis instance to model
        analysis = analysis_cls()
        self.analysis.append(analysis)

    @Property
    def model_info(self):
        """
        Generate a report dict containing model information.
        """
        info = {}
        info['gas_names'] = self.gas_names
        info['adsorbate_names'] = self.adsorbate_names
        info['steady_state_coverages'] = \
            [float(cvg) for cvg in self.steady_state_coverages]
        info['TOFs'] = [float(tof) for tof in self.TOFs]
        info['reversibilities'] = self.reversibilities
        info['rxn_expressions'] = self.rxn_expressions

        return info

    @Property
    def data_file(self):
        ''' Get the name of file where serialzed data stored.
        '''
        if mpi.size == 1:
            return "data.pkl"
        else:
            return "./data/data_{}.pkl".format(mpi.rank)

    @Property
    def log_allowed(self):
        """
        Flag for if log output is allowed.
        """
        # All processors can output log information.
        return True

    @Property
    def steady_state_coverages(self):
        try:
            return self.__ss_cvgs
        except AttributeError:
            raise AttributeError("Unsolved model has no steady state coverages.")

    @Property
    def TOFs(self):
        try:
            return self.__tofs
        except AttributeError:
            raise AttributeError("Unsolved model has no turnover frequencies.")

    @Property
    def reversibilities(self):
        try:
            return self.__reversibilities
        except AttributeError:
            raise AttributeError("Unsolved model has no reversibilities.")

    @Property
    def error(self):
        try:
            return self.__error
        except AttributeError:
            raise AttributeError("Unsolved model has no error.")


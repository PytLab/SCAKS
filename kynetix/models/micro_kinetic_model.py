import cPickle as cpkl
import logging
import os

from kynetix.mpicommons import mpi
import kynetix.models.kinetic_model as km
import kynetix.descriptors.descriptors as dc
import kynetix.descriptors.component_descriptors as cpdc
from kynetix.utilities.profiling_utitlities import do_cprofile


class MicroKineticModel(km.KineticModel):

    # {{{
    # Data precision.
    decimal_precision = dc.Integer("decimal_precision", default=100)

    # Perturbation size for numerical jacobian matrix.
    perturbation_size = dc.Float("perturbation_size", default=0.01)

    # Direction of perturbation.
    perturbation_direction = dc.String("perturbation_direction", 
                                       default="right",
                                       candidates=["right", "left"])

    # Archived variables.
    # Candidates: 'initial_guess', 'steady_state_coverages', 'steady_state_error',
    #             'rates', 'net_rates', 'reversibilities', 'tofs'
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

    # Reference energies used to calculate formation energy.
    ref_energies = dc.RefEnergies("ref_energies", default={})
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
        >>> from kynetix.models.kinetic_model import MicroKineticModel
        >>> model = MicroKineticModel(setup_file="setup.mkm",
                                      logger_level=logging.WARNING)
        """
        super(MicroKineticModel, self).__init__(**kwargs)

        # Create data directory if need.
        if mpi.size != 1 and not os.path.exists("./data"):
            mpi.barrier()
            if mpi.is_master:
                os.mkdir("./data")

    def _set_logger(self, filename=None):
        super(MicroKineticModel, self)._set_logger(filename)
        # if not master processor, no INFO to console.
        if not mpi.is_master:
            self.set_logger_level("StreamHandler", logging.WARNING)

    @do_cprofile("./mkm_run.prof")
    def run(self, **kwargs):
        """
        Function to solve Micro-kinetic model using Steady State Approxmiation
        to get steady state coverages and turnover frequencies.

        Parameters:
        -----------
        init_cvgs: Initial guess for coverages, tuple of floats.

        relative_energies: Relative energies for all elementary reactions, dict.
            NOTE: keys "Gaf" and "Gar" must be in relative energies dict.
            e.g. {"Gaf": [...], "Gar": [...], "dG": [...]}

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
        relative_energies = kwargs.pop("relative_energies", None)
        solve_ode = kwargs.pop("solve_ode", False)
        fsolve = kwargs.pop("fsolve", False)
        coarse_guess = kwargs.pop("coarse_guess", True)
        XRC = kwargs.pop("XRC", False)
        product_name = kwargs.pop("product_name", None)
        data_file = kwargs.pop("data_file", "./rel_energy.py")

        if kwargs:
            for key in kwargs:
                msg = "Found redundant keyword argument: {}".format(key)
                self._logger.warning(msg)

        if self.log_allowed:
            self._logger.info('--- Solve Micro-kinetic model ---')

        # Get parser and solver.
        parser = self.__parser
        solver = self.__solver

        # Parse data.
        if self.log_allowed:
            self._logger.info('reading data...')
        parser.parse_data(filename=data_file)

        # -- solve steady state coverages --
        if self.log_allowed:
            self._logger.info('passing data to solver...')
        solver.get_data()

        # solve ODE
        # !! do ODE integration AFTER passing data to solver !!
        if solve_ode:
            if self.log_allowed:
                self._logger.info("initial coverages = %s", str(init_cvgs))
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
            if len(init_cvgs) != len(self.adsorbate_names):
                msg = "init_cvgs must have {} elements, but {} is supplied"
                msg = msg.format(len(self.adsorbate_names), len(init_cvgs))
                raise ParameterError(msg)

            if self.log_allowed:
                self._logger.info('use user-defined coverages as initial guess...')

        elif os.path.exists(self.data_file):
            with open(self.data_file, 'rb') as f:
                data = cpkl.load(f)
            init_guess = 'steady_state_coverage'
            if init_guess in data:
                if self.log_allowed:
                    msg = 'use coverages in {} as initial guess...'.format(self.data_file)
                    self._logger.info(msg)
                init_cvgs = data[init_guess]
                coarse_guess = False
            else:
                if self.log_allowed:
                    self._logger.info('use Boltzmann coverages as initial guess...')
                init_cvgs = solver.boltzmann_coverages()

        else:  # use Boltzmann coverage
            if self.log_allowed:
                self._logger.info('use Boltzmann coverages as initial guess...')
            init_cvgs = solver.boltzmann_coverages()

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
        reversibilities = solver.get_reversibilities(rf, rr)

        # Calculate XRC.
        if XRC:
            if product_name is None:
                raise ParameterError("production name must be provided to get XRC.")

            solver.get_single_XRC(product_name,
                                  epsilon=1e-5,
                                  relative_energies=relative_energies)

        return
        # }}}

    @dc.Property
    def data_file(self):
        if mpi.size == 1:
            return "data.pkl"
        else:
            return "./data/data_{}.pkl".format(mpi.rank)

    @dc.Property
    def log_allowed(self):
        """
        Flag for if log output is allowed.
        """
        # All processors can output log information.
        return True

    @dc.Property
    def steady_state_coverages(self):
        try:
            return self.__ss_cvgs
        except AttributeError:
            raise AttributeError("Unsolved model has no steady state coverages.")

    @dc.Property
    def TOFs(self):
        try:
            return self.__tofs
        except AttributeError:
            raise AttributeError("Unsolved model has no turnover frequencies.")

    @dc.Property
    def reversibilities(self):
        try:
            return self.__reversibilities
        except AttributeError:
            raise AttributeError("Unsolved model has no reversibilities.")


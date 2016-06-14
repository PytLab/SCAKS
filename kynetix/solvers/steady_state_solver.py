import copy
import logging
import random
import re
import signal

from scipy.integrate import odeint, ode
from scipy.linalg import norm
from scipy.optimize import fsolve

from kynetix import file_header
from kynetix.errors.error import *
from kynetix.functions import get_list_string
from kynetix.parsers.rxn_parser import *
from kynetix.solvers.rootfinding_iterators import *
from kynetix.solvers.solver_base import *


class SteadyStateSolver(SolverBase):
    def __init__(self, owner):
        # {{{
        super(SteadyStateSolver, self).__init__(owner)

        # set logger
        self.__logger = logging.getLogger('model.solvers.SteadyStateSolver')

        # Set default parameter dict
        defaults = dict(rootfinding='MDNewton',
                        tolerance=1e-8,
                        max_rootfinding_iterations=100,
                        residual_threshold=1.0,
                        initial_guess_scale_factor=100,
                        stable_criterion=1e-10)
        defaults = self.update_defaults(defaults)

        # Set varibles in defaults protected attributes of solver.
        protected_defaults = {"_{}".format(key): value
                              for key, value in defaults.iteritems()}
        self.__dict__.update(protected_defaults)
        # }}}

    def __constrain_coverages(self, cvgs_tuple):
        """
        Private function to constrain coverages of absorbates
        between 0.0 and 1.0 or total number.
        """
        # {{{
        species_definitions = self._owner.species_definitions()
        adsorbate_names = self._owner.adsorbate_names()
        site_names = self._owner.site_names()

        # Convert tuple to dict
        cvgs_dict = self._cvg_tuple2dict(cvgs_tuple)

        # Enforce explicit maxima, cannot be larger than 1.0, smaller than 0.0
        for adsorbate_name in adsorbate_names:
            if cvgs_dict[adsorbate_name] > 1.0:
                cvgs_dict[adsorbate_name] = self._mpf('1.0')
            if cvgs_dict[adsorbate_name] < 0.0:
                cvgs_dict[adsorbate_name] = self._mpf('0.0')

        # Enforce explicit maxima, cannot be larger than site's total number
        for site_name in site_names:
            total_cvg = species_definitions[site_name]['total']
            for adsorbate_name in self._classified_adsorbates[site_name]:
                if cvgs_dict[adsorbate_name] > total_cvg:
                    cvgs_dict[adsorbate_name] = self._mpf(total_cvg)

        def constrain_in_total(cvgs_dict, max_cvg):
            "Make sure sum of cvgs in cvgs_dict is not larger than max_cvg"
            total_cvg = sum(cvgs_dict.values())
            if total_cvg > max_cvg:
                for key in cvgs_dict:
                    cvgs_dict[key] = cvgs_dict[key]/total_cvg*max_cvg
            return cvgs_dict

        # Sum of cvgs on one type of surface <= site total e.g 1.0
        for site_name in site_names:
            max_cvg = species_definitions[site_name]['total']
            sub_cvgs_dict = {}

            # Add cvgs of the site to sub_cvgs_dict
            for adsorbate_name in self._classified_adsorbates[site_name]:
                sub_cvgs_dict.setdefault(adsorbate_name,
                                         cvgs_dict[adsorbate_name])
            # Add free site coverage
            sub_cvgs_dict = constrain_in_total(sub_cvgs_dict, max_cvg)
            cvgs_dict.update(sub_cvgs_dict)
        # Convert dict to tuple, and return
        constrained_cvgs_tuple = self._cvg_dict2tuple(cvgs_dict)

        def compare_cvgs(cvgs1, cvgs2):
            "Compare two coverage tuples."
            if len(cvgs1) != len(cvgs2):
                self.__logger.warning('coverage length inconsistency is detected.')
                return False
            for cvg1, cvg2 in zip(cvgs1, cvgs2):
                if abs(cvg1 - cvg2) > 10e-20:
                    return False
            return True

        consistant = compare_cvgs(constrained_cvgs_tuple, cvgs_tuple)
        # log if constraint has been carried out
        if not consistant:
            self.__logger.warning('coverage constraining...\n')
            self.__logger.debug('    initial coverage: %s', str(map(float, cvgs_tuple)))
            self.__logger.debug('constrained coverage: %s\n',
                                str(map(float, constrained_cvgs_tuple)))

        return constrained_cvgs_tuple
        # }}}

    def get_elementary_dtheta_dt_expression(self,
                                            adsorbate_name,
                                            rxn_expression):
        """
        Function to get dtheta_dt of the corresponding adsorbate in single
        elementary equation.

        Parameters:
        -----------
        adsorbate_name: The adsorbate name whose coverage is derived wrt time, str.

        rxn_expression: The string expression of an elementary reaction.

        Returns:
        --------
        The dtheta/dt expression, str.

        Example:
        --------
        >>> s.get_elementary_dtheta_dt_expression("O_s", 'O2_g + 2*_s -> 2O_s')
        >>> "2*kf[1]*p['O2_g']*theta['*_s']**2 - 2*kr[1]*theta['O_s']**2"
        """
        # {{{
        # Check adsorbate name.
        if adsorbate_name not in self._owner.adsorbate_names():
            raise ValueError("'{}' is not an adsorbate!".format(adsorbate_name))

        # Get formula object list of the rxn_expression.
        rxn_equation = RxnEquation(rxn_expression)
        elementary_rxn_list = rxn_equation.to_formula_list()

        # Find formula list index.
        for formula_list in elementary_rxn_list:
            for formula in formula_list:
                stoichiometry = formula.stoichiometry()
                current_species = formula.species_site()
                if current_species == adsorbate_name:
                    break

            if current_species == adsorbate_name:
                # Get state idx to get direction info.
                state_idx = elementary_rxn_list.index(formula_list)
                break

        # If adsorbate name not in reaction, stop.
        if current_species != adsorbate_name:
            return

        # Get rate expression of the elementary equation.
        f_expr, r_expr = self.get_elementary_rate_expression(rxn_expression)

        # Adsorbate is consumed
        if state_idx == 0:
            if stoichiometry == 1:
                increase_rate, decrease_rate = r_expr, f_expr
            else:
                increase_rate, decrease_rate = [str(stoichiometry) + '*' + rate_expr
                                                for rate_expr in [r_expr, f_expr]]
        # Adsorbate is produced.
        else:
            if stoichiometry == 1:
                increase_rate, decrease_rate = f_expr, r_expr
            else:
                increase_rate, decrease_rate = [str(stoichiometry) + '*' + rate_expr
                                                for rate_expr in [f_expr, r_expr]]

        # Return.
        return increase_rate + ' - ' + decrease_rate
        # }}}

    def get_adsorbate_dtheta_dt_expression(self, adsorbate_name):
        """
        Function to get dtheta/dt expression over all elementary reactions wrt an adsorbate.

        Parameters:
        -----------
        adsorbate_name: The name of the adsorbate whose dtheta/dt expression
                        would be returned.

        Returns:
        --------
        dtheta_dt_expression: String of dtheta/dt expression.
        """
        # {{{
        # Collect all dtheta/dt exprression for an elementary reaction.
        dtheta_dt_expression_list = []

        for rxn_expression in self._owner.rxn_expressions():
            single_dtheta_dt = self.get_elementary_dtheta_dt_expression(adsorbate_name,
                                                                        rxn_expression)
            if single_dtheta_dt:
                dtheta_dt_expression_list.append(single_dtheta_dt)

        # Join the expressions.
        dtheta_dt_expression = ' + '.join(dtheta_dt_expression_list)

        return dtheta_dt_expression
        # }}}

    def get_dtheta_dt_expressions(self):
        """
        Function to get dtheta/dt expression strings for all adsorbatets.
        """
        # {{{
        dtheta_dt_expressions_list = []
        for idx, adsorbate_name in enumerate(self._owner.adsorbate_names()):
            dtheta_dt_expression = "dtheta_dt[" + str(idx) + "] = "
            dtheta_dt_expression += self.get_adsorbate_dtheta_dt_expression(adsorbate_name)
            dtheta_dt_expressions_list.append(dtheta_dt_expression)

        dtheta_dt_expressions_tup = tuple(dtheta_dt_expressions_list)

        return dtheta_dt_expressions_tup
        # }}}

    def steady_state_function(self, cvgs_tuple, relative_energies=None):
        """
        Recieve a coverages tuple containing coverages of adsorbates,
        return a tuple of dtheta_dts of corresponding adsorbates.
        """
        # {{{
        # Set theta, kf, kr, p, dtheta_dt
        # Coverages(theta).
        theta = self._cvg_tuple2dict(cvgs_tuple)

        # Rate constants(kf, kr).
        kf, kr = self.get_rate_constants(relative_energies)

        # Pressure.
        p = self._p

        # Concentration.
        c = self._c

        # Rate of coverage change(dtheta_dt).
        dtheta_dt = [0.0]*len(self._owner.adsorbate_names())

        dtheta_dt_expressions = '\n'.join(self.get_dtheta_dt_expressions())
        exec dtheta_dt_expressions in locals()

        return tuple(dtheta_dt)
        # }}}

    @staticmethod
    def __term_adsorbate_derivation(adsorbate_name, term_expression):
        """
        Expect a single expression and an adsorbate_name
        e.g. "kf[2]*theta['CO_s']*theta['*_s']" 'CO_s',
        return a derivation expression wrt adsorbate_name.
        """
        # {{{
        # Escape.
        if '*' in adsorbate_name:
            adsorbate_name = '\\' + adsorbate_name

        regex = r"((.*)\*|)(theta\['" + adsorbate_name + r"'\])(\*{2}(\d)|)(\*(.*)|)"
        #r"(.*)\*(theta\['CO_s'\])(\*\*(\d)|)(\*(.*)|)"
        ###########################################################
        # group(1) -> ((.*)\*|), group(2) -> (.*)                 #
        # group(3) -> (theta\['"+adsorbate_name+"'\])             #
        # group(4) -> (\*{2}(\d)|), group(5) -> \*{2}(\d) or None #
        # group(6) -> (\*(.*)|), group(7) -> (.*) or None         #
        ###########################################################

        # Coefficient
        m = re.search(regex, term_expression)
        if m.group(7) and m.group(2):
            coefficient = m.group(2) + '*' + m.group(7)
        elif m.group(2):
            coefficient = m.group(2)
        elif m.group(7):
            coefficient = m.group(7)
        else:
            coefficient = '1'
        # Power of cvg
        if m.group(4):
            power = int(m.group(5))
        else:
            power = 1

        if power == 1:  # e.g. "kf[2]*theta['CO_s']*theta['*_s']"
            derivation_expression = coefficient
        else:
            if coefficient == '1':
                coefficient = str(power)
            else:
                coefficient = str(power) + '*' + coefficient
            cvg_expression = m.group(3) + '**' + str(power-1)
            derivation_expression = coefficient + '*' + cvg_expression
#        else:
#            derivation_expression = '0'

        return derivation_expression
        # }}}

    def __total_term_adsorbate_derivation(self,
                                          adsorbate_name,
                                          term_expression):
        """
        Private function to get derivation expression taking FREE SITE into consideration.

        NOTE: the coverage of free site can be expressed as (1 - theta_CO_s - ...),
              so the derivation must take coverages of free site into consideration.

        Parameters:
        -----------
        adsorbate_name: The adsorbate name which derivation expression wrt
                        would be returned, str.
        term_expression: The string of term expression, str.

        Returns:
        --------
        The derivation expression, str.
        """
        # {{{
        if adsorbate_name not in self._owner.adsorbate_names():
            raise ValueError("'" + adsorbate_name + "' is not in adsorbate_names")

        def theta(sp_name):
            return "theta['{}']".format(sp_name)

        # Get all sites in term_expression.
        site_cvg_regex = r"theta\['\*_(\w*)'\]"
        sites_list = re.findall(site_cvg_regex, term_expression)

        # Get site info of adsorbate.
        formula = ChemFormula(adsorbate_name)
        site_name = formula.site()
        site_cvg_expr = theta('*_' + site_name)
        site_total = self._owner.species_definitions()[site_name]['total']

        # Get derivation expression wrt free site.
        def deriv_site_part(site_name, term_expression):
            initial_expr = self.__term_adsorbate_derivation('*_'+site_name, term_expression)

            # Convert site expression to adsobate expression.
            if site_cvg_expr in initial_expr:
                # Get substitute expression.
                substitute_expr = str(site_total)
                for adsorbate_name in self._classified_adsorbates[site_name]:
                    substitute_expr += ' - ' + theta(adsorbate_name)
                substitute_expr = '(' + substitute_expr + ')'

                # Do substitution
                site_cvg_regex = r"theta\['\*_" + site_name + r"'\]"
                final_expr = re.sub(site_cvg_regex, substitute_expr, initial_expr)

                # Add minus.
                final_expr = "-" + final_expr
            else:
                # Just add a minus before expression.
                final_expr = '-' + initial_expr

            return final_expr

        # Get derivation expression wrt adsorbate.
        def deriv_adsorbate_part(adsorbate_name, term_expression):
            return self.__term_adsorbate_derivation(adsorbate_name, term_expression)

        # If contains both.
        if (site_name in sites_list) and (theta(adsorbate_name) in term_expression):
            #split two parts
            regex = r"(.*)\*(theta\['" + adsorbate_name + r"'\])(\*{2}(\d)|)(\*(.*)|)"
            m = re.search(regex, term_expression)
            if m.group(6):
                site_part = m.group(1) + '*' + m.group(6)
            else:
                site_part = m.group(1)
            adsorbate_part = m.group(2) + m.group(3)
            # Get derivation expression
            derivation_expression = (deriv_adsorbate_part(adsorbate_name, adsorbate_part) +
                                     '*' + site_part + ' + ' +
                                     deriv_site_part(site_name, site_part) +
                                     '*' + adsorbate_part)
        # Derive empty site coverage only.
        elif site_name in sites_list:
            derivation_expression = deriv_site_part(site_name, term_expression)
        # Derive adsorbate coverage only.
        elif theta(adsorbate_name) in term_expression:
            derivation_expression = deriv_adsorbate_part(adsorbate_name, term_expression)
        else:
            derivation_expression = '0'

        return derivation_expression
        # }}}

    def poly_adsorbate_derivation(self, adsorbate_name, poly_expression):
        """
        Expect a polynomial expression of dtheta_dt and an adsorbate_name,
        return a derivation expression about the adsorbate.
        Function to the derivation expression wrt an adsorbate.

        Parameters:
        -----------
        adsorbate_name: The adsorbate name, str.
        poly_expression: A polynomial expression of dtheta/dt, str.

        Returns:
        --------
        The derivation expression, str.

        Example:
        --------
        >>> adsorbate = "CO_s"
        >>> poly_expression = "dtheta_dt[0] = kf[0]*p['CO_g']*theta['*_s'] - kr[0]*theta['CO_s']"
        >>> solver.poly_adsorbate_derivation(adsorbate, poly_expression)
        """
        # {{{
        # Split poly_expression.
        poly_list = poly_expression.split()
        operators, terms = poly_list[3::2], poly_list[2::2]

        # Generate derived term expressions.
        derived_terms = [self.__total_term_adsorbate_derivation(adsorbate_name,
                                                                term_expression)
                         for term_expression in terms]

        # Combine 2 lists.
        derived_poly_list = []
        for combination in zip(derived_terms, operators):
            derived_poly_list += list(combination)
        derived_poly_list.append(derived_terms[-1])

        # Join and return.
        return ' '.join(derived_poly_list)
        # }}}

    def analytical_jacobian(self, cvgs_tuple, relative_energies=None):
        """
        Function to get analytical Jacobian matrix of the nonlinear equation system.

        Parameters:
        -----------
        cvgs_tuple: A tuple of adsorbate coverages, tuple of float.

        relative_energies: A dict of relative eneriges of elementary reactions.
            NOTE: keys "Gaf" and "Gar" must be in relative energies dict.

        Returns:
        --------
        The analytical Jacobian matrix, N x N matrix of float.
        N is the number of adsorbates.
        """
        # {{{
        # Check input parameter.
        dtheta_dt_expressions = self.get_dtheta_dt_expressions()
        m, n = len(dtheta_dt_expressions), len(cvgs_tuple)

        if m != n:
            msg = "{} coverages are expected, but {} are provided.".format(m, n)
            raise ParameterError(msg)

        # Set theta, kf, kr, p, dtheta_dt
        # Coverages(theta).
        theta = self._cvg_tuple2dict(cvgs_tuple)

        # Rate constants(kf, kr).
        kf, kr = self.get_rate_constants(relative_energies)

        # Pressure.
        p = self._p

        # Concentration.
        c = self._c

        # Generate Jacobian matrix.
        J = self._matrix(m, n)
        adsorbate_names = self._owner.adsorbate_names()

        # Fill matrix.
        for i in xrange(m):
            poly_expression = dtheta_dt_expressions[i]
            for j in xrange(n):
                # Get adsorbate_name.
                adsorbate_name = adsorbate_names[j]
                # Fill the matrix
                derivation = self.poly_adsorbate_derivation(adsorbate_name, poly_expression)
                J[i, j] = eval(derivation)

        # Return Jacobian matrix.
        return J
        # }}}

    ######################################################
    ######                                          ######
    ###### calculate micro kinetic model with Sympy ######
    ######                                          ######
    ######################################################

    def get_elementary_dtheta_dt_sym(self, adsorbate_name, rxn_expression):
        """
        Function to get dtheta/dt expression symbol
        for an adsorbate and an elementary reaction.
        """
        # {{{
        # Check adsorbate name.
        if adsorbate_name not in self._owner.adsorbate_names():
            msg = "'{}' is not an adsorbate".format(adsorbate_name)
            raise ParameterError(msg)

        # Check rxn expression.
        if rxn_expression not in self._owner.rxn_expressions():
            msg = "'{}' not in model's rxn_expressions.".format(rxn_expression)
            raise ParameterError(msg)

        # Get formula list.
        rxn_equation = RxnEquation(rxn_expression)
        elementary_rxn_list = rxn_equation.to_formula_list()

        for state_list in elementary_rxn_list:
            for formula in state_list:
                stoichiometry = formula.stoichiometry()
                species_site = formula.species_site()

                # If find it, jump out.
                if species_site == adsorbate_name:
                    break

            # Get state idx to get direction info.
            if species_site == adsorbate_name:
                state_idx = elementary_rxn_list.index(state_list)
                break

        # If adsorbate name not in elementary_rxn list, stop.
        if species_site != adsorbate_name:
            return

        # Get dtheta_dt sym according rate symbols.
        rf_sym, rr_sym = self.get_single_rate_sym(rxn_expression)
        if state_idx == 0:
            dtheta_dt_sym = -rf_sym + rr_sym
        else:
            dtheta_dt_sym = rf_sym - rr_sym

        return dtheta_dt_sym
        # }}}

    def get_adsorbate_dtheta_dt_sym(self, adsorbate_name):
        """
        Function to get dtheta/dt for an adsorbate.

        Parameters:
        -----------
        adsorbate_name: The adsorbate name, str.
        """
        total_dtheta_dt_sym = 0

        # Loop over all elementary reactions.
        for rxn_expression in self._owner.rxn_expressions():
            dtheta_dt_sym = self.get_elementary_dtheta_dt_sym(adsorbate_name,
                                                              rxn_expression)
            # Adsorbate not in this reaction expression.
            if dtheta_dt_sym is None:
                continue
            else:
                total_dtheta_dt_sym += dtheta_dt_sym

        return total_dtheta_dt_sym

    def get_dtheta_dt_syms(self, log_latex=False):
        """
        Function to get dtheta/dt expressions for all adsorbates.
        """
        dtheta_dt_syms = []
        for adsorbate_name in self._owner.adsorbate_names():
            dtheta_dt_sym = self.get_adsorbate_dtheta_dt_sym(adsorbate_name)
            dtheta_dt_syms.append(dtheta_dt_sym)

        dtheta_dt_syms = tuple(dtheta_dt_syms)

        # Latex strings.
        dtheta_dt_latexs = self.get_latex_strs(part1=r'\frac{d\theta_{', part2=r'}}{dt}} ',
                                               symbols=dtheta_dt_syms)

        # Log it.
        if log_latex:
            self.log_latex(self.dtheta_dt_latex)

        return dtheta_dt_syms

    def steady_state_function_by_sym(self, cvgs_tuple):
        """
        Recieve a coverages tuple containing coverages of adsorbates,
        return a tuple of dtheta_dts of corresponding adsorbates.
        """
        # Get dtheta/dt expressions.
        dtheta_dt_syms = self.get_dtheta_dt_syms()

        # Get substitution dict.
        subs_dict = self.get_subs_dict(coverages=cvgs_tuple)

        # Loop to get values of dtheta/dt.
        dtheta_dts = []
        for dtheta_dt_sym in dtheta_dt_syms:
            dtheta_dt = self._mpf(dtheta_dt_sym.evalf(subs=subs_dict))
            dtheta_dts.append(dtheta_dt)

        return tuple(dtheta_dts)

    def analytical_jacobian_sym(self):
        """
        Function to get the jacobian matrix symbol expressions of
        the dtheta/dt nonlinear equations.
        """
        # Get dtheta/dt expressions.
        dtheta_dt_syms = self.get_dtheta_dt_syms()

        # Allocate memories for jacobian matrix.
        m = n = len(dtheta_dt_syms)
        sym_jacobian = [[0.0]*m, [0.0]*n]

        # dtheta/dt (row).
        for i in xrange(m):
            dthe_dt_sym = dtheta_dt_syms[i]
            # adsorbate (column).
            for j in xrange(n):
                ads_name = self._owner.adsorbate_names()[j]
                theta_sym = self._extract_symbol(ads_name, 'ads_cvg')
                sym_jacobian[i][j] = sym.Derivative(dthe_dt_sym, theta_sym).doit()

        return sym_jacobian

    def analytical_jacobian_by_sym(self, cvgs_tuple):
        """
        Get the jacobian matrix of the dtheta/dt nonlinear equations.
        Return a jacobian matrix(in self._matrix form).
        """
        # NOTE: precision may lose here.

        # Get substitution dicts.
        subs_dict = self.get_subs_dict(coverages=cvgs_tuple)

        # Get symbol jacobian matrix.
        sym_jacobian = self.analytical_jacobian_sym()

        m = n = len(sym_jacobian)
        num_jacobian = [[0.0]*m, [0.0]*n]

        for i in xrange(m):
            for j in xrange(n):
                entry_value = sym_jacobian[i][j].evalf(subs=subs_dict)
                num_jacobian[i][j] = self._mpf(entry_value)

        return self._matrix(num_jacobian)

#    def get_rate_control_by_sym(self, RDS):
#        """
#        RDS: int, Rate Determining Step number.
#        """
#        #get quasi_quilibrium_solver instance
#        _temp = __import__('quasi_equilibrium_solver',
#                           globals(), locals(), ['QuasiEquilibriumSolver'])
#        qe_solver = _temp.QuasiEquilibriumSolver(owner=self._owner)
#        qe_solver.RDS = RDS  # set Rate Determining Step
#        XTRCs = qe_solver.get_XTRCs()
#        self.qe_solver = qe_solver
#
#        return XTRCs

    ##########################################################
    ###### calculate micro kinetic model with Sympy END ######
    ##########################################################

    def get_residual(self, cvgs_tuple, relative_energies=None):
        """
        Function to get residual value of equations(the max value of dthe/dt).

        Parameters:
        -----------
        cvgs_tuple: A tuple of adsorbate coverages.

        relative_energies: A dict of relative eneriges of elementary reactions.
            NOTE: keys "Gaf" and "Gar" must be in relative energies dict.

        Returns:
        --------
        The max value of dtheta/dt wrt the coverages.
        """
        #constrain cvgs
        #cvgs_tuple = self.__constrain_coverages(cvgs_tuple)
        dtheta_dts = self.steady_state_function(cvgs_tuple, relative_energies)
        residual = max([abs(dtheta_dt) for dtheta_dt in dtheta_dts])

        return residual

    # -- Use scipy.optimize.fsolve function to solve equations --

    def coarse_steady_state_cvgs(self, c0, relative_energies=None):
        '''
        Use scipy.optimize.fsolve to solve non-linear equations
        with fast speed and low-precison.

        Parameters:
        -----------
        c0: initial coverages, a list or tuple of float

        relative_energies: A dict of relative eneriges of elementary reactions.
            NOTE: keys "Gaf" and "Gar" must be in relative energies dict.

        Return:
        -------
        steady_state_coverages: in the order of self._owner.adsorbate_names,
                                tuple of float
        '''

        def get_jacobian(c0, relative_energies=None):
            # Jacobian matrix.
            jm = self.analytical_jacobian(c0, relative_energies).tolist()
            # Convert to floats.
            jm = [[float(df) for df in dfs] for dfs in jm]

            return jm

        #fprime = lambda x: get_jacobian(c0, relative_energies=relative_energies)

        # Main hotpot.
        c0 = map(float, c0)  # covert to float
        converged_cvgs = fsolve(self.steady_state_function, c0,
                                (relative_energies, ),
                                fprime=get_jacobian)

        return converged_cvgs

    def data_archive(func):
        '''
        Decorator for **_steady_state_coverages,
        archive data when getting converged converages.
        '''
        def wrapped_func(*args, **kwargs):
            '''
            Use scipy.optimize.fsolve to solve non-linear equations.

            Parameters:
            -----------
            c0: initial coverages, a list or tuple of float

            Return:
            -------
            steady_state_coverages: in the order of self._owner.adsorbate_names,
                                    tuple of float
            '''
            converged_cvgs = func(*args, **kwargs)
            # Archive data.
            # Get error.
            self = args[0]
            errors = self.steady_state_function(converged_cvgs, kwargs["relative_energies"])
            error = norm(errors)
            self._error = error

            self._coverage = converged_cvgs
            # log steady state coverages.
            self.__log_sscvg(converged_cvgs, self._owner.adsorbate_names())
            self.__logger.info('error = %e', error)

            # Archive converged root and error.
            self.archive_data('steady_state_coverage', converged_cvgs)
            self.archive_data('steady_state_error', error)
            self.good_guess = args[1]
            # Archive initial guess.
            self.archive_data('initial_guess', args[1])

            return converged_cvgs

        return wrapped_func

    @data_archive
    def fsolve_steady_state_cvgs(self, c0, relative_energies=None):
        """
        Use scipy.optimize.fsolve to get steady state coverages.
        """
        return self.coarse_steady_state_cvgs(c0, relative_energies)

    # -- fsolve END --

    def get_steady_state_cvgs(self, c0, single_pt=False, relative_energies=None):
        """
        Function to get steady state coverages.

        Expect an inital coverages tuple, use Newton Method to
        solving nonlinear equations, return steady state coverages,
        if converged.

        Parameters
        ----------
        c0: initial coverages, tuple of float.

        single_pt : bool
            if True, no initial guess check will be done.

        relative_energies: Optional, dict of relative eneriges of elementary reactions.
            NOTE: keys "Gaf" and "Gar" must be in relative energies dict.

        """
        # Intial coverage must have physical meaning.
        c0 = self.__constrain_coverages(c0)
        self.__initial_guess = c0

        # Start root finding algorithm.
        f = lambda x: self.steady_state_function(x, relative_energies=relative_energies)
        f_resid = lambda x: self.get_residual(x, relative_energies=relative_energies)
        constraint = self.__constrain_coverages
        J = lambda x: self.analytical_jacobian(x, relative_energies=relative_energies)

        ############    Main Loop with changed initial guess   ##############
        self.__logger.info('Entering main loop...')
        icvg_counter = 1  # initial coverage counter, outer
        cancel = False

        converged = False  # Flag for convergence.
        while not cancel:  # outer loop
            # Good initial coverages.
            if f_resid(c0) <= self._tolerance and not single_pt:
                converged = True
                self._coverage = c0
                self.__logger.info('Good initial guess: \n%s', str(map(float, c0)))

                # log steady state coverages
                self.__log_sscvg(c0, self._owner.adsorbate_names())

                # Get error.
                fx = self.steady_state_function(c0, relative_energies)  # dtheta/dts
                norm = self._norm(fx)
                resid = self.get_residual(c0)
                error = min(norm, resid)
                self._error = error
                self.__logger.info('error = %e', error)
                break

            # Instantiate rootfinding iterator
            # ConstrainedNewton iterator
            if self._rootfinding == 'ConstrainedNewton':
                iterator_parameters = dict(J=J,
                                           constraint=constraint,
                                           norm=self._norm,
                                           mpfloat=self._mpf,
                                           matrix=self._matrix,
                                           Axb_solver=self._Axb_solver)
                newton_iterator = ConstrainedNewton(f, c0, **iterator_parameters)
            # MDNewton iterator
            elif self._rootfinding == 'MDNewton':
                iterator_parameters = dict(J=J, verbose=False)
                newton_iterator = MDNewton(f, c0, **iterator_parameters)
            else:
                msg='Unrecognized rootfinding iterator name [{}]'.format(self._rootfinding)
                raise ParameterError(msg)

            self.__logger.info('{} Iterator instantiation - success!'.format(self._rootfinding))

            x = c0
            old_error = 1e99
            if c0:
                # log initial guess
                self.__logger.info('initial guess coverage - success')
                self.__logger.debug(str(map(float, c0)))

            #####    Sub LOOP for a c0    #####

            nt_counter = 0    # newton loop counter, inner
            self.__logger.info('entering Newton Iteration( %d )...', icvg_counter)
            # log title
            self.__logger.info('  %-10s   %5s  %18s  %18s',
                               'status', 'N', 'residual', 'norm')
            self.__logger.info('-'*60)

            for x, error, fx in newton_iterator:  # inner loop
                nt_counter += 1
                resid = f_resid(x)
                self.__logger.info('%-10s%10d%23.10e%23.10e', 'in_process',
                                   nt_counter, float(resid), float(error))
                # Less than tolerance
                if error < self._tolerance:
                    if resid < self._tolerance:
                        # Check whether there is minus value in x
                        for cvg in x:
                            if cvg < 0.0:
                                lt_zero = True  # less than 0
                                break
                            else:
                                lt_zero = False
                        # check END #
                        if not lt_zero:
                            self.__logger.info('%-10s%10d%23.10e%23.10e', 'success',
                                               nt_counter, float(resid), float(error))
                            # log steady state coverages
                            self.__log_sscvg(x, self._owner.adsorbate_names())
                            self._coverage = x
                            self._error = error
                            self.__logger.info('error = %e', min(error, resid))

                            # Update flags.
                            cancel = True
                            converged = True
                            break
                        else:  # bad root, iteration continue...
                            self.__logger.warning('bad root: %s', str(map(float, x)))
                            self.__logger.warning('root finding continue...\n')
                    else:
                        error = f_resid(x)  # use residual as error and continue

                # if convergence is slow when the norm is larger than 0.1
                elif (nt_counter > self._max_rootfinding_iterations or
                      abs(error - old_error) < 1e-4) and error > 1e-10:
                    self.__logger.info('%-10s%10d%23.10e%23.10e', 'break',
                                     nt_counter, float(resid), float(error))
                    self.__logger.warning('slow convergence rate!')
                    self.__logger.warning('root finding break for this initial guess...\n')
                    # Jump out of loop for this c0
                    cancel = False
                    break

                # residual is almost stagnated
                elif abs(error - old_error) < self._stable_criterion:
                    self.__logger.info('%-10s%10d%23.10e%23.10e', 'stable',
                                       nt_counter, float(resid), float(error))
                    self.__logger.warning('stable root: %s', str(map(float, x)))
                    self.__logger.debug(' difference: %-24.16e', abs(error - old_error))
                    # Jump out of loop for this c0.
                    cancel = False
                    break

                old_error = error  # set old error to be compared in next loop
                #self._coverage = x
                #self._error = error

                # archive data every 100 steps
                if nt_counter % 100 == 0:
                    self.archive_data('iter_coverage', x)
                    self.archive_data('iter_error', error)
            #####    Sub loop for a c0 END    #####

            # Change the initial guess(c0).
            if not cancel:
                # get a new initial guess coverage
                c0 = self.modify_init_guess(x, fx)
                icvg_counter += 1

        ##############    main loop end   #################

        if converged:
            # Archive converged root and error.
            self.archive_data('steady_state_coverage', self._coverage)
            self.archive_data('steady_state_error', self._error)
            self.good_guess = c0

            # Archive initial guess.
            self.archive_data('initial_guess', c0)

            return self._coverage

    def __log_sscvg(self, cvgs_tuple, ads_names):
        """
        Private helper function to log steady state coverage of every species.
        """
        head_str = "\n {:<10s}{:<25s}{:<30s}\n".format("index",
                                                       "intermediate name",
                                                       "steady state coverage")
        head_str = "Steady State Coverages:\n" + head_str
        line_str = '-'*60 + '\n'

        all_data = ''
        all_data += head_str + line_str
        for idx, (ads_name, cvg) in enumerate(zip(ads_names, cvgs_tuple)):
            idx = str(idx).zfill(2)
            data = " {:<10s}{:<25s}{:<30.16e}\n".format(idx, ads_name, float(cvg))
            all_data += data
        all_data += line_str

        self.__logger.info(all_data)

        return all_data

    def __get_Gs_tof(self, Gs, gas_name=None):  # Gs -> free energies
        """
        Private function to get TOF for a given formation energies.

        Parameters:
        -----------
        Gs: free energies of all intermediates(adsorbates and transition states).

        gas_name: The gas whose TOF would be returned.

        Returns:
        --------
        If gas name is not specified, a list of TOF for all gas species would be returned.
        If gas name is specified, the TOF of the gas would be returned.
        """
        # Get initial guess
        if hasattr(self, "_coverage"):
            init_guess = self._coverage
        else:
            msg = ("Converged coverages are needed to calculate TOF for a list " +
                   "formation energies, so try to get steady state coverages first.")
            raise AttributeError(msg)

        Gs_order = self._owner.adsorbate_names() + self._owner.transition_state_names()

        # Copy the original energies.
        G_copy = copy.deepcopy(self._G)

        # Update formation energies of solver.
        for intermediate, G in zip(Gs_order, Gs):
            self._G[intermediate] = G

        # Get new relative energies.
        relative_energies = self.get_relative_from_absolute()

        # Calculate the new steady state coverages.
        steady_state_cvg = self.get_steady_state_cvgs(init_guess, relative_energies)

        # Get turnover frequencies.
        tof_list = self.get_tof(steady_state_cvg, relative_energies)

        # Recover solver's formation energy dict.
        self._G = G_copy

        # Return.
        if gas_name is None:
            return tof_list
        else:
            # Check gas name.
            gas_names = self._owner.gas_names()
            if gas_name not in gas_names:
                msg = "'{}' is not a gas species in model".format(gas_name)
                raise ParameterError(msg)
            idx = gas_names.index(gas_name)
            return tof_list[idx]

    def __get_intermediates_Gs(self):
        """
        Private helper function to get formation energies of intermediates.
        """
        all_intermediates = (self._owner.adsorbate_names() +
                             self._owner.transition_state_names())
        Gs = []
        for intermediates_name in all_intermediates:
            Gs.append(self._G[intermediates_name])

        return Gs

    def get_rate_control(self, epsilon=None):
        """
        Function to get DTRC matrix for all gas species.

        Parameters:
        -----------
        epsilon: The perturbation size for numerical jacobian matrix.

        Returns:
        --------
        An DTRC mpmath.matrix for all gas species.
            - rows for gas species.
            - columns for intermediates and transition states.
        """
        # Get intermediates formation energies.
        Gs = self.__get_intermediates_Gs()

        kT = self._owner.kB()*self._owner.temperature()

        # Get perturbation size.
        if epsilon is None:
            epsilon = self._mpf(self._perturbation_size)

        # Get dr/dG matrix.
        drdG = numerical_jacobian(f=self.__get_Gs_tof,
                                  x=Gs,
                                  num_repr=self._numerical_representation,
                                  matrix=self._matrix,
                                  h=epsilon,
                                  direction=self._perturbation_direction)
        r = self.__get_Gs_tof(Gs)

        # multiply 1/r to drdG matrix.
        diag_matrix = self._linalg.diag([-kT/tof for tof in r])
        DTRC = diag_matrix*drdG

        # Covert it to list.
        DTRC_list = DTRC.tolist()

        # Archive
        self.archive_data('DTRC', DTRC_list)

        # Log it.
        self.__log_DTRC(DTRC_list)

        return DTRC

    def __log_DTRC(self, DTRC_matrix):
        """
        Private helper function to log DTRC for all gas species.
        """
        gas_names = self._owner.gas_names()
        intermediate_names = (self._owner.adsorbate_names() +
                              self._owner.transition_state_names())

        head_str = "\n{:<15s}{:<30s}{:<30s}\n".format("gas", "intermediate", "DTRC")
        head_str = "Degree of Thermodynamic Rate Control:\n" + head_str
        line_str = "-"*70 + "\n"

        all_data = ""
        all_data += head_str + line_str

        for gas_name, DTRC_vect in zip(gas_names, DTRC_matrix):
            for intermediate, DTRC in zip(intermediate_names, DTRC_vect):
                data = "{:<15s}{:<30s}{:<30.16e}\n".format(gas_name, intermediate, float(DTRC))
                all_data += data

        all_data += line_str

        self.__logger.info(all_data)

        return all_data

    def modify_init_guess(self, *args):
        '''
        return a list of random coverages.
        '''
        n_adsorbates = len(self._owner.adsorbate_names())
        random_cvgs = []

        sum_cvgs = 0.0
        for i in xrange(n_adsorbates):
            cvg = random.random()*(1.0 - sum_cvgs)
            random_cvgs.append(cvg)
            sum_cvgs += cvg
        #add to log
        self.__logger.info('modify initial coverage - success')
        self.__logger.debug(str(random_cvgs))

        return tuple(random_cvgs)

    def modify_init_guess_old(self, c0, dtheta_dts):
        "Return a new initial guess according to dthe_dts."
#        max_dtheta_dt = np.max(np.abs(dtheta_dts))
        base_coefficient = self._initial_guess_scale_factor
        coefficients = []
        for dtheta_dt in np.abs(dtheta_dts):
            if abs(dtheta_dt) >= self._tolerance:
                #coefficients.append(dtheta_dt/max_dtheta_dt*base_coefficient)
                coefficients.append(base_coefficient)
            else:
                coefficients.append(1.0)
        self.__logger.debug('coeff: %s', str(coefficients))
        #if coeffs are all 1.0, break!
        #add later...

        #create a diagnol matrix
        c0_diag = np.matrix(np.diag(c0))
        #convert coeffients to column vector
        coefficients = np.matrix(coefficients).reshape(-1, 1)
        new_c0 = (c0_diag*coefficients).reshape(1, -1)
        new_c0 = tuple(new_c0.tolist()[0])
        #add to log
        self.__logger.info('modify initial coverage - success')
        self.__logger.debug(str(map(float, c0)))

        #return self.__constrain_coverage__constrain_coverages(new_c0)
        return new_c0

    def modify_init_guess_new(self, c0, dtheta_dts):
        "Return a new initial guess according to dthe_dts."
        max_dtheta_dt = np.max(np.abs(dtheta_dts))
        base_coefficient = self._initial_guess_scale_factor
        coefficients = []
        for idx, dtheta_dt in enumerate(np.abs(dtheta_dts)):
            if abs(dtheta_dt) >= self._tolerance:
                if dtheta_dt < 0:
                    coefficients.append(dtheta_dt/max_dtheta_dt*base_coefficient)
                elif dtheta_dt > 0:
                    coefficients.append(dtheta_dt/max_dtheta_dt/base_coefficient)
                #coefficients.append(base_coefficient)
            else:
                coefficients.append(1.0)
        print coefficients
        #create a diagnol matrix
        c0_diag = np.matrix(np.diag(c0))
        #convert coeffients to column vector
        coefficients = np.matrix(coefficients).reshape(-1, 1)
        new_c0 = (c0_diag*coefficients).reshape(1, -1)
        new_c0 = tuple(new_c0.tolist()[0])
        #add to log
        self.__logger.info('modify initial coverage - success')
        self.__logger.debug(str(map(float, c0)))

        return new_c0

    ####################################
    ## solve model by ODE integration ##
    ####################################

    def solve_ode(self, algo='lsoda', time_start=0.0, time_end=100.0,
                  time_span=0.1, initial_cvgs=None, relative_energies=None):
        """
        Solve the differetial equations, return points of coverages.

        Parameters:
        -----------
        algo: algorithm for ODE solving, optional, str.
              'vode' | 'zvode' | 'lsoda' | 'dopri5' | 'dop853'

        time_span: time span for each step, float, default to be 0.1

        time_start: time when begin integration, float.

        time_end: time when stop integration, float.

        initial_cvgs: initial coverages at time_start, tuple of float

        relative_energies: A dict of relative eneriges of elementary reactions.
            NOTE: keys "Gaf" and "Gar" must be in relative energies dict.

        Returns:
        --------
        ts: time points, list of float.
        ys: integrated function values, list of list of float.

        Examples:
        ---------
        >>> m.solver.solve_ode(initial_cvgs=(0.0, 0.0))

        """
        # set timr variables
        t_start = time_start
        t_end = time_end
        t_step = time_span

        adsorbate_names = self._owner.adsorbate_names()
        nads = len(adsorbate_names)

        # set initial points
        if not initial_cvgs:
            try:
                initial_cvgs = self.boltzmann_coverages()
            except IOError:
                initial_cvgs = [0.0]*nads

        # differential equation, solve over t for initial coverages cvgs_tuple
        def f(t, cvgs_tuple):
            return list(self.steady_state_function(cvgs_tuple, relative_energies))

        # ode solver object
        r = ode(f)
        r.set_integrator(algo, method='bdf')
        r.set_initial_value(initial_cvgs, t_start)

        ts, ys = [], []

        # integration loop
        self.__logger.info('entering %s ODE integration loop...\n', algo)
        self.__logger.info('%10s%20s' + '%20s'*nads, 'process',
                           'time(s)', *adsorbate_names)
        self.__logger.info('-'*(20*nads + 30))

        try:
            while r.t < t_end:
                self.__logger.info('%10.2f%%%20f' + '%20.8e'*nads, r.t/t_end*100,
                                   r.t, *r.integrate(r.t + t_step))
                ts.append(r.t)
                ys.append(r.y.tolist())
            self.__logger.info('%10s\n', 'stop')

        finally:
            # write to archive file
            with open('auto_ode_coverages.py', 'w') as f:
                time_str = get_list_string('times', ts)
                cvgs_str = get_list_string('coverages', ys)
                content = file_header + time_str + cvgs_str
                f.write(content)
                self.__logger.info('ODE integration trajectory is written' +
                                   ' to auto_ode_coverages.py.')

        return ts, ys

    def error(self):
        """
        Query function for converged error.
        """
        return self._error

    def coverages(self):
        """
        Query function for converaged coverages.
        """
        return self._coverages


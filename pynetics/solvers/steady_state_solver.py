import re

from solver_base import *


class SteadyStateSolver(SolverBase):
    def __init__(self, owner):
        SolverBase.__init__(self, owner)

        #set default parameter dict
        defaults = dict(
            tolerance=1e-8,
            max_rootfinding_iterations=100,
            residual_threshold=1.0,
            initial_guess_scale_factor=100,
            stable_criterion=1e-10,
            )
        defaults = self.update_defaults(defaults)
        self.__dict__.update(defaults)

        #update steady_state_solver's own logger template
        self.logger_template_dict = {
            'coverages_constraint_warning': '${initial_cvgs} --> ' +
                                            '${constrained_cvgs}\n',

                    'inital_guess_success': 'initial point (${initial_cvgs})\n',
                  'rootfinding_in-process': 'iteration ${n_iter} with ' +
                                            'residual ${resid}, norm ${norm}, '
                                            '\n' + ' '*38 + 'dtheta/dt ${fx}' +
                                            '\n' + ' '*38 + 'coverage ${x}',

           'rootfinding_in-process(brief)': 'iteration ${n_iter} with ' +
                                            'residual ${resid}, norm ${norm}',

                     'rootfinding_success': 'iteration ${n_iter} with residual ' +
                                            '${resid}\n\nsteady_state_root:\n${root}',

                        'rootfinding_fail': '\niteration ${n_iter} with ' +
                                            'residual ${resid}, norm ${norm}',

                       'rootfinding_break': 'iteration ${n_iter} with residual ' +
                                            '${resid}, \nstable_root (${root})\n',

                     'rootfinding_maxiter': '\nout of iteration at iteration ' +
                                            '{$n_iter} - residual ${resid}',

            'modify_initial_guess_success': '\n${old_init_cvg} --> ' +
                                            '${new_init_cvg}\n',

                    'rootfinding_continue': 'Bad steady_state_root : ${bad_root}, ' +
                                            'iteration continue.',

                      'rootfinding_stable': 'iteration ${n_iter} with residual ' +
                                            '${resid}, \nstable_root (${root})\n',
        }
        self.logger._templates_dict.update(self.logger_template_dict)

    def cvg_tuple2dict(self, cvgs_tuple):
        "Convert coverages list to coverages dict."
        #there are some error when convert tuple to dict
        #which is so small that we can ignore it

        #create cvgs_dict containing adsorbates
        cvgs_dict = {}
        for adsorbate_name in self._owner.adsorbate_names:
            idx = self._owner.adsorbate_names.index(adsorbate_name)
            cvgs_dict.setdefault(adsorbate_name, cvgs_tuple[idx])
        #add free site coverages
        for site_name in self._owner.site_names:
            total_cvg = \
                self._mpf(self._owner.species_definitions[site_name]['total'])
            #print total_cvg
            sum_cvg = 0.0
            for sp in self.classified_adsorbates[site_name]:
                #print cvgs_dict[sp]
                sum_cvg += cvgs_dict[sp]
            free_site_cvg = total_cvg - sum_cvg
            cvgs_dict.setdefault('*_'+site_name, free_site_cvg)
        return cvgs_dict

    def cvg_dict2tuple(self, cvgs_dict):
        "Convert coverages dict to coverages list."
        cvgs_list = []
        for adsorbate_name in self._owner.adsorbate_names:
            cvgs_list.append(cvgs_dict[adsorbate_name])
        return tuple(cvgs_list)

    def constrain_converage(self, cvgs_tuple):
        """
        Constrain coverages of absorbates between 0.0 and 1.0/total number.
        """
        #convert tuple to dict
        cvgs_dict = self.cvg_tuple2dict(cvgs_tuple)
        #enforce explicit maxima, cannot be larger than 1.0, smaller than 0.0
        for adsorbate_name in self._owner.adsorbate_names:
            if cvgs_dict[adsorbate_name] > 1.0:
                cvgs_dict[adsorbate_name] = self._mpf('1.0')
            if cvgs_dict[adsorbate_name] < 0.0:
                cvgs_dict[adsorbate_name] = self._mpf('0.0')

        #enforce explicit maxima, cannot be larger than site's total number
        for site_name in self._owner.site_names:
            total_cvg = self._owner.species_definitions[site_name]['total']
            for adsorbate_name in self.classified_adsorbates[site_name]:
                if cvgs_dict[adsorbate_name] > total_cvg:
                    cvgs_dict[adsorbate_name] = self._mpf(total_cvg)

        def constrain_in_total(cvgs_dict, max_cvg):
            "Make sure sum of cvgs in cvgs_dict is not larger than max_cvg"
            total_cvg = sum(cvgs_dict.values())
            if total_cvg > max_cvg:
                for key in cvgs_dict:
                    cvgs_dict[key] = cvgs_dict[key]/total_cvg*max_cvg
            return cvgs_dict

        #enforce site conservation,
        #sum of cvgs on one type of surface <= site total e.g 1.0
        for site_name in self._owner.site_names:
            max_cvg = self._owner.species_definitions[site_name]['total']
            sub_cvgs_dict = {}

            #add cvgs of the site to sub_cvgs_dict
            for adsorbate_name in self.classified_adsorbates[site_name]:
                sub_cvgs_dict.setdefault(adsorbate_name,
                                         cvgs_dict[adsorbate_name])
            #add free site coverage
#            sub_cvgs_dict.setdefault('*_'+site_name, cvgs_dict['*_'+site_name])
            sub_cvgs_dict = constrain_in_total(sub_cvgs_dict, max_cvg)
            cvgs_dict.update(sub_cvgs_dict)
        #convert dict to tuple, and return
        constrained_cvgs_tuple = self.cvg_dict2tuple(cvgs_dict)
        #log if constraint has been carried out
        if constrained_cvgs_tuple != cvgs_tuple:
            self.logger.log(
                log_type='event',
                event='coverages_constraint_warning',
                initial_cvgs=tuple(map(float, cvgs_tuple)),
                constrained_cvgs=tuple(map(float, constrained_cvgs_tuple))
            )
        return constrained_cvgs_tuple

    def get_elementary_dtheta_dt_expression(self, adsorbate_name,
                                            elementary_rxn_list):
        """
        Expect elementary_rxn_list and an adsorbate_name in it,
        return dtheta_dt of the corresponding adsorbate in single
        elementary equation.
        """
        #species must be adsorbate
        if adsorbate_name not in self._owner.adsorbate_names:
            raise ValueError("'"+adsorbate_name+"' is not an adsorbate!")
        for state_list in elementary_rxn_list:
            for species_str in state_list:
                stoichiometry, site_name = self.split_species(species_str)
                if site_name == adsorbate_name:
                    break
            if site_name == adsorbate_name:
                #get state idx to get direction info
                state_idx = elementary_rxn_list.index(state_list)
                break
        #if adsorbate name not in elementary_rxn list, stop
        if site_name != adsorbate_name:
            return
        #get rate expression of the elementary equation
        f_expr, r_expr = \
            self.get_elementary_rate_expression(elementary_rxn_list)
        if state_idx == 0:  # adsorbate is consumed
            if stoichiometry == 1:
                increase_rate, decrease_rate = r_expr, f_expr
            else:
                increase_rate, decrease_rate = \
                    [str(stoichiometry)+'*'+rate_expr
                     for rate_expr in [r_expr, f_expr]]
        else:  # adsorbate is producted
            if stoichiometry == 1:
                increase_rate, decrease_rate = f_expr, r_expr
            else:
                increase_rate, decrease_rate = \
                    [str(stoichiometry)+'*'+rate_expr
                     for rate_expr in [f_expr, r_expr]]
        #return dtheta_dt expression
        return increase_rate + ' - ' + decrease_rate

    def get_adsorbate_dtheta_dt_expression(self, adsorbate_name):
        """
        Expect a adsorbate_name, and go through self.rxns_list,
        return dtheta_dt of the adsorbate.
        """
        dtheta_dt_expression_list = []
        for elementary_rxn_list in self.rxns_list:
            single_dtheta_dt = \
                self.get_elementary_dtheta_dt_expression(adsorbate_name,
                                                         elementary_rxn_list)
            if single_dtheta_dt:
                dtheta_dt_expression_list.append(single_dtheta_dt)
        dtheta_dt_expression = ' + '.join(dtheta_dt_expression_list)

        return dtheta_dt_expression

    def get_dtheta_dt_expressions(self):
        """
        Go through adsorbate_names,
        return a tuple of dtheta_dt_expressions.
        """
        dtheta_dt_expressions_list = []
        for adsorbate_name in self._owner.adsorbate_names:
            adsorbate_idx = self._owner.adsorbate_names.index(adsorbate_name)
            dtheta_dt_expression = "dtheta_dt[" + str(adsorbate_idx) + "] = "
            dtheta_dt_expression += \
                self.get_adsorbate_dtheta_dt_expression(adsorbate_name)
            dtheta_dt_expressions_list.append(dtheta_dt_expression)

        dtheta_dt_expressions_tup = tuple(dtheta_dt_expressions_list)
        setattr(self, 'dtheta_dt_expressions',
                dtheta_dt_expressions_tup)

        return dtheta_dt_expressions_tup

    def steady_state_function(self, cvgs_tuple):
        """
        Recieve a coverages tuple containing coverages of adsorbates,
        return a tuple of dtheta_dts of corresponding adsorbates.
        """
        #set theta, kf, kr, p, dtheta_dt
        #coverages(theta)
        theta = self.cvg_tuple2dict(cvgs_tuple)
        #self.constrain_converage(theta)
        #rate constants(kf, kr)
        kf, kr = self.get_rate_constants()
        #pressure
        p = self.p
        #concentration
        c = self.c
        #rate of coverage change(dtheta_dt)
        dtheta_dt = [0.0]*len(self._owner.adsorbate_names)

        dtheta_dt_expressions = '\n'.join(self.get_dtheta_dt_expressions())
        exec dtheta_dt_expressions in locals()

        return tuple(dtheta_dt)

    @staticmethod
    def term_adsorbate_derivation(adsorbate_name, term_expression):
        """
        Expect a single expression and an adsorbate_name
        e.g. "kf[2]*theta['CO_s']*theta['*_s']" 'CO_s',
        return a derivation expression wrt adsorbate_name.
        """
        if '*' in adsorbate_name:
            adsorbate_name = '\\' + adsorbate_name
#        cvg_name = "theta['"+adsorbate_name.strip('\\')+"']"
#        if cvg_name in term_expression:
        regex = \
            "((.*)\*|)(theta\['"+adsorbate_name+"'\])(\*{2}(\d)|)(\*(.*)|)"
        #r"(.*)\*(theta\['CO_s'\])(\*\*(\d)|)(\*(.*)|)"
        ###########################################################
        # group(1) -> ((.*)\*|), group(2) -> (.*)                 #
        # group(3) -> (theta\['"+adsorbate_name+"'\])             #
        # group(4) -> (\*{2}(\d)|), group(5) -> \*{2}(\d) or None #
        # group(6) -> (\*(.*)|), group(7) -> (.*) or None         #
        ###########################################################
        #coefficient
        m = re.search(regex, term_expression)
        if m.group(7) and m.group(2):
            coefficient = m.group(2) + '*' + m.group(7)
        elif m.group(2):
            coefficient = m.group(2)
        elif m.group(7):
            coefficient = m.group(7)
        else:
            coefficient = '1'
        #for power of cvg
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

    def total_term_adsorbate_derivation(self, adsorbate_name,
                                        term_expression):
        """
        Expect a single expression and an adsorbate_name
        return a total derivation expression about adsorbate_name
        taking free site into consideration.
        """
        if adsorbate_name not in self._owner.adsorbate_names:
            raise ValueError("'"+adsorbate_name+"' is not in adsorbate_names")

        def theta(sp_name):
            return "theta['" + sp_name + "']"

        site_cvg_regex = r"theta\['\*_(\w*)'\]"
        sites_list = re.findall(site_cvg_regex, term_expression)
        site_name = self._owner.species_definitions[adsorbate_name]['site']
        site_cvg_expr = theta('*_' + site_name)
        site_total = self._owner.species_definitions[site_name]['total']

        #get derivation expression wrt free site
        def deriv_site_part(site_name, term_expression):
            initial_expr = \
                self.term_adsorbate_derivation('*_'+site_name, term_expression)
            if site_cvg_expr in initial_expr:
                #convert site expression to adsobate expression
                #get substitute expression
                substitute_expr = str(site_total)
                for adsorbate_name in self.classified_adsorbates[site_name]:
                    substitute_expr += ' - ' + theta(adsorbate_name)
                substitute_expr = '(' + substitute_expr + ')'
                #do substitution
                site_cvg_regex = "theta\['\*_"+site_name+"'\]"
                final_expr = \
                    re.sub(site_cvg_regex, substitute_expr, initial_expr)
            else:
                #just add a minus before expression
                final_expr = '-' + initial_expr
            return final_expr

        #get derivation expression wrt adsorbate
        def deriv_adsorbate_part(adsorbate_name, term_expression):
            return self.term_adsorbate_derivation(adsorbate_name,
                                                  term_expression)

        #if contains both
        if site_name in sites_list and \
           theta(adsorbate_name) in term_expression:
            #split two parts
            regex = \
                "(.*)\*(theta\['" + adsorbate_name + "'\])(\*{2}(\d)|)(\*(.*)|)"
            m = re.search(regex, term_expression)
            if m.group(6):
                site_part = m.group(1) + '*' + m.group(6)
            else:
                site_part = m.group(1)
            adsorbate_part = m.group(2) + m.group(3)
            #get derivation expression
            derivation_expression = \
                deriv_adsorbate_part(adsorbate_name, adsorbate_part) + \
                '*' + site_part + ' + ' + \
                deriv_site_part(site_name, site_part) + '*' + adsorbate_part
        elif site_name in sites_list:
            derivation_expression = deriv_site_part(site_name, term_expression)
        elif theta(adsorbate_name) in term_expression:
            derivation_expression = \
                deriv_adsorbate_part(adsorbate_name, term_expression)
        else:
            derivation_expression = '0'

        return derivation_expression

    def poly_adsorbate_derivation(self, adsorbate_name, poly_expression):
        """
        Expect a polynomial expression of dtheta_dt and an adsorbate_name,
        return a derivation expression about the adsorbate.
        """
        #split poly_expression
        poly_list = poly_expression.split()
        operators, terms = poly_list[3::2], poly_list[2::2]
        #generate derived term expressions
        derived_terms = [self.total_term_adsorbate_derivation(adsorbate_name,
                                                              term_expression)
                         for term_expression in terms]
        #combine 2 lists
        derived_poly_list = []
        for combination in zip(derived_terms, operators):
            derived_poly_list += list(combination)
        derived_poly_list.append(derived_terms[-1])

        return ' '.join(derived_poly_list)

    def analytical_jacobian(self, dtheta_dt_expressions, cvgs_tuple):
        "Get the jacobian matrix of the steady_state_function."
        #set theta, kf, kr, p, dtheta_dt
        #coverages(theta)
        theta = self.cvg_tuple2dict(cvgs_tuple)
        #rate constants(kf, kr)
        kf, kr = self.get_rate_constants()
        #pressure
        p = self.p
        #concentration
        c = self.c

        #generate Jacobian matrix
        m, n = len(dtheta_dt_expressions), len(cvgs_tuple)
        J = self._matrix(m, n)
        for i in xrange(m):
            poly_expression = dtheta_dt_expressions[i]
            for j in xrange(n):
                #get adsorbate_name
                adsorbate_name = self._owner.adsorbate_names[j]
                J[i, j] = \
                    eval(self.poly_adsorbate_derivation(adsorbate_name,
                                                        poly_expression))
        return J

    ######################################################
    ######                                          ######
    ###### calculate micro kinetic model with Sympy ######
    ######                                          ######
    ######################################################

    def get_elementary_dtheta_dt_sym(self, adsorbate_name,
                                     elementary_rxn_list):
        """
        Expect elementary_rxn_list and an adsorbate_name in it,
        return dtheta_dt symbols of the corresponding adsorbate
        in single elementary equation.
        """
        #species must be adsorbate
        if adsorbate_name not in self._owner.adsorbate_names:
            raise ValueError("'"+adsorbate_name+"' is not an adsorbate!")
        for state_list in elementary_rxn_list:
            for species_str in state_list:
                stoichiometry, site_name = self.split_species(species_str)
                if site_name == adsorbate_name:
                    break
            if site_name == adsorbate_name:
                #get state idx to get direction info
                state_idx = elementary_rxn_list.index(state_list)
                break
        #if adsorbate name not in elementary_rxn list, stop
        if site_name != adsorbate_name:
            return

        #get dtheta_dt sym according rate symbols
        rf_sym, rr_sym = self.get_single_rate_sym(elementary_rxn_list)
        if state_idx == 0:
            dtheta_dt_sym = -rf_sym + rr_sym
        else:
            dtheta_dt_sym = rf_sym - rr_sym

        return dtheta_dt_sym

    def get_adsorbate_dtheta_dt_sym(self, adsorbate_name):
        """
        Expect a adsorbate_name, and go through self.rxns_list,
        return dtheta_dt of the adsorbate.
        """
        #total_dtheta_dt_sym = sym.Symbol('0', is_real=True)
        total_dtheta_dt_sym = 0
        for elementary_rxn_list in self.rxns_list:
            dtheta_dt_sym = \
                self.get_elementary_dtheta_dt_sym(adsorbate_name,
                                                  elementary_rxn_list)
            if not dtheta_dt_sym:  # rxn equation do not contain the adsorbate
                continue
            total_dtheta_dt_sym = total_dtheta_dt_sym + dtheta_dt_sym

        return total_dtheta_dt_sym

    def get_dtheta_dt_syms(self, log_latex=False):
        "Go through adsorbate_names to get dtheta_dts list."
        dtheta_dt_syms = []
        for adsorbate_name in self._owner.adsorbate_names:
            dtheta_dt_sym = self.get_adsorbate_dtheta_dt_sym(adsorbate_name)
            dtheta_dt_syms.append(dtheta_dt_sym)

        dtheta_dt_syms = tuple(dtheta_dt_syms)
        self.dtheta_dt_syms = dtheta_dt_syms

        #latex strings
        dtheta_dt_latexs = self.get_latex_strs(part1=r'\frac{d\theta_{', part2=r'}}{dt} ',
                                               symbols=dtheta_dt_syms)
        self.dtheta_dt_latex = tuple(dtheta_dt_latexs)

        if log_latex:
            #log it
            self.logger.log_latex(self.dtheta_dt_latex)

        return dtheta_dt_syms

    def steady_state_function_by_sym(self, cvgs_tuple):
        """
        Recieve a coverages tuple containing coverages of adsorbates,
        return a tuple of dtheta_dts of corresponding adsorbates.
        """
        if not hasattr(self, 'dtheta_dt_syms'):
            self.get_dtheta_dt_syms()
        #get substitution dict
        subs_dict = self.get_subs_dict(cvgs_tuple=cvgs_tuple)
        #loop to get values of dtheta/dt
        dtheta_dts = []
        for dtheta_dt_sym in self.dtheta_dt_syms:
            dtheta_dt = self._mpf(dtheta_dt_sym.evalf(subs=subs_dict))
            dtheta_dts.append(dtheta_dt)

        return tuple(dtheta_dts)

    def analytical_jacobian_sym(self, dtheta_dt_syms):
        """
        Get the jacobian matrix symbol expressions of
        the dtheta/dt nonlinear equations.
        Return a jacobian matrix(in self._matrix form).
        """
        m = n = len(dtheta_dt_syms)
        sym_jacobian = self._matrix(m, n)
        for i in xrange(m):
            dthe_dt_sym = dtheta_dt_syms[i]
            for j in xrange(n):
                ads_name = self._owner.adsorbate_names[j]
                theta_sym = self.extract_symbol(ads_name, 'ads_cvg')
                sym_jacobian[i, j] = \
                    sym.Derivative(dthe_dt_sym, theta_sym).doit()

        return sym_jacobian

    def analytical_jacobian_by_sym(self, dtheta_dt_syms, cvgs_tuple):
        """
        Get the jacobian matrix of the dtheta/dt nonlinear equations.
        Return a jacobian matrix(in self._matrix form).
        """
        #get substitution dicts
        subs_dict = self.get_subs_dict(cvgs_tuple=cvgs_tuple)
        #get symbol jacobian matrix
        sym_jacobian = self.analytical_jacobian_sym(dtheta_dt_syms)
        #get numerial jacobian matrix
        num_jacobian = sym_jacobian.evalf(subs=subs_dict)
        #keep precision
#        m, n = num_jacobian.shape
#        for i in xrange(m):
#            for j in xrange(n):
#                num_jacobian[i, j] = self._mpf(num_jacobian[i, j])

        return num_jacobian  # may lose precision

    def get_rate_control_by_sym(self, RDS):
        """
        RDS: int, Rate Determining Step number.
        """
        #get quasi_quilibrium_solver instance
        _temp = __import__('quasi_equilibrium_solver',
                           globals(), locals(), ['QuasiEquilibriumSolver'])
        qe_solver = _temp.QuasiEquilibriumSolver(owner=self._owner)
        qe_solver.RDS = RDS  # set Rate Determining Step
        XTRCs = qe_solver.get_XTRCs()
        self.qe_solver = qe_solver

        return XTRCs

    ##########################################################
    ###### calculate micro kinetic model with Sympy END ######
    ##########################################################

    def get_residual(self, cvgs_tuple):
        "Return the minimum cvg rate wrt coverage."
        #constrain cvgs
        #cvgs_tuple = self.constrain_converage(cvgs_tuple)
        dtheta_dts = self.steady_state_function(cvgs_tuple)
        residual = max([abs(dtheta_dt) for dtheta_dt in dtheta_dts])
        return residual

    def get_steady_state_cvgs(self, c0, single_pt=False, full_info=False):
        """
        Expect an inital coverages tuple,
        use Newton Method to solving nonlinear equations,
        return steady state coverages, if converged.

        Parameters
        ----------
        single_pt : bool
            if True, no initial guess check.
        full_info : bool
            if True, show full info instead of brief info.
        """
        #Oh, intial coverage must have physical meaning!
        c0 = self.constrain_converage(c0)
        self.initial_guess = c0  # set self.initial_guess

        #start root finding algorithm
        f = self.steady_state_function
        f_resid = self.get_residual
        J = self.analytical_jacobian
        constraint = self.constrain_converage
        if hasattr(self, 'dtheta_dt_expressions'):
            f_expression = self.dtheta_dt_expressions
        else:
            f_expression = self.get_dtheta_dt_expressions()

        ############    Main Loop with changed initial guess   ##############

        cancel = False

        while not cancel:
            if f_resid(c0) <= self.tolerance and not single_pt:
                self._coverage = converged_cvgs = c0
                print 'Good initial guess.'
                #return c0
                break

            newton_iterator = NewtonRoot(
                f=f, J=J, x0=c0, constraint=constraint,
                norm=self._norm, mpfloat=self._mpf,
                matrix=self._matrix, Axb_solver=self._Axb_solver,
                dtheta_dt_expressions=f_expression
            )
            i = 0
            x = c0
            old_error = 1e99
            #lambda func used when log
            mpf2float = lambda x: [float(i) for i in x]
            if c0:
                #log initial guess
                self.logger.log(log_type='iteration',
                                event='inital_guess_success',
                                initial_cvgs=mpf2float(c0))

            #####    Sub LOOP for a c0    #####
            for x, error, fx in newton_iterator:
                i += 1  # counter for loop
                #if iterations is larger than 100, log every 10 steps, full info
                if i > 100 and i % 30 == 0:
                    self.logger.log(
                        log_type='iteration',
                        event='rootfinding_in-process',
                        n_iter=i, resid=float(f_resid(x)),
                        norm=float(error), fx=tuple(map(float, fx)),
                        x=tuple(map(float, x))
                    )
                elif i < 100:
                    if full_info:  # show full information on screen
                        self.logger.log(
                            log_type='iteration',
                            event='rootfinding_in-process',
                            n_iter=i, resid=float(f_resid(x)),
                            norm=float(error), fx=tuple(map(float, fx)),
                            x=tuple(map(float, x))
                        )
                    else:  # show brief information on screen
                        self.logger.log(
                            log_type='iteration',
                            event='rootfinding_in-process(brief)',
                            n_iter=i, resid=float(f_resid(x)),
                            norm=float(error)
                        )

                #less than tolerance
                if error < self.tolerance:
                    if f_resid(x) < self.tolerance:
                        #check whether there is minus value in x
                        for cvg in x:
                            if cvg < 0.0:
                                lt_zero = True  # less than 0
                                break
                            else:
                                lt_zero = False
                        # check END #
                        if not lt_zero:
                            self.logger.log(
                                log_type='iteration',
                                event='rootfinding_success',
                                n_iter=i,
                                resid=float(f_resid(x)),
                                root=mpf2float(x)
                            )
                            converged_cvgs = x
                            cancel = True
                            break
                        else:  # bad root, iteration continue...
                            self.logger.log(
                                log_type='iteration',
                                event='rootfinding_continue',
                                bad_root=mpf2float(x)
                            )
                    else:
                        #No constraint when iterations less than 6
#                        if x >= 5:
#                            x = self.constrain_converage(x)
                        error = f_resid(x)  # use residual as error

                #if convergence is slow when the norm is larger than 0.1
                elif (i > self.max_rootfinding_iterations or
                      abs(error - old_error) < 1e-4) and error > 1e-1:
                #elif i > self.max_rootfinding_iterations:
                    #print "root stable"
                    self.logger.log(
                        log_type='iteration',
                        event='rootfinding_break',
                        n_iter=i, resid=float(f_resid(x)),
                        root=mpf2float(x)
                    )
                    #jump out of loop for this c0
                    cancel = False
                    break

                #residual is almost stagnated
                elif abs(error - old_error) < self.stable_criterion:

                    self.logger.log(
                        log_type='iteration',
                        event='rootfinding_stable',
                        n_iter=i, resid=float(f_resid(x)),
                        root=mpf2float(x)
                    )
#                    print 'difference: ' + str(abs(error - old_error))
                    #jump out of loop for this c0
                    cancel = False
                    break

                old_error = error  # set old error to be compared in next loop
                self._coverage = x
                self._error = error

                #archive
                if i % 100 == 0:
                    self.logger.archive_data('iter_coverage', x)
                    self.logger.archive_data('iter_error', error)

            #####    Sub loop for a c0 END    #####
            #change the initial guess(c0)
            if not cancel:
                #get a new initial guess coverage
                c0 = self.modify_init_guess(x, fx)

        ##############    main loop end   #################

        if converged_cvgs:
            self._coverage = converged_cvgs
            #archive converged root and error
            self.logger.archive_data('steady_state_coverage',
                                     converged_cvgs)
            if error:
                self.logger.archive_data('steady_state_error', error)
            self.good_guess = c0
            #archive initial guess
            self.logger.archive_data('initial_guess', c0)
            #log steady state coverages
            self.logger.log_sscvg(converged_cvgs, self._owner.adsorbate_names)
            return converged_cvgs

    def modify_init_guess(self, c0, dtheta_dts):
        "Return a new initial guess according to dthe_dts."
#        max_dtheta_dt = np.max(np.abs(dtheta_dts))
        base_coefficient = self.initial_guess_scale_factor
        coefficients = []
        for dtheta_dt in np.abs(dtheta_dts):
            if abs(dtheta_dt) >= self.tolerance:
                #coefficients.append(dtheta_dt/max_dtheta_dt*base_coefficient)
                coefficients.append(base_coefficient)
            else:
                coefficients.append(1.0)
        print 'coeff: ' + str(coefficients)
        #if coeffs are all 1.0, break!
        #add later...

        #create a diagnol matrix
        c0_diag = np.matrix(np.diag(c0))
        #convert coeffients to column vector
        coefficients = np.matrix(coefficients).reshape(-1, 1)
        new_c0 = (c0_diag*coefficients).reshape(1, -1)
        new_c0 = tuple(new_c0.tolist()[0])
        #add to log
        #lambda func used when log
        mpf2float = lambda x: [float(i) for i in x]
        self.logger.log(
            log_type='iteration',
            event='modify_initial_guess_success',
            old_init_cvg=mpf2float(c0),
            new_init_cvg=mpf2float(new_c0),
        )

        #return self.constrain_converage(new_c0)
        return new_c0

    def modify_init_guess_new(self, c0, dtheta_dts):
        "Return a new initial guess according to dthe_dts."
        max_dtheta_dt = np.max(np.abs(dtheta_dts))
        base_coefficient = self.initial_guess_scale_factor
        coefficients = []
        for idx, dtheta_dt in enumerate(np.abs(dtheta_dts)):
            if abs(dtheta_dt) >= self.tolerance:
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
        #lambda func used when log
        mpf2float = lambda x: [float(i) for i in x]
        self.logger.log(
            log_type='iteration',
            event='modify_initial_guess_success',
            old_init_cvg=mpf2float(c0),
            new_init_cvg=mpf2float(new_c0),
        )

        return new_c0

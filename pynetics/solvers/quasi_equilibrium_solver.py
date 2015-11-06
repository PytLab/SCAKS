# -*- coding: utf-8 -*-
import copy
import logging

import sympy as sym

from solver_base import SolverBase


class QuasiEquilibriumSolver(SolverBase):
    def __init__(self, owner):
        SolverBase.__init__(self, owner)
        self.logger = logging.getLogger('model.solvers.QuasiEquilibriumSolver')

        #set default parameter dict
        defaults = dict(
            RDS=2,  # rate determining step number
            )
        defaults = self.update_defaults(defaults)
        self.__dict__.update(defaults)

        #update quasi_equilibrium_solver's own logger template
        self.logger_template_dict = {}
        self.logger._templates_dict.update(self.logger_template_dict)

        #species names that has been represented by theta*
        self.represented_species = []

        #get data symbols
        if not self.has_symbols:
            self.get_data_symbols()

    def get_theta_f_sym(self):  # 请找个安静的地方读这段代码
        """
        Return a solved analytical expression of theta_f and
        theta_f symbol object.
        """
        #refresh self.represented_species
        self.represented_species = []
        #operate copy later
        rxns_list_copy = copy.copy(self.rxns_list)
        Ks_list_copy = list(copy.copy(self.K_sym))
        #remove rate determinating step
        RDS_rxn_list = rxns_list_copy[self.RDS]
        rxns_list_copy.remove(RDS_rxn_list)
        #remove rate determinating step's K
        RDS_K = Ks_list_copy[self.RDS]
        Ks_list_copy.remove(RDS_K)
        #create a free site theta symbol used in all rxns
        theta_f = sym.Symbol('theta_f', positive=True, real=True)
        syms_sum = 0  # sum expression of all adsorbates' thetas

        #set equivalent dict
        #A complete substitution dict for each adsorbate symbols
        self.eq_dict = {}

        rxns_num = len(self.rxns_list)

        loop_counter = 0
        while rxns_list_copy:
#            print rxns_list_copy
            origin_num = len(rxns_list_copy)  # number of rxns

            for K_sym, rxn_list in zip(Ks_list_copy, rxns_list_copy):
#                print rxns_list_copy
#                print Ks_list_copy
                loop_counter += 1
                #get adsorbate name that will be represented
                target_adsorbate = self.check_repr(rxn_list)

                if target_adsorbate and target_adsorbate != 'all_represented':
                    #get target adsorbate theta symbol
                    theta_target = self.extract_symbol(target_adsorbate, 'ads_cvg')
                    #represented by theta_f
#                    print "target_adsorbate: %s" % target_adsorbate
                    theta_target_subs = self.represent(rxn_list, target_adsorbate,
                                                       theta_f, K_sym)
#                    print theta_target_subs
                    theta_target_subs = theta_target_subs.subs(self.eq_dict)
#                    print theta_target_subs

                    #add it to self.eq_dict
                    if theta_target in self.eq_dict:
                        self.eq_dict[theta_target] = theta_target_subs
                    else:
                        self.eq_dict.setdefault(theta_target, theta_target_subs)
                    #add theta to sym_sum
#                    print "added theta %s: %s" % (str(target_adsorbate), str(theta_target_subs))
                    syms_sum += theta_target_subs
                    #add this good species name to self.represented_species
                    self.represented_species.append(target_adsorbate)

                    #Ok, come on bros!
                    bros, bros_syms = \
                        self.get_related_adsorbates_of_adsorbate(target_adsorbate)

                    if bros and bros_syms:
                        if not hasattr(self, 'related_theta_subs_dict'):
                            self.get_related_theta_subs_dict()

                        #add bros to self.eq_dict, syms_sum, self.represented_species
                        for rel_ads, rel_ads_sym in zip(bros, bros_syms):
                            key = rel_ads_sym
                            value = self.related_theta_subs_dict[key].subs(self.eq_dict)
                            #add this bro to self.eq_dict
                            self.eq_dict[key] = value
                            #add this bro to syms_sum
                            rel_ads_sym = rel_ads_sym.subs(self.eq_dict)  # bro symbol needs substitution
#                            print "added theta %s: %s" % (str(rel_ads), str(rel_ads_sym))
                            syms_sum += rel_ads_sym
                            #add this bro to self.represented_species
                            self.represented_species.append(rel_ads)

                    rxns_list_copy.remove(rxn_list)
                    Ks_list_copy.remove(K_sym)
                elif target_adsorbate and target_adsorbate == 'all_represented':
                    #just remove it
                    rxns_list_copy.remove(rxn_list)
                    Ks_list_copy.remove(K_sym)
                else:
                    #move the rxn_list to the end of rxns_list_copy
                    rxns_list_copy.remove(rxn_list)  # remove it
                    rxns_list_copy.append(rxn_list)  # insert to the end
                    #move the K_sym to the end of Ks_list_copy
                    Ks_list_copy.remove(K_sym)
                    Ks_list_copy.append(K_sym)

#                print rxns_list_copy  # see what left

            remaining_num = len(rxns_list_copy)  # number of rxn remaining in list
#            print "remain: %d" % remaining_num
#            print "loop count: %d" % loop_counter

            if remaining_num == origin_num and loop_counter > rxns_num:
#                print "In merging part..."
#                print rxns_list_copy
                #insert K for merged rxn list to head of Ks_list_copy
                merged_K = self.get_merged_K(rxns_list_copy)
                Ks_list_copy.insert(0, merged_K)
                #merge all remaining elementary_rxn lists
                merged_rxn_list = \
                    self._owner.parser.merge_elementary_rxn_list(*rxns_list_copy)
                #insert new merged list to the head of rxns_list_copy
                rxns_list_copy.insert(0, merged_rxn_list)

        #get theta_f expression
#        print syms_sum
        normalization_expr = syms_sum + theta_f - 1
        theta_f_expr = sym.solve(normalization_expr, theta_f, check=0)[0]
#        print theta_f_expr
        return theta_f, theta_f_expr

    def get_related_adsorbates_of_adsorbate(self, adsorbate_name):
        "Return a list containing names of related adsorbate of adsorbate in parameter."
        if not hasattr(self._owner, 'related_adsorbate_names'):
            self._owner.parser.get_related_adsorbates()

        bros, bros_syms = [], []
        for rel_ads_tup in self._owner.related_adsorbate_names:
            if adsorbate_name in rel_ads_tup:
                bros.extend(list(rel_ads_tup))
                break
        if bros:
            bros.remove(adsorbate_name)
            #get corresponding theta symbols
            for ads in bros:
                ads_sym = self.extract_symbol(ads, 'ads_cvg')
                bros_syms.append(ads_sym)

        return bros, bros_syms

    def get_simplified_tof_sym(self, theta_f, theta_f_expr):
        "No substitution second time to get a simplifed expression with K."
        #get theta_f, theta_f_expr
        theta_f, theta_f_expr = self.get_theta_f_sym()
        subs_dict = dict({theta_f: theta_f_expr}, **self.eq_dict)
        #get net rate of rate determinating step
        if not hasattr(self, 'net_rate_syms'):
            self.get_net_rate_syms()
        tof_sym = self.net_rate_syms[self.RDS]
        #substitute thetas of adsorbates
        tof_sym = tof_sym.subs(subs_dict)

        return tof_sym

    def get_tof_sym(self):
        "Do substitution twice to get complete expanded expression."
        #get theta_f, theta_f_expr
        theta_f, theta_f_expr = self.get_theta_f_sym()
        #get complete equivalent dict
        complete_eq_dict = self.get_complete_eq_dict(theta_f, theta_f_expr)
#        print complete_eq_dict
        #get net rate of rate determinating step
        if not hasattr(self, 'net_rate_syms'):
            self.get_net_rate_syms()
        tof_sym = self.net_rate_syms[self.RDS]
        #substitute thetas of adsorbates
        complete_tof_sym = tof_sym.subs(complete_eq_dict)
        #substitute again to subs Ks in complete_eq_dict itself!
        #e.g. complete_eq_dict[theta_H_s]
        complete_tof_sym = complete_tof_sym.subs(complete_eq_dict)

        return complete_tof_sym

    def get_tof(self):
        #get substitution dict
        subs_dict = self.get_subs_dict()
        tof_sym = self.get_tof_sym()
        tof = tof_sym.evalf(subs=subs_dict)

        return tof

    def get_merged_K(self, rxns_list):
        "Merged K is the multiplication of K for corresponding rxns."
        merged_K = 1
        for rxn_list in rxns_list:
            #get index
            idx = self.rxns_list.index(rxn_list)  # could use rxns_list.index(rxn_list) here
                                                  # if the first merging dosen't work
                                                  # an exception will be raised here
            #get corresponding K
            K = self.K_sym[idx]
            merged_K *= K
        return merged_K

    def get_complete_eq_dict(self, theta_f, theta_f_expr):
        "substitute theta_f and K, to get complete equivalent dict."
        #check number of elements in eq_dict
        ads_num = len(self._owner.adsorbate_names)
        if len(self.eq_dict) != ads_num:
            raise ValueError('eq_dict is illegal.')
        for ads_sym in self.eq_dict:
            self.eq_dict[ads_sym] = \
                self.eq_dict[ads_sym].subs({theta_f: theta_f_expr})
        #get equilibrium constant subs dict
        if not hasattr(self, 'K_expr_syms'):
            self.get_K_syms()
        K_subs_dict = {}
        for i, K_sym in enumerate(self.K_sym):
            K_subs_dict.setdefault(K_sym, self.K_expr_syms[i])

        #merge two dicts
        self.eq_dict = dict(self.eq_dict, **K_subs_dict)

        return self.eq_dict  # Note: K still in it!

    def represent(self, rxn_list, target_adsorbate, theta_f, K):
        """
        Expect a rxn_list which can be represented by theta_f,
        return the symbol expression of theta_target_adsorbate.

        Parameters
        ----------
        rxn_list : list of states
            e.g. [['*_s', 'HCOOH_g'], ['HCOOH_s']]
        target_adsorbate : str
            adsorbate_name which will be represented
        theta_f : sympy.core.symbol.Symbol
            coverage of free sites

        Example
        -------
        >>> m.solver.represent([['O2_s', 'NO_s'], ['*_s', 'ONOO_s']], 'NO_s', theta_f, K)
        >>> theta_O2_s**(-1.0)*theta_ONOO_s**1.0*(theta_f/K)**1.0

        """
        left_syms, right_syms = [], []  # syms to be multipled later

        #go thtough rxn_list's head and tail
        ends_list = [rxn_list[0], rxn_list[-1]]
        for state_idx, state_list in enumerate(ends_list):
            #go through sp_list to locate theta
            for sp_str in state_list:
                stoichiometry, species_name = self.split_species(sp_str)
                #free site
                if '*' in species_name:
                    sym_term = theta_f**stoichiometry  # add exponential
                #represented adsorbate
                elif species_name in self._owner.adsorbate_names:
                    #extract symbol of the adsorbate
                    theta_sym = self.extract_symbol(species_name, 'ads_cvg')
                    sym_term = theta_sym**stoichiometry  # add exponential
                #gas pressure
                elif species_name in self._owner.gas_names:
                    p_sym = self.extract_symbol(species_name, 'pressure')
                    sym_term = p_sym**stoichiometry  # add exponential
                #liquid concentration
                elif species_name in self._owner.liquid_names:
                    c_sym = self.extract_symbol(species_name, 'concentration')
                    sym_term = c_sym**stoichiometry
                #classify the symbol term to left or right list
                if state_idx == 0:
                    left_syms.append(sym_term)
                else:
                    right_syms.append(sym_term)

        def get_multi_sym_list(syms_list):
            "product all elements in symbols list, left or right list."
            multi_sym = 1
            for symbol in syms_list:
                multi_sym *= symbol
            return multi_sym

        #get equation to be solved
        L = get_multi_sym_list(left_syms)  # left multipled symbols
        R = get_multi_sym_list(right_syms)  # right multipled symbols
        equation = R/L - K  # R/L - K = 0
        #get target theta symbol
        theta_t = self.extract_symbol(target_adsorbate, 'ads_cvg')
        #do substitution
        #get substitution dict
        if not hasattr(self, 'related_theta_subs_dict'):
            self.get_related_theta_subs_dict()
        related_theta_subs_dict = self.related_theta_subs_dict
        #substitute the other adsorbates theta with theta_t in equation
        equation = equation.subs(related_theta_subs_dict)
        #solve the equation to get theta_t
        represented_theta_t = sym.solve(equation, theta_t)  # list
        if len(represented_theta_t) == 1:
            ans = represented_theta_t[0]
        else:
            raise ValueError('No unique solution! Solutions: %s' %
                             str(represented_theta_t))

        return ans  # theta_t actually

    def get_related_theta_subs_dict(self):
        """
        Use the first adsorbate theta symbol to represented the other theta symbols.
        """
        #add related adsorbates theta symbol substitution to eq_dict
        if not hasattr(self._owner, 'related_adsorbates'):
            self._owner.parser.get_related_adsorbates()
        #related_adsorbates is like [{}, {'H_s': 1, 'OH_s': 1}]
        related_theta_subs_dict = {}
        for rel_ads_dict in self._owner.related_adsorbates:
            if rel_ads_dict:
                related_adsorbate_names = sorted(rel_ads_dict.keys())
                # the first adsorbate is the pivot adsorbate
                pivot_adsorbate = related_adsorbate_names[0]
                denominator = rel_ads_dict[pivot_adsorbate]
                pivot_theta_sym = self.extract_symbol(pivot_adsorbate, 'ads_cvg')
                for adsorbate_name in related_adsorbate_names[1:]:
                    theta_sym = self.extract_symbol(adsorbate_name, 'ads_cvg')
                    key = theta_sym
                    numerator = rel_ads_dict[adsorbate_name]
                    ratio = numerator/float(denominator)
                    value = ratio*pivot_theta_sym
                    related_theta_subs_dict.setdefault(key, value)
        self.related_theta_subs_dict = related_theta_subs_dict

        return related_theta_subs_dict

    def check_repr(self, rxn_list):
        """
        Expect a rxn_list, e.g. [['*_s', 'HCOOH_g'], ['HCOOH_s']].
        Check if there is only one adsorbate cvg
        that can be represented by cvg of free site theta_f.

        e.g. [['*_s', 'HCOOH_g'], ['HCOOH_s']] is good, return 'HCOOH_s'
        [['HCOOH_s', '*_s'], ['H-COOH_s', '*_s'], ['COOH_s', 'H_s']] not,
        return None.
        """
        #merge IS and FS sp_list
        merged_list = rxn_list[0] + rxn_list[-1]
        #initial number of theta that hasn't been represented by theta_*
        free_num = 0
        #adsorbate name that will be archived in loop
        archived_adsorbates = []
        for sp_str in merged_list:
            stoichiometry, species_name = self.split_species(sp_str)
            if (species_name in self._owner.adsorbate_names and
                    species_name not in self.represented_species):
                free_num += 1
                archived_adsorbates.append(species_name)

        archived_adsorbates = tuple(sorted(archived_adsorbates))
#        print "archived_adsorbates: %s" % str(archived_adsorbates)

        if not hasattr(self._owner, 'related_adsorbate_names'):
            self._owner.parser.get_related_adsorbates()

        if free_num == 1:
            return archived_adsorbates[0]
        elif free_num == 0:
            return 'all_represented'
        elif archived_adsorbates in self._owner.related_adsorbate_names:  # if the archived adsorbates are related
            return archived_adsorbates[0]                                 # return them all in a tuple
        else:
            return

    def get_XTRC(self, intermediate_name):
        r = self.get_tof_sym()  # tof symbol expression
        G = self.extract_symbol(intermediate_name, 'free_energy')  # free energy symbol
        k_B, T = self.k_B_sym, self.T_sym
        XTRC = -k_B*T/r*(sym.Derivative(r, G).doit())

        subs_dict = self.get_subs_dict()
        XTRC_value = XTRC.evalf(subs=subs_dict)

        return XTRC_value

    def get_XTRCs(self):
        sp_list = self._owner.adsorbate_names + \
            self._owner.transition_state_names

        XTRCs = [self.get_XTRC(sp) for sp in sp_list]

        return XTRCs

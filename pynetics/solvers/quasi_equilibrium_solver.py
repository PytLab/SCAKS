from solver_base import *
import sympy as sym


class QuasiEquilibriumSolver(SolverBase):
    def __init__(self, owner):
        SolverBase.__init__(self, owner)

        #set default parameter dict
        defaults = dict(
            RDS=3,  # rate determining step number
            )
        defaults = self.update_defaults(defaults)
        self.__dict__.update(defaults)

        #update quasi_equilibrium_solver's own logger template
        self.logger_template_dict = {}
        self.logger._templates_dict.update(self.logger_template_dict)

        self.represented_thetas = []  # thetas that has been represented by theta*

    def get_tof_sym(self):
        waiting_room = []
        for rxn_idx, rxn_list in enumerate(self.rxns_list):
            #get adsorbate name that will be represented
            target_adsorbate = self.check_repr(rxn_list)
            #get target adsorbate theta symbol
            theta_target = self.extract_symbol(target_adsorbate, 'ads_cvg')
            #create a free site theta symbol used in all rxns
            theta_f = sym.Symbol('theta_f', is_positive=True, is_real=True)

            if target_adsorbate:
                pass

            else:  # throw it to waiting room
                waiting_room.append(rxn_list)

    def represent(self, rxn_list, target_adsorbate, theta_f):
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
        """
        #get equilibrium constant symbol
        rxn_idx = self.rxns_list.index(rxn_list)
        K = self.K_sym[rxn_idx]

        left_syms, right_syms = [], []  # syms to be multipled later

        #go thtough rxn_list's head and tail
        ends_list = [rxn_list[0], rxn_list[-1]]
        for state_idx, state_list in enumerate(ends_list):
            #go through sp_list to locate theta
            for sp_str in state_list:
                stoichiometry, species_name = self.split_species(sp_str)
                #free site
                if '*' in species_name:
                    # get location of theta_f
                    if state_idx == 0:
                        theta_f_loc = 'left'
                    else:
                        theta_f_loc = 'right'
                    #add exponential
                    theta_f_term = theta_f**stoichiometry
                #target species theta
                elif species_name == target_adsorbate:
                    # get location of theta_target_ads
                    if state_idx == 0:
                        theta_t_loc = 'left'
                    else:
                        theta_t_loc = 'right'
                    #get theta_t term exponential
                    theta_t_exp = stoichiometry
                #represented adsorbate
                elif species_name in self._owner.adsorbate_names:
                    #extract symbol of the adsorbate
                    theta_sym = self.extract_symbol(species_name, 'ads_cvg')
                    #add exponential
                    theta_sym_term = theta_sym**stoichiometry
                    #get location of theta represented already
                    if state_idx == 0:
                        left_syms.append(theta_sym_term)
                    else:
                        right_syms.append(theta_sym_term)
                #gas pressure
                elif species_name in self._owner.gas_names:
                    p_sym = self.extract_symbol(species_name, 'pressure')
                    #add exponential
                    p_sym_term = p_sym**stoichiometry
                    if state_idx == 0:
                        left_syms.append(p_sym_term)
                    else:
                        right_syms.append(p_sym_term)

        #use theta_f to represent theta_target_adsorbate
        '''
        --------------------------------------------------------------------------
        |theta_* location | theta_t_loc | theta_target_adsorbate expression      |
        --------------------------------------------------------------------------
        |    left         |    left     |            R/(K*L*theta_f_term)        |
        --------------------------------------------------------------------------
        |    left         |    right    |            K*L*theta_f_term/R          |
        --------------------------------------------------------------------------
        |    right        |    left     |            R*theta_f_term/(L*K)        |
        --------------------------------------------------------------------------
        |    right        |    right    |            K*L/(R*theta_f_term)        |
        --------------------------------------------------------------------------
        '''
        def get_multi_sym_list(syms_list):
            multi_sym = 1
            for symbol in syms_list:
                multi_sym *= symbol
            return multi_sym

        L = get_multi_sym_list(left_syms)  # left multipled symbols
        R = get_multi_sym_list(right_syms)  # right multipled symbols
        total_exp = 1.0/theta_t_exp

        if theta_f_loc == 'left' and theta_t_loc == 'left':
            theta_t = (R/(K*L*theta_f_term))**total_exp
        elif theta_f_loc == 'left' and theta_t_loc == 'right':
            theta_t = (K*L*theta_f_term/R)**total_exp
        elif theta_f_loc == 'right' and theta_t_loc == 'left':
            theta_t = (R*theta_f_term/(L*K))**total_exp
        elif theta_f_loc == 'right' and theta_t_loc == 'right':
            theta_t = (K*L/(R*theta_f_term))**total_exp

        return theta_t

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
        archived_ads = ''
        for sp_str in merged_list:
            stoichiometry, species_name = self.split_species(sp_str)
            if (species_name in self._owner.adsorbate_names and
                    species_name not in self.represented_thetas):
                free_num += 1
                archived_ads = species_name

        if free_num == 1:
            return archived_ads
        else:
            return

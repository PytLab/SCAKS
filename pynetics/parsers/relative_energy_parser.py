import os

import sympy as sym

from parser_base import *


class RelativeEnergyParser(ParserBase):
    def __init__(self, owner):
        '''
        A class to parse relative energies
        covert them to generalize formation energies.
        '''
        ParserBase.__init__(self, owner)

        #intialize generalized formation energy dict
        self.G_dict = {}

        #parse in data file
        if os.path.exists('rel_energy.py'):
            globs, locs = {}, {}
            execfile('rel_energy.py', globs, locs)
            #set variables in data file as attr of parser
            for key in locs:
                setattr(self, key, locs[key])
        else:
            raise ValueError("No rel_energy.py in current path.")

    def chk_data_validity(self):
        #check data validity
        if (not self.Ea or
                not (len(self.Ea) == len(self._owner.elementary_rxns_list))):
            raise ValueError("No Ea or invalid shape of Ea.")
        elif (not self.G0 or
                not (len(self.G0) == len(self._owner.elementary_rxns_list))):
            raise ValueError("No G0 or invalid shape of G0.")

    def get_energy_symbols(self):
        """
        Get energy symbols in order of site_names +
        gas_names + adsorbate_names + transition_state_names.
        """
        self.all_sp = self._owner.site_names + self._owner.gas_names + \
            self._owner.adsorbate_names + self._owner.transition_state_names

        energy_symbols = []
        for sp_name in self.all_sp:
            if sp_name in self._owner.ref_species:
                #set reference energies as 0
                self.G_dict.setdefault(sp_name, 0.0)
            symbol = sym.Symbol('G_' + sp_name, real=True, positive=True)
            energy_symbols.append(symbol)

        self.energy_symbols = energy_symbols

        return energy_symbols

    def str_sym(self, obj):
        "Return corresponding string of the symbols and symbol of the string."
        if not hasattr(self, 'energy_symbols'):
            self.get_energy_symbols()

        if obj in self.energy_symbols:
            idx = self.energy_symbols.index(obj)
            return self.all_sp[idx]
        elif obj in self.all_sp:
            idx = self.all_sp.index(obj)
            return self.energy_symbols[idx]
        else:
            raise ValueError("%s not in energy_symbols or all_sp." % obj)

    def get_poly_term(self, state_list):
        term = 0.0
        for sp_name in state_list:
            if '*' in sp_name:
                star, sp_name = sp_name.split('_')
            if sp_name in self._owner.ref_species:
                symbol = 0.0  # float number actually
            else:
                symbol = self.str_sym(sp_name)
            term += symbol

        return term

    def get_single_polys(self, elementary_rxn_list):
        """
        Expect a elementary_rxn_list,
        e.g. [['*_s', 'NO_g'], ['NO-_s'], ['NO_s']]
        return corresponding polynomial wrt Ea and G0.
        e.g. [G_NO-_s - 0.655, G_NO_s + 0.455]
        """
        idx = self._owner.elementary_rxns_list.index(elementary_rxn_list)
        Ea, G0 = self.Ea[idx], self.G0[idx]

        polys = []  # Ea equation and G0 equation
        if Ea != 0 and len(elementary_rxn_list) != 2:  # has barrier
            #get Ea equation
            is_list, ts_list, fs_list = elementary_rxn_list
            is_term = self.get_poly_term(is_list)
            ts_term = self.get_poly_term(ts_list)
            ts_poly = ts_term - is_term - Ea
            polys.append(ts_poly)
        #get G0 equation
        is_list, fs_list = elementary_rxn_list[0], elementary_rxn_list[-1]
        is_term = self.get_poly_term(is_list)
        fs_term = self.get_poly_term(fs_list)
        fs_poly = fs_term - is_term - G0
        polys.append(fs_poly)

        return polys

    def solve_equations(self):
        """
        Solve a system of equations to get values of generalized formation energy.
        """
        #get polynomials list
        equations_list = []
        for elementary_rxn_list in self._owner.elementary_rxns_list:
            polys_list = self.get_single_polys(elementary_rxn_list)
            equations_list.extend(polys_list)

        return equations_list

    ####### Use matrix to get generalized formation energy ########
    def get_unknown_species(self):
        "Get species whose free energy is unknown."
        all_sp = self._owner.site_names + self._owner.gas_names + \
            self._owner.adsorbate_names + self._owner.transition_state_names
        all_sp = list(all_sp)

        for known_sp in self._owner.ref_species:
            all_sp.remove(known_sp)
            self.G_dict.setdefault(known_sp, 0.0)

        self.unknowns = all_sp

        return all_sp

    def list2dict(self, state_list):
        """
        Expect a state list, e.g. ['*_s', 'NO_g']
        return a corresponding dict, e.g. {'*_s': 1, 'NO_g': 1}.
        """
        state_dict = {}
        for sp_str in state_list:
            stoichiometry, species_name = self.split_species(sp_str)
            state_dict.setdefault(species_name, stoichiometry)

        return state_dict

    def get_unknown_coeff_vector(self, elementary_rxn_list):
        """
        Expect a elementary_rxn_list,
        e.g. [['O2_s', 'NO_g'], ['*_s', 'ON-OO_s'], ['*_s', 'ONOO_s']]
        return coefficient vectors, Ea, G0.
        e.g. ([[0, 0, -1, 0, 0, 0, 0, 0, 1], [0, 0, -1, 1, 0, 0, 0, 0, 0]], 0.655, -0.455)
        """
        idx = self._owner.elementary_rxns_list.index(elementary_rxn_list)
        Ea, G0 = self.Ea[idx], self.G0[idx]

        if not hasattr(self, 'unknowns'):
            self.get_unknown_species()

        coeff_vects = []
        if Ea != 0 and len(elementary_rxn_list) != 2:  # has barrier
            #get ts coefficient vector
            is_list, ts_list = elementary_rxn_list[0], elementary_rxn_list[1]
            is_dict, ts_dict = self.list2dict(is_list), self.list2dict(ts_list)
            coeff_vect = []
            for unknown in self.unknowns:
                if unknown in is_dict:
                    coeff = -is_dict[unknown]
                elif unknown in ts_dict:
                    coeff = ts_dict[unknown]
                else:
                    coeff = 0
                coeff_vect.append(coeff)
            coeff_vects.append(coeff_vect)

        #coefficient vector for G0
        is_list, fs_list = elementary_rxn_list[0], elementary_rxn_list[-1]
        is_dict, fs_dict = self.list2dict(is_list), self.list2dict(fs_list)
        coeff_vect = []
        for unknown in self.unknowns:
            if unknown in is_dict:
                coeff = -is_dict[unknown]
            elif unknown in fs_dict:
                coeff = fs_dict[unknown]
            else:
                coeff = 0
            coeff_vect.append(coeff)
        coeff_vects.append(coeff_vect)

        if Ea:
            return coeff_vects, [Ea, G0]
        else:
            return coeff_vects, [G0]

    def parse_data(self):  # correspond with parse_data() in csv_parser.py
        """
        Solve Axb equation to get value of generalized free energies.
        """
        A, b = [], []
        for rxn_list in self._owner.elementary_rxns_list:
            coeff_vects, value = self.get_unknown_coeff_vector(rxn_list)
            A.extend(coeff_vects)
            b.extend(value)

        A, b = np.matrix(A), np.matrix(b).reshape(-1, 1)
#        print A
#        self.A = A
#        print b
#        self.b = b

        x = A.I*b  # values for unknowns
        #convert column vector to list
        x = x.reshape(1, -1).tolist()[0]

        #put values to G_dict
        for sp_name, G in zip(self.unknowns, x):
            self.G_dict.setdefault(sp_name, G)

        #put generalized formation energy into species_definition
        for sp_name in self.G_dict:
            sp_dict = self._owner.species_definitions
            sp_dict[sp_name].setdefault('formation_energy', self.G_dict[sp_name])

        setattr(self._owner, 'hasdata', True)

        return

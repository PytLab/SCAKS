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

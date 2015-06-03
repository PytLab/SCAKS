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
        if obj in self.energy_symbols:
            idx = self.energy_symbols.index(obj)
            return self.all_sp[idx]
        elif obj in self.all_sp:
            idx = self.all_sp.index(obj)
            return self.energy_symbols[idx]
        else:
            raise ValueError("%s not in energy_symbols or all_sp." % obj)

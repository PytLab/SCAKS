import sympy as sym

from parser_base import *


class RelativeEnergyParser(ParserBase):
    def __init__(self, owner):
        '''
        A class to parse relative energies
        covert them to generalize formation energies.
        '''
        ParserBase.__init__(self, owner)

    def get_energy_symbols(self):
        """
        Get energy symbols in order of site_names +
        gas_names + adsorbate_names + transition_state_names.
        """
        all_sp = self._owner.site_names + self._owner.gas_names + \
            self._owner.adsorbate_names + self._owner.transition_state_names

        energy_symbols = []
        for sp_name in all_sp:
            if sp_name in self._owner.ref_species:
                continue
            else:
                symbol = sym.Symbol('G_' + sp_name, real=True,
                                    positive=True)
                energy_symbols.append(symbol)

        self.energy_symbols = energy_symbols

        return energy_symbols

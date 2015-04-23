from corrector_base import *
from parameter_data import *
import numpy as np
import mpmath as mp


class ThermodynamicCorrector(CorrectorBase):
    def __init__(self, owner):
        CorrectorBase.__init__(self, owner)

        defaults = dict(
            _kJmol2eV=0.01036427,
            )
        self.__dict__.update(defaults)

        self._zpe_dict = {}
        self._enthalpy_dict = {}
        self._entropy_dict = {}

    def shomate_gas(self):
        """
        Correct free energies of gases with shomate_gas,
        then update E dict of solver.
        """
        gas_names = self._owner.gas_names
        temperature = self._owner.temperature
        temperature_ref = 298.15
        #shomate_params = shomate_params

        def H(T, params):
            A, B, C, D, E, F, G, H = params
            t = T/1000.0
            H = A*t + (B/2.0)*t**2 + (C/3.0)*t**3 + (D/4.0)*t**4 - E/t + F - H
            #kJ/mol
            return H

        def S(T, params):
            A, B, C, D, E, F, G, H = params
            t = T/1000.0
            S = A*np.log(t) + B*t + (C/2.0)*t**2 + (D/3.0)*t**3 - E/(2.0*t**2)\
                + G  # J/mol*K
            return S

        def Cp(T, params):
            A, B, C, D, E, F, G, H = params
            t = T/1000.0
            Cp = A + B*t + C*t**2 + D*t**3 + E/(t**2)
            return Cp

        thermo_dict = {}
        for gas in gas_names:
            for key in shomate_params.keys():
                gas_key, T_range = key.split(':')
                T_min, T_max = [float(t) for t in T_range.split('-')]
                if (gas == gas_key
                        and temperature >= T_min
                        and temperature <= T_max):
                    params = shomate_params[key]
                    Cp_ref = Cp(temperature_ref, params)
                    dH = H(temperature, params) - H(temperature_ref, params)
                    #deltaH(298-T) = shomate(T) - shomate(298)
                    dS = S(temperature, params)
                    dH = (temperature_ref*Cp_ref/1000.0 + dH)*(self._kJmol2eV)  # eV
                    #dH = 298*Cp(298) + dH(298-T)
                    dS = dS*(self._kJmol2eV/1e3)  # eV/K
                    #ZPE = sum(self.frequencies[gas])/2.0
                    #free_energy = ZPE +  dH - temperature*dS
                    free_energy = dH - temperature*dS
                    #self._zpe_dict[gas] = ZPE
                    self._enthalpy_dict[gas] = dH
                    self._entropy_dict[gas] = dS
                    thermo_dict[gas] = mp.mpf(free_energy)
                elif temperature < T_min and T_min < 300:
                    params = shomate_params[key]
                    Cp_ref = Cp(T_min, params)
                    dS = S(T_min, params)
                    dH = (temperature*Cp_ref/1000.0)*(self._kJmol2eV)  # eV
                    dS = dS*(self._kJmol2eV/1e3)  # eV/K
                    #ZPE = sum(self.frequencies[gas])/2.0
                    free_energy = dH - temperature*dS
                    #self._zpe_dict[gas] = ZPE
                    self._enthalpy_dict[gas] = dH
                    self._entropy_dict[gas] = dS
                    thermo_dict[gas] = mp.mpf(free_energy)
        for key in gas_names:
            not_there = []
            if key not in thermo_dict:
                not_there.append(key)
            if not_there:
                raise ValueError('No Shomate parameters specified for '+' '.join(not_there))

        return thermo_dict

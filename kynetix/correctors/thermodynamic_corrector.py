import logging
from math import pi, exp, log

import numpy as np
import mpmath as mp

from kynetix.correctors.corrector_base import *
from kynetix.database.thermo_data import *
from kynetix.errors.error import *


# Constants.
_kJmol2eV = 0.01036427


class ThermodynamicCorrector(CorrectorBase):
    def __init__(self, owner):
        super(ThermodynamicCorrector, self).__init__(owner)

        self._zpe_dict = {}
        self._enthalpy_dict = {}
        self._entropy_dict = {}

        # set logger object
        self.logger = logging.getLogger('model.correctors.ThermodynamicCorrector')

    def shomate_correction(self):
        """
        Function to get energy correction dict for gas in model using Shomate equation.
        """
        # {{{
        gas_names = self._owner.gas_names()
        temperature = self._owner.temperature()
        temperature_ref = 298.15

        def H(T, params):
            A, B, C, D, E, F, G, H = params
            t = T/1000.0
            H = A*t + (B/2.0)*t**2 + (C/3.0)*t**3 + (D/4.0)*t**4 - E/t + F - H
            #kJ/mol
            return H

        def S(T, params):
            A, B, C, D, E, F, G, H = params
            t = T/1000.0
            S = A*np.log(t) + B*t + (C/2.0)*t**2 + (D/3.0)*t**3 - E/(2.0*t**2) + G  # J/mol*K
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
                if (gas == gas_key and temperature >= T_min and temperature <= T_max):
                    params = shomate_params[key]
                    Cp_ref = Cp(temperature_ref, params)
                    dH = H(temperature, params) - H(temperature_ref, params)
                    #deltaH(298-T) = shomate(T) - shomate(298)
                    dS = S(temperature, params)
                    dH = (temperature_ref*Cp_ref/1000.0 + dH)*(_kJmol2eV)  # eV
                    #dH = 298*Cp(298) + dH(298-T)
                    dS = dS*(_kJmol2eV/1e3)  # eV/K
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
                    dH = (temperature*Cp_ref/1000.0)*(_kJmol2eV)  # eV
                    dS = dS*(_kJmol2eV/1e3)  # eV/K
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
                msg = 'No Shomate parameters specified for {}'.format(not_there)
                raise ValueError(msg)

        return thermo_dict
        # }}}

    def entropy_correction(self, species_name, m=None, p=None, T=None):
        """
        Function to get free energy constributions from
        translational and internal modes of species in **GAS**.

        Parameters:
        -----------
        species_name: gas molecular formula, str.

        m: absolute molecular mass, species mass by default, float.

        p: partial pressure, model's pressure by default, float.

        T: temperature, model's temperature by default, float.

        Example:
        --------
        >>> m.corrector.entropy_correction('CO')
        >>> -1.1538116935108251
        """
        # {{{
        # Match species in database
        rotation_included = species_name in rotation_temperatures
        vibration_included = species_name in vibration_temperatures

        if not (rotation_included and vibration_included):
            msg_template = '[ {} ] is not in database (thermodynamic_corrector.py)'
            msg = msg_template.format(species_name)
            raise SpeciesError(msg)

        # Set default parameter values.
        if m is None:
            parser = self._owner.parser()
            m = parser.get_molecular_mass(species_name, absolute=True)
        if p is None:
            full_name = species_name + '_g'
            species_definitions = self._owner.species_definitions()
            p = species_definitions[full_name]['pressure']
        if T is None:
            T = self._owner.temperature()

        # Calculate partition functions.

        # Translation partition functions.
        try:
            V = kB_J*T/p
        except ZeroDivisionError:
            msg_template = "Pressure of '{}' is 0.0, please check your input file."
            msg = msg_template.format(species_name)
            raise ZeroDivisionError(msg)

        qt = V*(2*pi*m*kB_J*T/(h_J**2))**(3/2.0)

        # Rotation partition function.
        sigma = rotation_temperatures[species_name]['sigma']
        thetas = rotation_temperatures[species_name]['theta']

        # Linear molecule.
        if len(thetas) == 1:
            theta, = thetas
            ratio = theta/T
            if ratio <= 0.01:
                qr = T/(sigma*theta)
            else:
                qr = T/(sigma*theta)*(1 + theta/(3*T) + theta**2/(15*T**2))
                if ratio >= 0.3:
                    msg_template = "T/theta = {:.3e} is larger than 0.3, big error may be expected"
                    msg = msg_template.format(ratio)
                    self.logger.warning(msg)
        # Nonlinear molecule.
        elif len(thetas) == 3:
            product = reduce(lambda x, y: x*y, thetas)
            qr = (pi)**0.5/sigma*(T**3/(product))**0.5

        # Vibration partition functions.
        thetas = vibration_temperatures[species_name]

        # Linear molecule.
        if len(thetas) == 1:
            theta, = thetas
            qv = exp(-theta/(2*T)) / (1 - exp(-theta/T))
        # Nonlinear molecule.
        else:
            temp_list = [1./(1 - exp(-theta/T)) for theta in thetas]
            qv = reduce(lambda x, y: x*y, temp_list)

        # Molecular partition function.
        q = qt*qr*qv

        return -kB_eV*T*log(q)  # eV
        # }}}


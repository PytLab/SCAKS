from math import exp, pi, sqrt

from kynetix import ModelShell
from kynetix.database.thermo_data import kB_J, kB_eV, h_eV
from kynetix.functions import *


class SolverBase(ModelShell):
    """
    Abstract base class to be herited by other solver classes.
    """

    def __init__(self, owner):
        super(SolverBase, self).__init__(owner)

    @staticmethod
    def get_kTST(Ga, T):
        """
        Static function to get rate constants according to Transition State Theory.

        Parameters:
        -----------
        Ga: free energy barrier, float.

        T: thermodynamics constants, floats.
        """

        kTST = kB_eV*T/h_eV*exp(-Ga/(kB_eV*T))

        return kTST

    @staticmethod
    def get_kCT(Ea, Auc, act_ratio, p, m, T, f=1.0):
        """
        Static function to get rate constant/collision rate
        according to Collision Theory.

        Parameters:
        -----------
        Ea: energy barrier( NOT free energy barrier), float.

        Auc: area of unitcell (m^-2), float.

        act_ratio: area of active sites/area of unitcell, float(<= 1.0).

        p: partial pressure of gas, float.

        m: absolute mass of molecule (kg), float.

        f: factor accounts for a further reduction in the sticking probability,
           if particle with certain initial states are not efficiently steered
           along the MEP, and reflected by a higher barrier, float(<= 1.0).

        T: temperature (K), float.
        """
        # Check parameters.
        if act_ratio > 1.0:
            msg = "active area ratio must be less than 1.0"
            raise ParameterError(msg)

        if f > 1.0:
            msg = "factor f must be less than 1.0"
            raise ParameterError(msg)

        # Sticking coefficient.
        S = f*act_ratio*exp(-Ea/(kB_eV*T))

        # Rate constant.
        kCT = S*(p*Auc)/(sqrt(2*pi*m*kB_J*T))

        return kCT


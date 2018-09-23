import logging
from math import exp, pi, sqrt, log
from operator import mul

from .. import ModelShell
from ..compatutil import reduce
from ..database.thermo_data import kB_J, kB_eV, h_eV, P0
from ..functions import *
from ..parsers.rxn_parser import *
from ..parsers.parser_base import ParserBase


class SolverBase(ModelShell):
    """
    Abstract base class to be herited by other solver classes.
    """

    def __init__(self, owner):
        super(SolverBase, self).__init__(owner)

        # Set logger.
        self.__logger = logging.getLogger("model.solver.SolverBase")

    def get_kTST(self, Ga, T):
        """
        Calculate rate constants according to Transition State Theory.

        Parameters:
        -----------
        Ga: Free energy barrier, float.

        T: Temperature(K), floats.
        """
        if self.__class__.__name__ == "SteadyStateSolver":
            # Convert float numbers to mpfloat.
            kB_eV_mp = self._mpf(kB_eV)
            T_mp = self._mpf(T)
            h_eV_mp = self._mpf(h_eV)
            Ga_mp = self._mpf(Ga)
            kTST = kB_eV_mp*T_mp/h_eV_mp*self._math.exp(-Ga_mp/(kB_eV_mp*T_mp))
        else:
            kTST = kB_eV*T/h_eV*exp(-Ga/(kB_eV*T))

        return kTST

    @staticmethod
    def get_kCT(Ea, Auc, act_ratio, m, T, p=P0, f=1.0):
        """
        Static function to get rate constant/collision rate
        according to Collision Theory.

        Parameters:
        -----------
        Ea: energy barrier( NOT free energy barrier ), float.

        Auc: area of unitcell (m^-2), float.

        act_ratio: area of active sites/area of unitcell, float(<= 1.0).

        p: pressure of gas, float, default value is 101325 Pa (atm).

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

    @staticmethod
    def get_TST_barrier_from_CT(kCT, T):
        """
        Static method to get TST barrier from the rate constant calculated by Collision Theory.

        Parameters:
        -----------
        kCT: Rate constant calculated by Collision Theory, float.

        T: Temperature(K), float.
        """
        Ga = -log(kCT*h_eV/(kB_eV*T))*(kB_eV*T)

        return Ga

    def get_rxn_rates_TST(self, rxn_expression, relative_energies, include_pressure=False):
        """
        Function to get rate constants for an elementary reaction
        using Transition State Theory.

        Parameters:
        -----------
        rxn_expression: The expression of an elementary reaction, str.

        relative_energies: The relative energies for all elementary reactions.

        include_pressure: The flag for whether to include the actual gas pressure
                          (not the pressure of standard condition, 101325.0 Pa).
                          The default value is False, meaning that we calculate the
                          RATE CONSTANT by default.
        """
        # {{{
        # Get the condition for log info output.
        cls_name = self.__class__.__name__
        log_allowed = (self._owner.log_allowed and cls_name == "KMCSolver")

        Gaf, Gar, dG = self._get_relative_energies(rxn_expression, relative_energies)
        if log_allowed:
            self.__logger.info("{} (Gaf={}, Gar={}, dG={})".format(rxn_expression, Gaf, Gar, dG))

        T = self._owner.temperature

        # Use Transition State Theory to get rate constants.
        rf, rr = [self.get_kTST(Ga, T) for Ga in [Gaf, Gar]]
        if log_allowed:
            self.__logger.info("Rates without partial pressure:")
            self.__logger.info("R(forward) = {} s^-1 (Transition State Theory)".format(rf))
            self.__logger.info("R(reverse) = {} s^-1 (Transition State Theory)".format(rr))

        # Include pressures.
        if include_pressure:
            rxn_equation = RxnEquation(rxn_expression)
            f = lambda gas: self._owner.species_definitions[gas]["pressure"]

            # Forward rate.
            adsorption_gases = rxn_equation.adsorption_gases()
            pressures = [f(gas.formula()) for gas in adsorption_gases]
            p = reduce(mul, pressures) if pressures else 1.0
            rf = p*rf

            # Reverse rate.
            desorption_gases = rxn_equation.desorption_gases()
            pressures = [f(gas.formula()) for gas in desorption_gases]
            p = reduce(mul, pressures) if pressures else 1.0
            rr = p*rr

            if log_allowed:
                self.__logger.info("Rates including partial pressure:")
                self.__logger.info("R(forward) = {} s^-1 (Transition State Theory)".format(rf))
                self.__logger.info("R(reverse) = {} s^-1 (Transition State Theory)".format(rr))

        return rf, rr
        # }}}

    def get_rxn_rates_CT(self, rxn_expression, relative_energies, include_pressure=False):
        """
        Function to get rate constants for an elementary reaction
        using Collision Theory wrt adsorption process.

        Parameters:
        -----------
        rxn_expression: The expression of an elementary reaction, str.

        relative_energies: The relative energies for all elementary reactions.

        include_pressure: The flag for whether to include the actual gas pressure
                          (not the pressure of standard condition, 101325.0 Pa).
                          The default value is False, meaning that we calculate the
                          RATE CONSTANT by default.
        """
        # {{{
        # Get the condition for log info output.
        cls_name = self.__class__.__name__
        log_allowed = (self._owner.log_allowed and cls_name == "KMCSolver")

        # Get raw relative energies.
        Gaf, Gar, dG = self._get_relative_energies(rxn_expression, relative_energies)
        if log_allowed:
            self.__logger.info("{} (Gaf={}, Gar={}, dG={})".format(rxn_expression, Gaf, Gar, dG))

        # Get reactants and product types.
        rxn_equation = RxnEquation(rxn_expression)
        formula_list = rxn_equation.to_formula_list()
        istate, fstate = formula_list[0], formula_list[-1]
        is_types = [formula.type() for formula in istate]
        fs_types = [formula.type() for formula in fstate]
        if log_allowed:
            self.__logger.info("species type: {} -> {}".format(is_types, fs_types))

        # Get rate constant.
        T = self._owner.temperature
        Auc = self._owner.unitcell_area

        if not Auc > 0.0:
            raise AttributeError("Invalid unitcell area: {}".format(Auc))

        act_ratio = self._owner.active_ratio

        # Get model corrector.
        corrector = self._owner.corrector
        # Check.
        if type(corrector) == str:
            msg = "No instantialized corrector, try to modify '{}'"
            msg = msg.format(self._owner.setup_file)
            raise SetupError(msg)

        # Adsorption process.
        if "gas" in is_types:
            # Get gas pressure.
            idx = is_types.index("gas")
            formula = istate[idx]
            gas_name = formula.formula()
            m = ParserBase.get_molecular_mass(formula.species(), absolute=True)

            if include_pressure:
                p = P0*self._owner.species_definitions[gas_name]["pressure"]
            else:
                p = P0

            # Use Collision Theory to get forward rate.
            Ea = Gaf
            rf = self.get_kCT(Ea=Ea, Auc=Auc, act_ratio=act_ratio, p=p, m=m, T=T)
            if log_allowed:
                self.__logger.info("R(forward) = {} s^-1 (Collision Theory)".format(rf))

            # Use equilibrium condition to get reverse rate.
#            correction_energy = corrector.entropy_correction(gas_name, m, p, T)
#            stoichiometry = formula.stoichiometry()
#            dG -= stoichiometry*correction_energy
#
#            # Info output.
#            if log_allowed:
#                msg = "Correct dG: {} -> {}".format(dG+correction_energy, dG)
#                self.__logger.info(msg)

#            # Use Equilibrium condition to get reverse rate.
#            K = exp(-dG/(kB_eV*T))
#            rr = rf/K
#            if log_allowed:
#                self.__logger.info("R(reverse) = {} s^-1 (Equilibrium Condition)".format(rr))

            # Use Transition State Theory to get reverse rate.
            # Get equivalent free energy barrier for TST.
            Gaf_TST = self.get_TST_barrier_from_CT(rf, T)

            # Correct the reaction energy dG.
            correction_energy = corrector.entropy_correction(gas_name, m, p, T)
            stoichiometry = formula.stoichiometry()
            dG -= stoichiometry*correction_energy
            Gar_TST = Gaf_TST - dG

            if Gar_TST < 0:
                Gar_TST = 0.0

            rr = self.get_kTST(Gar_TST, T);
            if log_allowed:
                self.__logger.info("Get equivalent TST barrier from CT: {} eV".format(Gaf_TST))
                self.__logger.info("Get equivalent reverse barrier: {} eV".format(Gar_TST))
                self.__logger.info("R(reverse) = {} s^-1 (Transition State Theory)".format(rr))

        # Desorption process.
        elif "gas" in fs_types:
            # Get gas pressure.
            idx = fs_types.index("gas")
            formula = fstate[idx]
            gas_name = formula.formula()

            if include_pressure:
                p = P0*self._owner.species_definitions[gas_name]["pressure"]
            else:
                p = P0

            # Use Collision Theory to get reverse rate.
#            Ea = Gar
            m = ParserBase.get_molecular_mass(formula.species(), absolute=True)
#            rr = self.get_kCT(Ea=Ea, Auc=Auc, act_ratio=act_ratio, p=p, m=m, T=T)
#            if log_allowed:
#                self.__logger.info("R(reverse) = {} s^-1 (Collision Theory)".format(rr))
            correction_energy = corrector.entropy_correction(gas_name, m, p, T)
            stoichiometry = formula.stoichiometry()
            Gar -= stoichiometry*correction_energy
            rr = self.get_kTST(Gar, T)

            # Use Transition State theory to get forward rate.
#            if Gar < 1.0e-10:
#                # NOTE: If the reverse barrier is 0, that means the
#                #       forward barrier depends on the final state energy.
#                correction_energy = corrector.entropy_correction(gas_name, m, p, T)
#                stoichiometry = formula.stoichiometry()
#                Gaf += stoichiometry*correction_energy
            rf = self.get_kTST(Gaf, T)

            if log_allowed:
                self.__logger.info("R(forward) = {} s^-1 (Transition State Theory)".format(rf))

#            # Use equilibrium condition to get forward rate.
#            correction_energy = corrector.entropy_correction(gas_name, m, p, T)
#            stoichiometry = formula.stoichiometry()
#            dG -= stoichiometry*correction_energy
#
#            # Info output.
#            if log_allowed:
#                msg = "Correct dG: {} -> {}".format(dG+correction_energy, dG)
#                self.__logger.info(msg)
#
#            K = exp(dG/(kB_eV*T))
#            rf = rr/K
#
#            if log_allowed:
#                self.__logger.info("R(forward) = {} s^-1 (Equilibrium Condition)".format(rf))

        # Reaction of intermediates.
        else:
            # Use Transition State Theory.
            rf = self.get_kTST(Gaf, T)
            rr = self.get_kTST(Gar, T)
            if log_allowed:
                self.__logger.info("R(forward) = {} s^-1 (Transition State Theory)".format(rf))
                self.__logger.info("R(reverse) = {} s^-1 (Transition State Theory)".format(rr))

        return rf, rr
        # }}}

    def _get_relative_energies(self, rxn_expression, relative_energies):
        """
        Private helper function to get relative energies for an elementary reaction.
        """
        # Get raw relative energies.
        rxn_expressions = self._owner.rxn_expressions
        idx = rxn_expressions.index(rxn_expression)

        Gaf = relative_energies["Gaf"][idx]
        Gar = relative_energies["Gar"][idx]
        dG = relative_energies["dG"][idx]

        return Gaf, Gar, dG


import logging

from kynetix.errors.error import *
from kynetix.database.lattice_data import grid_neighbor_offsets
from kynetix.parsers.rxn_parser import *
from kynetix.parsers.relative_energy_parser import RelativeEnergyParser
from kynetix.solvers.solver_base import SolverBase
from kynetix.utilities.check_utilities import check_process_dict


class KMCParser(RelativeEnergyParser):
    def __init__(self, owner):
        """
        Parser class for KMC simulation.
        """
        super(KMCParser, self).__init__(owner)

        # Set logger.
        self.__logger = logging.getLogger('model.parsers.KMCParser')

    def parse_processes(self, filename=None):
        """
        Function to read processes file and create KMCLibProcess objects.

        Parameters:
        -----------
        filename: The name of processes file, str.

        Returns:
        --------
        A list of KMCLibProcess objects.
        """
        pass

    def __parse_single_process(self, process_dict):
        """
        Private helper function to convert a process dict to KMCLibProcess object.
        """
        # Check process dict.
        process_dict = check_process_dict(process_dict)

        # Check if reaction in rxn_expressions.
        rxn_expressions = self._owner.rxn_expressions()

        if process_dict["reaction"] not in rxn_expressions:
            msg = "'{}' is not in model's rxn_expressions.".format(process_dict["reaction"])
            raise SetupError(msg)

        # Get rate constants.

    def __get_rxn_rates(self, rxn_expression):
        """
        Private helper function to get rate constants for an elementary reaction.
        """
        # {{{
        # Get raw relative energies.
        Gaf, Gar, dG = self.__get_relative_energies(rxn_expression)
        self.__logger.info("{} (Gaf={}, Gar={}, dG={})".format(rxn_expression, Gaf, Gar, dG))

        # Get reactants and product types.
        rxn_equation = RxnEquation(rxn_expression)
        formula_list = rxn_equation.to_formula_list()
        istate, fstate = formula_list[0], formula_list[-1]
        is_types = [formula.type() for formula in istate]
        fs_types = [formula.type() for formula in fstate]
        self.__logger.info("species type: {} -> {}".format(is_types, fs_types))

        # Get rate constant.
        T = self._owner.temperature()
        Auc = self._owner.unitcell_area()
        act_ratio = self._owner.active_ratio()

        # Get model corrector.
        corrector = self._owner.corrector()
        # Check.
        if type(corrector) == str:
            msg = "No instantialized corrector, try to modify '{}'"
            msg = msg.format(self._owner.setup_file())
            raise SetupError(msg)

        # Forward rate.
        self.__logger.info("Calculate forward rate...")

        # Gas participating.
        if "gas" in is_types:
            # Get gas pressure.
            idx = is_types.index("gas")
            formula = istate[idx]
            gas_name = formula.formula()
            p = self._owner.species_definitions()[gas_name]["pressure"]
            
            # Use Collision Theory.
            Ea = Gaf
            m = self.get_molecular_mass(formula.species(), absolute=True)
            rf = SolverBase.get_kCT(Ea, Auc, act_ratio, p, m, T)
            self.__logger.info("R(forward) = {} s^-1 (Collision Theory)".format(rf))
        # No gas participating.
        else:
            # ThermoEquilibrium and gas species in final state.
            if "gas" in fs_types and Gar < 1e-10:
                # Correction energy.
                idx = fs_types.index("gas")
                formula = fstate[idx]
                gas_name = formula.formula()
                p = self._owner.species_definitions()[gas_name]["pressure"]
                m = self.get_molecular_mass(formula.species(), absolute=True)
                correction_energy = corrector.entropy_correction(gas_name, m, p, T)
                Gaf += correction_energy

                # Info output.
                msg = "Correct forward barrier: {} -> {}".format(Gaf-correction_energy, Gaf)
                self.__logger.info(msg)

            rf = SolverBase.get_kTST(Gaf, T)
            self.__logger.info("R(forward) = {} s^-1 (Transition State Theory)".format(rf))

        # Reverse rate.
        self.__logger.info("Calculate reverse rate...")

        # Gas participating.
        if "gas" in fs_types:
            # Get gas pressure.
            idx = fs_types.index("gas")
            formula = fstate[idx]
            gas_name = formula.formula()
            p = self._owner.species_definitions()[gas_name]["pressure"]
            
            # Use Collision Theory.
            Ea = Gar
            m = self.get_molecular_mass(formula.species(), absolute=True)
            rr = SolverBase.get_kCT(Ea, Auc, act_ratio, p, m, T)
            self.__logger.info("R(reverse) = {} s^-1 (Collision Theory)".format(rr))
        # No gas participating.
        else:
            # ThermoEquilibrium and gas species in initial state.
            if "gas" in is_types and Gaf < 1e-10:
                # Correction energy.
                idx = is_types.index("gas")
                formula = istate[idx]
                gas_name = formula.formula()
                p = self._owner.species_definitions()[gas_name]["pressure"]
                m = self.get_molecular_mass(formula.species(), absolute=True)
                correction_energy = corrector.entropy_correction(gas_name, m, p, T)
                Gar += correction_energy

                # Info output.
                msg = "Correct reverse barrier: {} -> {}".format(Gar-correction_energy, Gar)
                self.__logger.info(msg)

            rr = SolverBase.get_kTST(Gar, T)
            self.__logger.info("R(reverse) = {} s^-1 (Transition State Theory)\n".format(rr))

        return rf, rr
        # }}}

    def __get_relative_energies(self, rxn_expression):
        """
        Private helper function to get relative energies for an elementary reaction.
        """
        # Check if parser has relative energies.
        if not self._has_relative_energy:
            msg = "Parser has no relative energies, try parse_data()."
            raise AttributeError(msg)

        # Get raw relative energies.
        rxn_expressions = self._owner.rxn_expressions()
        idx = rxn_expressions.index(rxn_expression)

        Gaf, dG = self._Ga[idx], self._dG[idx]
        Gar = Gaf - dG

        return Gaf, Gar, dG


import logging
import os
from math import exp
from operator import mul

import numpy as np

# KMCLibX.
try:
    from KMCLib import *
except ImportError:
    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    print "!!!                                                   !!!"
    print "!!!          WARNING: KMCLib is not installed         !!!"
    print "!!! Any kMC calculation using KMCLib will be disabled !!!"
    print "!!!                                                   !!!"
    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

from kynetix.database.thermo_data import kB_eV
from kynetix.errors.error import *
from kynetix.database.lattice_data import grid_neighbor_offsets
from kynetix.parsers.rxn_parser import *
from kynetix.parsers.relative_energy_parser import RelativeEnergyParser
from kynetix.solvers.solver_base import SolverBase
from kynetix.utilities.check_utilities import *


class KMCParser(RelativeEnergyParser):
    def __init__(self, owner):
        """
        Parser class for KMC simulation.
        """
        super(KMCParser, self).__init__(owner)

        # Set logger.
        self.__logger = logging.getLogger('model.parsers.KMCParser')

        # Set process reaction mapping.
        self.__process_mapping = []

    def construct_sitesmap(self, filename=None):
        """
        Function to read kmc_site file and create KMCLibSitesMap objects.

        Parameters:
        -----------
        filename: The name of configuration file, str.

        Returns:
        --------
        A KMCLibSitesMap objects.
        """
        # {{{
        # Load data.
        if filename is None:
            filename = "kmc_sites.py"

        possible_site_types = self._owner.possible_site_types

        # Get site number.
        repetitions = self._owner.repetitions
        basis_sites = self._owner.basis_sites
        nsite = reduce(mul, repetitions)*len(basis_sites)

        def init_default_types():
            """
            Nested function to get default site types according to
            lattice repetitions and basis sites.
            """
            default_type = possible_site_types[0]
            return [default_type]*nsite

        # Get site types.
        if not os.path.exists(filename):
            site_types = init_default_types()
        else:
            globs, locs = {}, {}
            execfile(filename, globs, locs)

            if "site_types" not in locs:
                site_types = init_default_types()
            else:
                site_types = locs["site_types"]

        # Check length of site types.
        if len(site_types) != nsite:
            msg = "'site_types' must have {} elements.".format(nsite)
            raise SetupError(msg)

        # Check element type.
        for site_type in site_types:
            if site_type not in possible_site_types:
                msg = "Element '{}' not in possible_site_types '{}'."
                msg = msg.format(site_type, possible_site_types)
                raise SetupError(msg)

        # Construct lattice.
        lattice = self.construct_lattice()

        # Construct sitemap.
        sitesmap = KMCSitesMap(lattice=lattice,
                               types=site_types,
                               possible_types=possible_site_types)

        return sitesmap
        # }}}

    def construct_lattice(self):
        """
        Function to construct KMCLattice object.
        """
        # Construct unitcell.
        cell_vectors = np.array(self._owner.cell_vectors)
        basis_sites = np.array(self._owner.basis_sites)
        unit_cell = KMCUnitCell(cell_vectors=cell_vectors, basis_points=basis_sites)

        # Construct lattice.
        repetitions = self._owner.repetitions
        periodic = self._owner.periodic
        lattice = KMCLattice(unit_cell=unit_cell,
                             repetitions=repetitions,
                             periodic=periodic)

        return lattice

    def parse_configuration(self, filename=None):
        """
        Function to read configuration file and create KMCLibConfiguration objects.

        Parameters:
        -----------
        lattice: The lattice of the configurartion as a KMCLattice.

        filename: The name of configuration file, str.

        Returns:
        --------
        A KMCLibConfiguration objects.
        """
        # {{{
        # Inner function to initialize emtpy lattice.
        def init_empty_types():
            repetitions = self._owner.repetitions
            basis_sites = self._owner.basis_sites
            nsite = reduce(mul, repetitions)*len(basis_sites)
            types = [empty_type for i in xrange(nsite)]
            return types

        possible_element_types = self._owner.possible_element_types

        # Check.
        empty_type = self._owner.empty_type
        check_string(empty_type, possible_element_types, "empty_type")

        # Load types data in file.
        if filename is None:
            filename = "kmc_configuration.py"

        # Use data in file.
        if os.path.exists(filename):
            globs, locs = {}, {}
            execfile(filename, globs, locs)
            if "types" in locs:
                types = locs["types"]
            else:
                types = init_empty_types()
        # Initialize as emtpy lattice.
        else:
            types = init_empty_types()

        # Construct lattice.
        lattice = self.construct_lattice()

        # Instantialize KMCLattice object.
        configuration = KMCConfiguration(lattice=lattice,
                                         types=types,
                                         possible_types=possible_element_types)

        return configuration
        # }}}

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
        if filename is None:
            filename = "kmc_processes.py"

        globs, locs = {}, {}
        execfile(filename, globs, locs)

        # Get all possible process objects.
        self.__process_mapping = []
        all_processes = []
        for process_dict in locs["processes"]:
            processes = self.__parse_single_process(process_dict)
            all_processes.extend(processes)

        if self._owner.log_allowed:
            self.__logger.info("Total {} processes instantalized.".format(len(all_processes)))

        return all_processes

    def __parse_single_process(self, process_dict):
        """
        Private helper function to convert a process dict to KMCLibProcess object.
        """
        # {{{
        # Check process dict.
        process_dict = check_process_dict(process_dict)

        # Check if reaction in rxn_expressions.
        rxn_expressions = self._owner.rxn_expressions

        if process_dict["reaction"] not in rxn_expressions:
            msg = "'{}' is not in model's rxn_expressions.".format(process_dict["reaction"])
            raise SetupError(msg)

        # Check if the elements are in possible elements.
        all_elements = process_dict["elements_before"] + process_dict["elements_after"]
        possible_elements = self._owner.possible_element_types
        for element in all_elements:
            if element not in possible_elements:
                msg = "Element '{}' in process not in possible types {}"
                msg = msg.format(element, possible_elements)
                raise SetupError(msg)

        # Get rate constants.
        rf, rr = self.__get_rxn_rates(process_dict["reaction"])

        # Get process fast flag, False by default.
        fast = process_dict.get("fast", False)

        # Get process redist flag, Falst by default.
        redist = process_dict.get("redist", False)

        # Get process redist species.
        redist_species = process_dict.get("redist_species", None)

        # Get KMCLibProcess objects.
        processes = []

        for basis_site in process_dict["basis_sites"]:
            for coordinates in process_dict["coordinates_group"]:
                if self._owner.log_allowed:
                    self.__logger.info("Coordinates = {}".format(coordinates))
                    self.__logger.info("Basis site = {}".format(basis_site))

                # Forward process.
                fprocess = KMCProcess(coordinates=coordinates,
                                      elements_before=process_dict["elements_before"],
                                      elements_after=process_dict["elements_after"],
                                      basis_sites=[basis_site],
                                      rate_constant=rf,
                                      fast=fast,
                                      redist=redist,
                                      redist_species=redist_species)
                processes.append(fprocess)

                # Add process reaction mapping.
                if not fast:
                    process_mapping = "{}(->)".format(process_dict["reaction"])
                    self.__process_mapping.append(process_mapping)

                # Info output.
                if self._owner.log_allowed:
                    self.__logger.info("Forward elements changes:")
                    self.__logger.info("    /{}".format(process_dict["elements_before"]))
                    self.__logger.info("    \{}".format(process_dict["elements_after"]))

                # --------------------------------------------------------------
                # NOTE: If the proess is a redistribution process which is only
                #       used to re-scatter the fast species, its reverse process
                #       would not be parsed.
                # --------------------------------------------------------------

                # Reverse process.
                if not redist:
                    rprocess = KMCProcess(coordinates=coordinates,
                                          elements_before=process_dict["elements_after"],
                                          elements_after=process_dict["elements_before"],
                                          basis_sites=[basis_site],
                                          rate_constant=rr,
                                          fast=fast)
                    processes.append(rprocess)

                # Add process reaction mapping.
                if not fast:
                    process_mapping = "{}(<-)".format(process_dict["reaction"])
                    self.__process_mapping.append(process_mapping)

                # Info output.
                if not redist and self._owner.log_allowed:
                    self.__logger.info("Reverse elements changes:")
                    self.__logger.info("    /{}".format(process_dict["elements_after"]))
                    self.__logger.info("    \{}".format(process_dict["elements_before"]))

        if self._owner.log_allowed:
            self.__logger.info("\n")

        return processes
        # }}}

    def __get_rxn_rates(self, rxn_expression):
        """
        Private helper function to get rate constants for an elementary reaction.
        """
        # {{{
        # Get raw relative energies.
        Gaf, Gar, dG = self.__get_relative_energies(rxn_expression)
        if self._owner.log_allowed:
            self.__logger.info("{} (Gaf={}, Gar={}, dG={})".format(rxn_expression, Gaf, Gar, dG))

        # Get reactants and product types.
        rxn_equation = RxnEquation(rxn_expression)
        formula_list = rxn_equation.to_formula_list()
        istate, fstate = formula_list[0], formula_list[-1]
        is_types = [formula.type() for formula in istate]
        fs_types = [formula.type() for formula in fstate]
        if self._owner.log_allowed:
            self.__logger.info("species type: {} -> {}".format(is_types, fs_types))

        # Get rate constant.
        T = self._owner.temperature
        Auc = self._owner.unitcell_area
        act_ratio = self._owner.active_ratio

        # Get model corrector.
        corrector = self._owner.corrector
        # Check.
        if type(corrector) == str:
            msg = "No instantialized corrector, try to modify '{}'"
            msg = msg.format(self._owner.setup_file)
            raise SetupError(msg)

        # Forward rate.

        # Gas participating.
        if "gas" in is_types:
            # Get gas pressure.
            idx = is_types.index("gas")
            formula = istate[idx]
            gas_name = formula.formula()
            p = self._owner.species_definitions[gas_name]["pressure"]
            
            # Use Collision Theory.
            Ea = Gaf
            m = self.get_molecular_mass(formula.species(), absolute=True)
            rf = SolverBase.get_kCT(Ea, Auc, act_ratio, p, m, T)
            if self._owner.log_allowed:
                self.__logger.info("R(forward) = {} s^-1 (Collision Theory)".format(rf))
        # No gas participating.
        else:
            # ThermoEquilibrium and gas species in final state.
            if "gas" in fs_types and Gar < 1e-10:
                # Correction energy.
                idx = fs_types.index("gas")
                formula = fstate[idx]
                gas_name = formula.formula()
                p = self._owner.species_definitions[gas_name]["pressure"]
                m = self.get_molecular_mass(formula.species(), absolute=True)
                correction_energy = corrector.entropy_correction(gas_name, m, p, T)
                Gaf += correction_energy

                # Info output.
                msg = "Correct forward barrier: {} -> {}".format(Gaf-correction_energy, Gaf)
                if self._owner.log_allowed:
                    self.__logger.info(msg)

            rf = SolverBase.get_kTST(Gaf, T)
            if self._owner.log_allowed:
                self.__logger.info("R(forward) = {} s^-1 (Transition State Theory)".format(rf))

        # Reverse rate.

        # Gas participating.
        if "gas" in fs_types:
            # Get gas pressure.
            idx = fs_types.index("gas")
            formula = fstate[idx]
            gas_name = formula.formula()
            p = self._owner.species_definitions[gas_name]["pressure"]
            
            # Use Collision Theory.
            Ea = Gar
            m = self.get_molecular_mass(formula.species(), absolute=True)
            rr = SolverBase.get_kCT(Ea, Auc, act_ratio, p, m, T)
            if self._owner.log_allowed:
                self.__logger.info("R(reverse) = {} s^-1 (Collision Theory)".format(rr))
        # No gas participating.
        else:
            # Check if it is an adsorption process.
            if "gas" in is_types and Gaf < 1e-10:
                # Correction energy.
                idx = is_types.index("gas")
                formula = istate[idx]
                gas_name = formula.formula()
                p = self._owner.species_definitions[gas_name]["pressure"]
                m = self.get_molecular_mass(formula.species(), absolute=True)
                correction_energy = corrector.entropy_correction(gas_name, m, p, T)
                dG -= correction_energy

                # Info output.
                if self._owner.log_allowed:
                    msg = "Correct dG: {} -> {}".format(dG+correction_energy, dG)
                    self.__logger.info(msg)

                # Use Equilibrium condition to get reverse rate.
                K = exp(-dG/(kB_eV*T))
                rr = rf/K
                if self._owner.log_allowed:
                    self.__logger.info("R(reverse) = {} s^-1 (Equilibrium Condition)".format(rr))
            else:
                # Use Transition State Theory.
                rr = SolverBase.get_kTST(Gar, T)
                if self._owner.log_allowed:
                    self.__logger.info("R(reverse) = {} s^-1 (Transition State Theory)".format(rr))

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
        rxn_expressions = self._owner.rxn_expressions
        idx = rxn_expressions.index(rxn_expression)

        Gaf, dG = self._Ga[idx], self._dG[idx]
        Gar = Gaf - dG

        return Gaf, Gar, dG

    def process_mapping(self):
        """
        Query function for process reaction type mapping.
        """
        return self.__process_mapping


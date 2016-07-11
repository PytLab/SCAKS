import copy
import re
import logging

import numpy as np

from kynetix import ModelShell
from kynetix import mpi_master
from kynetix.functions import *
from kynetix.errors.error import *
from kynetix.database.elements_data import *
from kynetix.parsers.rxn_parser import *


class ParserBase(ModelShell):
    '''
    class to operate and analyse rxn equations and rxn lists.
    '''
    def __init__(self, owner):
        """
        A class acts as a base class to be inherited by other
        parser classes, it is not functional on its own.
        """
        super(ParserBase, self).__init__(owner)

        # Set elementary parse regex(compiled)
        self.__regex_dict = {}

        # Parser's species definition
        # NOTE: parser's species definitions is the reference of model's.
        self.__species_definitions = owner._KineticModel__species_definitions

        # Set logger.
        self.__logger = logging.getLogger("model.parser.ParserBase")

    def parse_elementary_rxns(self, elementary_rxns):
        """
        Parse all elementary rxn equations.

        Parameters:
        -----------
        elementary_rxns: A list of elementary reaction strings.

        Return:
        -------
        A tuple of elementary related attributes:
            (adsorbate_names,
             gas_names,
             liquid_names,
             site_names,
             transition_state_names,
             elementary_rxns_list)

        Example:
        --------
        >>> elementary_rxns = ['CO_g + *_s -> CO_s',
                               'O2_g + 2*_s -> 2O_s',
                               'CO_s + O_s <-> CO-O_2s -> CO2_g + 2*_s']
        >>> parser.parse_elementary_rxns(elementary_rxns)
        """
        elementary_rxns_list = []
        adsorbate_names = []
        gas_names = []
        liquid_names = []
        site_names = []
        transition_state_names = []

        for equation_str in self._owner.rxn_expressions():
            # debug info
            if mpi_master:
                self.__logger.debug('parsing [ %s ]', equation_str)

            # Get RxnEquation object.
            equation = RxnEquation(equation_str)

            # Check conservation.
            equation.check_conservation()

            # Gather all information for a reaction expression.
            rxn_list = []
            state_list = equation.tolist()

            # For each state object.
            for state in state_list:
                rxn_list.append(state.tolist())

                # for each formula object.
                for formula in state.tolist():
                    site_dict = formula.get_sites_dict()
                    species_site = formula.species_site()

                    # Gas.
                    if "g" in site_dict:
                        gas_names.append(species_site)
                    # Liquid.
                    elif "l" in site_dict:
                        liquid_names.append(species_site)
                    # Transition state.
                    elif "-" in formula.species():
                        transition_state_names.append(species_site)
                    # If not empty site, then adsorbate.
                    elif "*" not in formula.species():
                        adsorbate_names.append(species_site)

                    # Site names.
                    site_names.extend([s for s in site_dict.keys() if s not in ("g", "l")])

            # Append reaction list.
            elementary_rxns_list.append(rxn_list)

        # Merge duplicates in lists
        adsorbate_names = tuple(sorted(list(set(adsorbate_names))))
        gas_names = tuple(sorted(list(set(gas_names))))
        liquid_names = tuple(sorted(list(set(liquid_names))))
        site_names = tuple(sorted(list(set(site_names))))
        transition_state_names = tuple(sorted(list(set(transition_state_names))))

        # Return.
        return (adsorbate_names,
                gas_names,
                liquid_names,
                site_names,
                transition_state_names,
                elementary_rxns_list)

    def get_stoichiometry_matrices(self):
        """
        Go through elementary_rxns_list, return sites stoichiometry matrix,
        reactants and products stoichiometry matrix.

        Returns:
        site_matrix: coefficients matrix for intermediates,
                     if species is on the left of arrow, the entry
                     will be positive, vice-versa.
                     row vector: [*, *self.adsorbate_names], numpy.matrix.

        reapro_matrix: coefficients matrix for reactants and product,
                       if species is on the left of arrow, the entry
                       will be positive, vice-versa.
                       row vector: [*self.gas_names], numpy.matrix.
        """
        # Site and adsorbate names.
        sites_names = (['*_'+site_name for site_name in self._owner.site_names()] +
                       list(self._owner.adsorbate_names()))

        # Reactant and product names.
        reapro_names = list(self._owner.gas_names() + self._owner.liquid_names())

        # Initialize matrices.
        m = len(self._owner.elementary_rxns_list())
        n_s, n_g = len(sites_names), len(reapro_names)
        site_matrix = np.matrix(np.zeros((m, n_s)))
        reapro_matrix = np.matrix(np.zeros((m, n_g)))

        rxns_list = self._owner.elementary_rxns_list()

        # Go through all elementary equations.
        for i, rxn_list in enumerate(rxns_list):
            # Initial state.
            for formula in rxn_list[0]:
                stoich = formula.stoichiometry()
                species_site = formula.species_site()

                # Empty site.
                if species_site in sites_names:
                    j = sites_names.index(species_site)
                    site_matrix[i, j] += stoich

                # Adsorbate.
                if species_site in reapro_names:
                    j = reapro_names.index(species_site)
                    reapro_matrix[i, j] += stoich

            # Final state.
            for formula in rxn_list[-1]:
                stoich = formula.stoichiometry()
                species_site = formula.species_site()

                # Empty site.
                if species_site in sites_names:
                    j = sites_names.index(species_site)
                    site_matrix[i, j] -= stoich

                # Adsorbate.
                if species_site in reapro_names:
                    j = reapro_names.index(species_site)
                    reapro_matrix[i, j] -= stoich

        return site_matrix, reapro_matrix

    def get_total_rxn_equation(self):
        """
        Function to get total reaction expression of the kinetic model.
        """
        site_matrix, reapro_matrix = self.get_stoichiometry_matrices()

        # Helper function to get null space of transposition of site_matrix.
        def null(A, eps=1e-10):
            u, s, vh = np.linalg.svd(A, full_matrices=1, compute_uv=1)
            null_space = np.compress(s <= eps, vh, axis=0)
            return null_space.T

        x = null(site_matrix.T)  # basis of null space
        if not x.any():  # x is not empty
            raise ValueError('Failed to get basis of nullspace.')
        x = map(abs, x.T.tolist()[0])
        #convert entries of x to integer
        min_x = min(x)
        x = [round(i/min_x, 1) for i in x]
        x = np.matrix(x)
        total_coefficients = (x*reapro_matrix).tolist()[0]

        # cope with small differences between coeffs
        abs_total_coefficients = map(abs, total_coefficients)
        min_coeff = min(abs_total_coefficients)
        total_coefficients = [int(i/min_coeff) for i in total_coefficients]

        # create total rxn expression
        reactants_list, products_list = [], []
        reapro_names = self._owner.gas_names() + self._owner.liquid_names()
        for sp_name in reapro_names:
            idx = reapro_names.index(sp_name)
            coefficient = total_coefficients[idx]
            if coefficient < 0:  # for products
                coefficient = abs(int(coefficient))
                if coefficient == 1:
                    coefficient = ''
                else:
                    coefficient = str(coefficient)
                products_list.append(coefficient + sp_name)
            else:  # for reactants
                coefficient = int(coefficient)
                if coefficient == 1:
                    coefficient = ''
                else:
                    coefficient = str(coefficient)
                reactants_list.append(coefficient + sp_name)

        # get total rxn list and set it as an attr of model
        total_rxn_list = [reactants_list, products_list]
        reactants_expr = ' + '.join(reactants_list)
        products_expr = ' + '.join(products_list)
        total_rxn_equation = reactants_expr + ' -> ' + products_expr

        # Check conservation.
        RxnEquation(total_rxn_equation).check_conservation()

        # If check passed, return.
        return total_rxn_equation

    @staticmethod
    def get_molecular_mass(species_name, absolute=False):
        """
        Static function to get relative/absolute molecular mass.

        Parameters:
        -----------
        species_name: name of the molecule species, str.

        absolute: return absolute mass or not(default), bool.

        Example:
        --------
        >>> m.parser.get_molecular_mass('CH4')
        >>> 16.04246
        >>> m.parser.get_molecular_mass('CH4', absolute=True)
        >>> 2.6639131127638393e-26
        """
        elements = string2symbols(species_name)

        # get molecular total relative mass
        molecular_mass = 0.0
        for element in elements:
            if not element in chem_elements:
                msg = 'Element [ %s ] not in database' % element
                raise ElementSearchingError(msg)
            element_mass = chem_elements[element]['mass']
            molecular_mass += element_mass

        if absolute:
            return amu*molecular_mass
        else:
            return molecular_mass

    def _get_state_energy(self, state):
        """
        Protected helper function to get state energy.

        Parameters:
        -----------
        state: An object of ChemState.

        Returns:
        --------
        Absolute free energy of the state, float.
        """
        # Extract species info.
        species_site_dict = state.get_species_site_dict()
        site_dict = state.get_sites_dict()
        formula_list = state.tolist()

        species_definitions = self._owner.species_definitions()
        energy = 0.0

        for formula in formula_list:
            n = formula.stoichiometry()
            species_site = formula.species_site()
            site = formula.site()

            # Adsorbate.
            if "*" not in species_site:
                energy += n*species_definitions[species_site]["formation_energy"]
            # Site.
            else:
                energy += n*species_definitions[site]["formation_energy"]

        return energy

    def get_single_relative_energies(self, rxn_expression):
        """
        Function to get relative energies for an elementary reaction:
            forward barrier,
            reverse barrier,
            reaction energy

        Parameters:
        -----------
        rxn_expression: elementary reaction expression, str.

        Returns:
        --------
        f_barrier: forward barrier.
        r_barrier: reverse barrier.
        reaction_energy: reaction energy.
        """
        # {{{
        # Check.
        if not self._owner.has_absolute_energy():
            msg = "Absolute energie are needed for getting barriers."
            raise AttributeError(msg)

        # Get RxnEquation object.
        rxn_equation = RxnEquation(rxn_expression)

        # Get free energy for states
        G_IS, G_TS, G_FS = 0.0, 0.0, 0.0

        # State list.
        states = rxn_equation.tolist()

        # IS energy.
        G_IS = self._get_state_energy(states[0])

        # FS energy.
        G_FS = self._get_state_energy(states[-1])

        # TS energy.
        if len(states) == 2:
            G_TS = max(G_IS, G_FS)

        if len(states) == 3:
            G_TS = self._get_state_energy(states[1])

        # Get relative energies.
        f_barrier = G_TS - G_IS
        r_barrier = G_TS - G_FS
        reaction_energy = G_FS - G_IS

        return f_barrier, r_barrier, reaction_energy
        # }}}

    def get_relative_from_absolute(self):
        """
        Function to set relative energies from absolute energies.
        """
        Gafs, Gars, dGs = [], [], []

        for rxn_expression in self._owner.rxn_expressions():
            Gaf, Gar, dG = self.get_single_relative_energies(rxn_expression)
            Gafs.append(Gaf)
            Gars.append(Gar)
            dGs.append(dG)

        relative_energies = dict(Gar=Gars, Gaf=Gafs, dG=dGs)

        return relative_energies

    def regex_dict(self):
        """
        Query function for regress expression dictionary.
        """
        return self.__regex_dict

    @return_deepcopy
    def species_definitions(self):
        """
        Query function for parser's species definitions.
        """
        # Use deep copy to avoid modification of the model's attribution.
        return self.__species_definitions

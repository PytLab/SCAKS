import os
import logging

from scipy.linalg import solve

from kynetix import mpi_master
from kynetix.parsers.parser_base import *
from kynetix.functions import *


class RelativeEnergyParser(ParserBase):
    def __init__(self, owner):
        '''
        A class to parse relative energies
        covert them to generalize formation energies.
        '''
        super(RelativeEnergyParser, self).__init__(owner)

        # Set logger.
        self.__logger = logging.getLogger('model.parsers.RelativeEnergyParser')

        #intialize generalized formation energy dict
        self.__G_dict = {}

        # Flags.
        self._has_relative_energy= False
        self._has_absolute_energy = False

    def __chk_data_validity(self):
        """
        Private helper function to check relative energy data validity.
        """

        if hasattr(self, "_Ga"):
            if len(self._Ga) != len(self._owner.rxn_expressions()):
                raise ValueError("Invalid shape of Ga.")

        if hasattr(self, "_dG"):
            if len(self._dG) != len(self._owner.rxn_expressions()):
                raise ValueError("Invalid shape of dG.")

        return

    ####### Use matrix to get generalized formation energy ########

    def __get_unknown_species(self):
        # {{{
        """
        Private function to get species whose free energy is unknown.
        """
        # Get all species.
        all_species = (self._owner.site_names() + self._owner.gas_names() +
                       self._owner.liquid_names() + self._owner.adsorbate_names() +
                       self._owner.transition_state_names())
        all_species = list(all_species)

        # Remove known species.
        for known_species in self._owner.ref_species():
            all_species.remove(known_species)
            self.__G_dict.setdefault(known_species, 0.0)

        unknown_species = all_species

        # Debug info output
        if mpi_master:
            self.__logger.debug('unknown species = {}'.format(unknown_species))
            self.__logger.debug('{} unknown species.'.format(len(unknown_species)))

        return unknown_species
        # }}}

    def __get_unknown_coeff_vector(self, rxn_expression):
        # {{{
        """
        Private helper function to get coefficient vector for unknown species.

        Parameters:
        -----------
        elementary_rxn_list: list of an elementary reaction.
                             e.g. [['O2_s', 'NO_g'], ['*_s', 'ON-OO_s'], ['*_s', 'ONOO_s']]

        Returns:
        --------
        coefficient vectors, Ga, dG
        e.g. ([[0, 0, -1, 0, 0, 0, 0, 0, 1], [0, 0, -1, 1, 0, 0, 0, 0, 0]], 0.655, -0.455)

        NOTE:
        -----
        The shape of coefficient vector is the same with that of unknown_species.
        """
        # Get relative energies for the elementary reaction.
        idx = self._owner.rxn_expressions().index(rxn_expression)
        Ga, dG = self._Ga[idx], self._dG[idx]

        # Get ChemState object list.
        rxn_equation = RxnEquation(rxn_expression)
        states = rxn_equation.tolist()

        # Unknown species names.
        unknown_species = self.__get_unknown_species()

        coeff_vects = []

        # Equation for E(TS) - E(IS) = Ea.
        if Ga != 0 and len(states) != 2:
            # Get ts coefficient vector
            istate, tstate = states[0], states[1]
            is_species_sites = istate.get_species_site_list()
            ts_species_sites = tstate.get_species_site_list()
            coeff_vect = []

            for unknown in unknown_species:
                if unknown in is_species_sites:
                    formula_idx = is_species_sites.index(unknown)
                    formula = istate.tolist()[formula_idx]
                    coeff = -formula.stoichiometry()
                elif unknown in ts_species_sites:
                    formula_idx = ts_species_sites.index(unknown)
                    formula = tstate.tolist()[formula_idx]
                    coeff = formula.stoichiometry()
                else:
                    coeff = 0
                coeff_vect.append(coeff)
            coeff_vects.append(coeff_vect)

        # Equation for E(FS) - E(IS) = dE.
        istate, fstate = states[0], states[-1]
        is_species_sites = istate.get_species_site_list()
        fs_species_sites = fstate.get_species_site_list()
        coeff_vect = []

        for unknown in unknown_species:
            if unknown in is_species_sites:
                formula_idx = is_species_sites.index(unknown)
                formula = istate.tolist()[formula_idx]
                coeff = -formula.stoichiometry()
            elif unknown in fs_species_sites:
                formula_idx = fs_species_sites.index(unknown)
                formula = fstate.tolist()[formula_idx]
                coeff = formula.stoichiometry()
            else:
                coeff = 0
            coeff_vect.append(coeff)

        coeff_vects.append(coeff_vect)

        # Debug info output
        if mpi_master:
            self.__logger.debug('elementary rxn: {}'.format(rxn_expression))
            self.__logger.debug('unknown species coeffs: {}'.format(coeff_vects))

        if Ga:
            if mpi_master:
                self.__logger.debug('Ga, dG: {}'.format([Ga, dG]))
            return coeff_vects, [Ga, dG]
        else:
            if mpi_master:
                self.__logger.debug('dG: {}'.format(dG))
            return coeff_vects, [dG]
        # }}}

    def __convert_data(self):
        # {{{
        """
        Solve Axb equation to get value of generalized free energies.
        Convert relative energies to absolute energies,
        then pass values to its owner.
        """

        # Get coefficients matrix A and values vector b
        A, b = [], []
        for rxn_expression in self._owner.rxn_expressions():
            coeff_vects, value = self.__get_unknown_coeff_vector(rxn_expression)
            A.extend(coeff_vects)
            b.extend(value)

        A, b = np.matrix(A), np.matrix(b).reshape(-1, 1)
        # Output debug info
        if mpi_master:
            self.__logger.debug('A = \n{}'.format(str(A)))
            self.__logger.debug('A.shape = {}'.format(str(A.shape)))
        row, col = A.shape
        if row != col:
            if mpi_master:
                self.__logger.warning('!!! %d equations for %d variables !!!' +
                                  'please check your [ ref_species ] in [ %s ]',
                                  row, col, self._owner.setup_file())
        if mpi_master:
            self.__logger.debug('b = \n{}'.format(str(b)))
            self.__logger.debug('b.shape = {}'.format(str(b.shape)))

        # Solve the equations.
        #x = A.I*b
        x = solve(A, b)

        # Output debug info
        if mpi_master:
            self.__logger.debug('x = \n{}'.format(str(x)))
            self.__logger.debug('x.shape = {}'.format(str(x.shape)))

        # Convert column vector to list
        x = x.reshape(1, -1).tolist()[0]

        # Put values to __G_dict
        unknown_species = self.__get_unknown_species()
        for sp_name, G in zip(unknown_species, x):
            self.__G_dict.setdefault(sp_name, G)

        return
        # }}}

    @staticmethod
    def __compare_relative_energies(dict1, dict2):
        """
        Private static method for compare two relative energies dicts.
        """
        for key in dict1:
            list1 = dict1[key]
            list2 = dict2[key]
            # Compare.
            for item1, item2 in zip(list1, list2):
                if (item1 - item2) > 10e-10:
                    return False

        return True

    def __get_relative_from_relative(self):
        """
        Private helper function to get relative energies dict from relative energy.
        """

        # Correct forward barrier.
        Gafs, dGs = [], []
        for Gaf, dG in zip(self._Ga, self._dG):
            # -------------------------------------------------------
            # NOTE: if the reaction has no barrier and is exothermic,
            #       the forward barrier should be positive and
            #       reaction heat should be the same with forward
            #       barrier.
            # -------------------------------------------------------
            if abs(Gaf - 0.0) < 1e-10 and dG > 0.0:
                Gaf = dG
            Gafs.append(Gaf)
            dGs.append(dG)

        relative_energies = dict(Gaf=Gafs, dG=dGs)

        # Get reverse barriers.
        Gars = [Ga - dG for Ga, dG in zip(Gafs, dGs)]
        relative_energies.setdefault("Gar", Gars)

        return relative_energies

    def parse_data(self, relative=False, filename="./rel_energy.py"):
        # {{{
        """
        Put generalized formation energy into species_definitions.

        Parameters:
        -----------
        relative: Only relative energies or not, default to be False.
                  If true, no absolute energies would be put to species definitions.

        filename: The filename of relative energies data.

        NOTE: This function will change the species definition of model.
              If relative is true, then the relative_energies of model will be set.
              If relative is false, formation energy in species_definitions
              and realtive energies will be set.
        """
        # Read relative energy data file.
        if os.path.exists(filename):
            globs, locs = {}, {}
            execfile(filename, globs, locs)

            # Set variables in data file as attr of parser
            for key in locs:
                setattr(self, "_" + key, locs[key])
        else:
            raise IOError("{} is not found.".format(filename))

        # Check data shape.
        self.__chk_data_validity()

        if not relative:
            # Get absolute energy for each species
            if not self.__G_dict:
                self.__convert_data()

            # NOTE: use the REFERENCE of model's species definitions.
            attribute_name = mangled_name(self._owner, "species_definitions")
            species_definitions = getattr(self._owner, attribute_name)

            # Update reference species formation energies.
            for species in self._owner.ref_species():
                if species in species_definitions:
                    species_definitions[species]["formation_energy"] = 0.0
                else:
                    species_definitions.setdefault(species, {"formation_energy": 0.0})

            # Add unknown species formation energies.
            for species in self.__G_dict:
                if species not in species_definitions:
                    energy_dict = {"formation_energy": self.__G_dict[species]}
                    species_definitions.setdefault(species, energy_dict)
                else:
                    species_definitions[species]["formation_energy"] = self.__G_dict[species]

            # Set flag.
            attribute_name = mangled_name(self._owner, "has_absolute_energy")
            setattr(self._owner, attribute_name, True)

        # Get relative energy.
        if '_dG' and '_Ga' in self.__dict__:
            self._has_relative_energy= True
        else:
            raise IOError(("No relative energy was read, " +
                           "please check the '{}'").format(filename))

        # Get relative energies.
        relative_energies = self.__get_relative_from_relative()
        attribute_name = mangled_name(self._owner, "relative_energies")
        setattr(self._owner, attribute_name, relative_energies)

        # Set flags.
        attribute_name = mangled_name(self._owner, "has_relative_energy")
        setattr(self._owner, attribute_name, True)

        # Check the consistency of relative energies.
        if self._owner.has_absolute_energy():
            relative_energies_from_absolute = self.get_relative_from_absolute()
            consistent = self.__compare_relative_energies(relative_energies_from_absolute,
                                                          self._owner.relative_energies())
            if not consistent:
                msg = ("relative energies from solving equations " +
                       "are not equal to that in data files.")
                raise ValueError(msg)

        return
        # }}}

    @return_deepcopy
    def Ga(self):
        """
        Query function for Ga.
        """
        return self._Ga

    @return_deepcopy
    def dG(self):
        """
        Query function for dG.
        """
        return self._dG

    @return_deepcopy
    def Ea(self):
        """
        Query function for Ea.
        """
        return self.__Ea

    @return_deepcopy
    def dE(self):
        """
        Query function for dE.
        """
        return self.__dE


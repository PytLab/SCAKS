import os
import logging

from scipy.linalg import solve

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
        self.__has_relative_energy = False
        self.__has_absolute_energy = False

    def __chk_data_validity(self):
        """
        Private helper function to check relative energy data validity.
        """

        if hasattr(self, "_RelativeEnergyParser__Ga"):
            if len(self.__Ga) != len(self._owner.elementary_rxns_list()):
                raise ValueError("Invalid shape of Ga.")

        if hasattr(self, "_RelativeEnergyParser__dG"):
            if len(self.__dG) != len(self._owner.elementary_rxns_list()):
                raise ValueError("Invalid shape of dG.")

        return

    ####### Use matrix to get generalized formation energy ########

    def __get_unknown_species(self):
        "Get species whose free energy is unknown."
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
        self.__logger.debug('unknown species = {}'.format(unknown_species))
        self.__logger.debug('{} unknown species.'.format(len(unknown_species)))

        return unknown_species

    def __list2dict(self, state_list):
        """
        Expect a state list, e.g. ['*_s', 'NO_g']
        return a corresponding dict, e.g. {'*_s': 1, 'NO_g': 1}.
        """
        state_dict = {}
        for sp_str in state_list:
            stoichiometry, species_name = self.split_species(sp_str)
            state_dict.setdefault(species_name, stoichiometry)

        return state_dict

    def __get_unknown_coeff_vector(self, elementary_rxn_list):
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
        idx = self._owner.elementary_rxns_list().index(elementary_rxn_list)
        Ga, dG = self.__Ga[idx], self.__dG[idx]

        unknown_species = self.__get_unknown_species()

        coeff_vects = []
        if Ga != 0 and len(elementary_rxn_list) != 2:  # has barrier
            # Get ts coefficient vector
            is_list, ts_list = elementary_rxn_list[0], elementary_rxn_list[1]
            is_dict, ts_dict = self.__list2dict(is_list), self.__list2dict(ts_list)
            coeff_vect = []
            for unknown in unknown_species:
                if unknown in is_dict:
                    coeff = -is_dict[unknown]
                elif unknown in ts_dict:
                    coeff = ts_dict[unknown]
                else:
                    coeff = 0
                coeff_vect.append(coeff)
            coeff_vects.append(coeff_vect)

        # Coefficient vector for dG
        is_list, fs_list = elementary_rxn_list[0], elementary_rxn_list[-1]
        is_dict, fs_dict = self.__list2dict(is_list), self.__list2dict(fs_list)
        coeff_vect = []
        for unknown in unknown_species:
            if unknown in is_dict:
                coeff = -is_dict[unknown]
            elif unknown in fs_dict:
                coeff = fs_dict[unknown]
            else:
                coeff = 0
            coeff_vect.append(coeff)
        coeff_vects.append(coeff_vect)

        # Debug info output
        self.__logger.debug('elementary rxn: {}'.format(elementary_rxn_list))
        self.__logger.debug('unknown species coeffs: {}'.format(coeff_vects))

        if Ga:
            self.__logger.debug('Ga, dG: {}'.format([Ga, dG]))
            return coeff_vects, [Ga, dG]
        else:
            self.__logger.debug('dG: {}'.format(dG))
            return coeff_vects, [dG]

    def __convert_data(self):
        '''
        Solve Axb equation to get value of generalized free energies.
        Convert relative energies to absolute energies,
        then pass values to its owner.
        '''

        # Get coefficients matrix A and values vector b
        A, b = [], []
        for rxn_list in self._owner.elementary_rxns_list():
            coeff_vects, value = self.__get_unknown_coeff_vector(rxn_list)
            A.extend(coeff_vects)
            b.extend(value)

        A, b = np.matrix(A), np.matrix(b).reshape(-1, 1)
        # Output debug info
        self.__logger.debug('A = \n{}'.format(str(A)))
        self.__logger.debug('A.shape = {}'.format(str(A.shape)))
        row, col = A.shape
        if row != col:
            self.__logger.warning('!!! %d equations for %d variables !!!' +
                                  'please check your [ ref_species ] in [ %s ]',
                                  row, col, self._owner.setup_file())
        self.__logger.debug('b = \n{}'.format(str(b)))
        self.__logger.debug('b.shape = {}'.format(str(b.shape)))

        # Solve the equations.
        #x = A.I*b
        x = solve(A, b)

        # Output debug info
        self.__logger.debug('x = \n{}'.format(str(x)))
        self.__logger.debug('x.shape = {}'.format(str(x.shape)))

        # Convert column vector to list
        x = x.reshape(1, -1).tolist()[0]

        # Put values to __G_dict
        unknown_species = self.__get_unknown_species()
        for sp_name, G in zip(unknown_species, x):
            self.__G_dict.setdefault(sp_name, G)

        return

    def __get_relative_from_relative(self):
        """
        Private helper function to get relative energies dict from relative energy.
        """
        relative_energies = dict(Gaf=self.__Ga, dG=self.__dG)

        # Get reverse barriers.
        Gars = [Ga - dG for Ga, dG in zip(self.__Ga, self.__dG)]
        relative_energies.setdefault("Gar", Gars)

        return relative_energies

    def parse_data(self, relative=False, filename="./rel_energy.py"):
        '''
        Put generalized formation energy into species_definitions.

        NOTE: This function will change the species definition of model.
        '''
        # Read relative energy data file.
        if os.path.exists(filename):
            globs, locs = {}, {}
            execfile(filename, globs, locs)

            # Set variables in data file as attr of parser
            for key in locs:
                attribute_name = mangled_name(self, key)
                setattr(self, attribute_name, locs[key])
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

            # Update formation energies in model's species definitions.
            for sp_name in self.__G_dict:
                species_definitions[sp_name].setdefault('formation_energy',
                                                        self.__G_dict[sp_name])

            # Set flag.
            attribute_name = mangled_name(self._owner, "has_absolute_energy")
            setattr(self._owner, attribute_name, True)

        # Get relative energy.
        if '_RelativeEnergyParser__dG' and '_RelativeEnergyParser__Ga' in self.__dict__:
            self.__has_relative_energy = True
        else:
            raise IOError(("No relative energy was read, " +
                           "please check the '{}'").format(self.__filename))

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

    @return_deepcopy
    def Ga(self):
        """
        Query function for Ga.
        """
        return self.__Ga

    @return_deepcopy
    def dG(self):
        """
        Query function for dG.
        """
        return self.__dG

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


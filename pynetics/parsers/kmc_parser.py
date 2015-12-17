import logging

from ..errors.error import *
from ..database.lattice_data import grid_neighbor_offsets
from .relative_energy_parser import RelativeEnergyParser


class KMCParser(RelativeEnergyParser):
    def __init__(self, owner):
        '''
        A class to parse KMC related data.
        '''
        RelativeEnergyParser.__init__(self, owner)
        # set logger object
        self.logger = logging.getLogger('model.parsers.KMCParser')

    @staticmethod
    def check_active_site_number(site_number):
        '''
        check whether this code support the site number.
        '''
        # I sincerely hope to remove this boring function in this code someday !-_-
        if site_number > 2:
            msg = 'Do not support multi-sites processing auto-generation' +\
                  'in this version.'
            raise ProcessParsingError(msg)
        if site_number < 1:
            msg = 'No lattice site in the rxn list: %s' % str(elementary_rxn_list)
            raise ProcessParsingError(msg)

    def get_elementary_elements_changes(self, elementary_rxn_list):
        '''
        Function to get local elements changes list of an elementary_rxn_list on grid.

        Parameters:
        -----------
        elementary_rxn_list: an elementary reaction list contains lists of states,
                             list of lists of str.

        Example:
        --------
        >>> l = [['CO_g', '*_s'], ['CO_s']]
        >>> m.parser.get_elementary_elements_changes(l)
        >>> [(['Vac', '*', '*', '*', '*'], ['CO', '*', '*', '*', '*'])]

        '''
        #------------------------------------------------------------
        #  !!! ignore multi-site here, site is 's' by default !!!
        #      may be added in the future
        #------------------------------------------------------------
        reactants, products = elementary_rxn_list[0], elementary_rxn_list[-1]

        # get all elements before and after
        elements_before_list = self.get_local_elements(reactants)
        elements_after_list = self.get_local_elements(products)

        return zip(elements_before_list, elements_after_list)

    def get_local_elements(self, reactants):
        '''
        Function to get all possible local configure elements.

        Parameters:
        -----------
        reactants: list of species, list of str.

        Example:
        --------
        >>> m.parser.get_local_elements(['CO_s', '*_s'])
        >>> [['CO', 'Vac', '*', '*', '*'],
             ['CO', '*', 'Vac', '*', '*'],
             ['CO', '*', '*', 'Vac', '*'],
             ['CO', '*', '*', '*', 'Vac']]
        '''
        # get neighbor offsets
        grid_type = self._owner.grid_type.lower().strip()
        neighbor_offsets = grid_neighbor_offsets[grid_type]

        local_elements_list = []
        active_elements = self.get_active_elements(reactants)

        # one element in the center
        if len(active_elements) == 1:
            local_elements = active_elements + ['*']*len(neighbor_offsets)
            local_elements_list.append(local_elements)

        # interact with one element in neighborhood
        elif len(active_elements) == 2:
            center_element, neighbor_element = active_elements
            for i in xrange(len(neighbor_offsets)):
                surrounds = ['*']*len(neighbor_offsets)
                surrounds[i] = neighbor_element
                local_elements = [center_element] + surrounds
                local_elements_list.append(local_elements)

        return local_elements_list

    def get_active_elements(self, reactants):
        '''
        Function to get active elements list.

        Parameters:
        -----------
        reactants: list of species, list of str.

        Example:
        --------
        >>> m.parser.get_active_elements(['CO_s', '*_s'])
        >>> ['CO', 'Vac']

        '''
        site_number = self.get_occupied_site_number(reactants)
        self.check_active_site_number(site_number)

        strip = lambda sp: sp.split('_')[0]

        # unimolecular adsorption and desorption
        if site_number == 1:
            for sp in reactants:
                if sp.endswith('_s'):
                    if '*' in sp:
                        center_element = 'Vac'
                    else:
                        center_element = strip(sp)
                    break

            return [center_element]

        # reaction, assiciative desorption, dessociative adsorption
        elif site_number == 2:
            active_elements = []
            for sp in reactants:
                if sp.endswith('_s'):
                    stoich, sp_name = self.split_species(sp)
                    if stoich == 2:  # for 2 same species
                        if '*' in sp_name:
                            active_elements = ['Vac', 'Vac']
                        else:
                            active_elements = [strip(sp_name)]*2
                        break
                    elif stoich == 1:
                        if '*' in sp_name:
                            active_elements.append('Vac')
                        else:
                            active_elements.append(strip(sp_name))

            return active_elements

    def get_occupied_site_number(self, reactants, target_site_type='s'):
        '''
        Function to get site number of an elementary reaction occupied.

        Parameters:
        -----------
        reactants: a list of reactant species in a elementary reaction,
                   list of str.

        target_site_type: site type/name, str

        Example:
        --------
        >>> m.parser.get_occupied_site_number(['*_s', 'CO_g'], 's')
        >>> 1

        '''

        site_number = 0
        for species in reactants:
            if '*' in species:  # site
                sites_dict = self.parse_site_expression(species)
                for site_type in sites_dict:
                    if site_type == target_site_type:
                        site_number += sites_dict[site_type]['number']
            else:  # not site
                species_dict = self.parse_species_expression(species)
                for species_type in species_dict:
                    site_type = species_dict[species_type]['site']
                    if site_type == target_site_type:
                        site_number += species_dict[species_type]['number']

        return site_number

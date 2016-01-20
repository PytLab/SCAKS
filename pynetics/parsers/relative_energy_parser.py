import os
import logging

from parser_base import *


class RelativeEnergyParser(ParserBase):
    def __init__(self, owner):
        '''
        A class to parse relative energies
        covert them to generalize formation energies.
        '''

        ParserBase.__init__(self, owner)
        # set tools logger as child of model's
        self.logger = logging.getLogger('model.parsers.RelativeEnergyParser')

        #intialize generalized formation energy dict
        self.G_dict = {}

        #parse in data file
        if os.path.exists('rel_energy.py'):
            globs, locs = {}, {}
            execfile('rel_energy.py', globs, locs)
            # set variables in data file as attr of parser
            for key in locs:
                setattr(self, key, locs[key])
            # pass relative energies to owner
            setattr(self._owner, 'relative_energies', {})
            for key in ['Ga', 'Ea', 'dG', 'dE']:
                if key in locs:
                    self._owner.relative_energies.setdefault(key, locs[key])
        else:
            raise IOError("No rel_energy.py in current path.")

    def chk_data_validity(self):
        '''check data validity.'''

        if (hasattr(self, 'Ga') and
                not (len(self.Ga) == len(self._owner.elementary_rxns_list))):
            raise ValueError("Invalid shape of Ga.")

        if (hasattr(self, 'Ea') and
                not (len(self.Ea) == len(self._owner.elementary_rxns_list))):
            raise ValueError("Invalid shape of Ea.")

        if (hasattr(self, 'dG') and
                not (len(self.dG) == len(self._owner.elementary_rxns_list))):
            raise ValueError("Invalid shape of dG.")

        if (hasattr(self, 'dE') and
                not (len(self.dE) == len(self._owner.elementary_rxns_list))):
            raise ValueError("Invalid shape of dE.")

        return

    ####### Use matrix to get generalized formation energy ########

    def get_unknown_species(self):
        "Get species whose free energy is unknown."
        all_sp = self._owner.site_names + self._owner.gas_names + \
            self._owner.liquid_names + self._owner.adsorbate_names + \
            self._owner.transition_state_names
        all_sp = list(all_sp)

        for known_sp in self._owner.ref_species:
            all_sp.remove(known_sp)
            self.G_dict.setdefault(known_sp, 0.0)

        self.unknowns = all_sp

        # debu info output
        self.logger.debug('unknown species = %s', str(self.unknowns))
        self.logger.debug('%d unknown species.', len(self.unknowns))

        return all_sp

    def list2dict(self, state_list):
        """
        Expect a state list, e.g. ['*_s', 'NO_g']
        return a corresponding dict, e.g. {'*_s': 1, 'NO_g': 1}.
        """
        state_dict = {}
        for sp_str in state_list:
            stoichiometry, species_name = self.split_species(sp_str)
            state_dict.setdefault(species_name, stoichiometry)

        return state_dict

    def get_unknown_coeff_vector(self, elementary_rxn_list):
        """
        Expect a elementary_rxn_list,
        e.g. [['O2_s', 'NO_g'], ['*_s', 'ON-OO_s'], ['*_s', 'ONOO_s']]
        return coefficient vectors, Ga, dG.
        e.g. ([[0, 0, -1, 0, 0, 0, 0, 0, 1], [0, 0, -1, 1, 0, 0, 0, 0, 0]], 0.655, -0.455)

        Note:
        -----
        The shape of coefficient vector is the same with that of unknowns.
        """
        idx = self._owner.elementary_rxns_list.index(elementary_rxn_list)
        Ga, dG = self.Ga[idx], self.dG[idx]

        if not hasattr(self, 'unknowns'):
            self.get_unknown_species()

        coeff_vects = []
        if Ga != 0 and len(elementary_rxn_list) != 2:  # has barrier
            #get ts coefficient vector
            is_list, ts_list = elementary_rxn_list[0], elementary_rxn_list[1]
            is_dict, ts_dict = self.list2dict(is_list), self.list2dict(ts_list)
            coeff_vect = []
            for unknown in self.unknowns:
                if unknown in is_dict:
                    coeff = -is_dict[unknown]
                elif unknown in ts_dict:
                    coeff = ts_dict[unknown]
                else:
                    coeff = 0
                coeff_vect.append(coeff)
            coeff_vects.append(coeff_vect)

        #coefficient vector for dG
        is_list, fs_list = elementary_rxn_list[0], elementary_rxn_list[-1]
        is_dict, fs_dict = self.list2dict(is_list), self.list2dict(fs_list)
        coeff_vect = []
        for unknown in self.unknowns:
            if unknown in is_dict:
                coeff = -is_dict[unknown]
            elif unknown in fs_dict:
                coeff = fs_dict[unknown]
            else:
                coeff = 0
            coeff_vect.append(coeff)
        coeff_vects.append(coeff_vect)

        # debug info output
        self.logger.debug('elementary rxn: %s', str(elementary_rxn_list))
        self.logger.debug('unknown species coeffs: %s', str(coeff_vects))

        if Ga:
            self.logger.debug('Ga, dG: %s', str([Ga, dG]))
            return coeff_vects, [Ga, dG]
        else:
            self.logger.debug('dG: %s', str(dG))
            return coeff_vects, [dG]

    def convert_data(self):
        '''
        Solve Axb equation to get value of generalized free energies.
        Convert relative energies to absolute energies,
        then pass values to its owner.
        '''

        # get coefficients matrix A and values vector b
        A, b = [], []
        for rxn_list in self._owner.elementary_rxns_list:
            coeff_vects, value = self.get_unknown_coeff_vector(rxn_list)
            A.extend(coeff_vects)
            b.extend(value)

        A, b = np.matrix(A), np.matrix(b).reshape(-1, 1)
        self.A = A
        self.b = b
        # output debug info
        self.logger.debug('A = \n%s', str(A))
        self.logger.debug('A.shape = %s', str(A.shape))
        row, col = A.shape
        if row != col:
            self.logger.warning('!!! %d equations for %d variables !!!' +
                                'please check your [ ref_species ] in [ %s ]',
                                row, col, self._owner.setup_file)
        self.logger.debug('b = \n%s', str(b))
        self.logger.debug('b.shape = %s', str(b.shape))

        x = A.I*b  # values for unknowns

        # output debug info
        self.logger.debug('x = \n%s', str(x))
        self.logger.debug('x.shape = %s', str(x.shape))
        #convert column vector to list
        x = x.reshape(1, -1).tolist()[0]

        #put values to G_dict
        for sp_name, G in zip(self.unknowns, x):
            self.G_dict.setdefault(sp_name, G)

        return

    def parse_data(self, relative=False):
        '''
        put generalized formation energy into species_definition.
        '''

        # get relative energy only
        if relative:
            if 'dG' and 'Ga' in self._owner.relative_energies:
                setattr(self._owner, 'has_relative_energy', True) 
                return
            elif 'dE' and 'Ea' in self._owner.relative_energies:
                setattr(self._owner, 'has_relative_energy', True)
                return
            else:
                raise IOError('No relative energy was read, ' +
                              'please check the \'rel_energy.py\'')

        # get absolute energy for each species
        if not self.G_dict:
            self.convert_data()

        for sp_name in self.G_dict:
            sp_dict = self._owner.species_definitions
            sp_dict[sp_name].setdefault('formation_energy', self.G_dict[sp_name])

        setattr(self._owner, 'has_absolute_energy', True)

        return

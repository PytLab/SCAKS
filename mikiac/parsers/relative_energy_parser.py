import os
import logging

from scipy.linalg import solve

from .parser_base import *
from ..functions import *


class RelativeEnergyParser(ParserBase):
    def __init__(self, owner):
        '''
        A class to parse relative energies
        covert them to generalize formation energies.
        '''
        super(RelativeEnergyParser, self).__init__(owner)

        # Set logger.
        self.__logger = logging.getLogger('model.parsers.RelativeEnergyParser')

        # Flags.
        self._has_relative_energy= False
        self._has_absolute_energy = False

    def __chk_data_validity(self, data_dict):
        """
        Private helper function to check relative energy data validity.
        """

        if "Ga" in data_dict:
            if len(data_dict["Ga"]) != len(self._owner.rxn_expressions):
                raise ValueError("Invalid shape of Ga.")

        if "dG" in data_dict:
            if len(data_dict["dG"]) != len(self._owner.rxn_expressions):
                raise ValueError("Invalid shape of dG.")

        return

    def __get_relative_energies(self, data_dict):
        """
        Private helper function to get relative energies from
        relative energies in data file.
        """
        # {{{
        # Correct forward barrier.
        Gafs, dGs = [], []
        for Gaf, dG in zip(data_dict["Ga"], data_dict["dG"]):
            # -------------------------------------------------------
            # NOTE: if the reaction has no barrier and is endothermic,
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
        # }}}

    def parse_data(self, filename="./rel_energy.py"):
        # {{{
        """
        Put generalized formation energy into species_definitions.

        Parameters:
        -----------
        filename: The filename of relative energies data.
        """
        # Read relative energy data file.
        if os.path.exists(filename):
            globs, locs = {}, {}
            exec(open(filename, "rb").read(), globs, locs)
        else:
            raise IOError("{} is not found.".format(filename))

        # Check data shape.
        self.__chk_data_validity(locs)

        # Get relative energies and pass it to model.
        relative_energies = self.__get_relative_energies(locs)
        self._owner._relative_energies = relative_energies

        # Set flags in model.
        self._owner._has_relative_energy = True

        return
        # }}}


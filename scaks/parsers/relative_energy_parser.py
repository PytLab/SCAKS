import os
import logging

from scipy.linalg import solve

from .parser_base import *
from ..functions import *


class RelativeEnergyParser(ParserBase):
    '''
    Parser for parsing relative energy data

    :param owner: The kinetic model that own this parser
    :type owner: KineticModel
    '''
    def __init__(self, owner):
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

    def parse_data(self, filename="./rel_energy.py", energy_data=None):
        """ Parse relative energy data and model information to add essential
        data to its owner, the kinetic model

        :param filename: The filename of relative energies data.
        :type filename: str

        :param energy_data: Relative energy data
        :type energy_data: dict
        """
        # {{{
        # Read relative energy data file.
        if os.path.exists(filename):
            globs, energy_data = {}, {}
            exec(open(filename, "rb").read(), globs, energy_data)
        elif energy_data is None:
            raise IOError("{} is not found.".format(filename))

        # Check data shape.
        self.__chk_data_validity(energy_data)

        # Get relative energies and pass it to model.
        relative_energies = self.__get_relative_energies(energy_data)
        self._owner._relative_energies = relative_energies

        # Set flags in model.
        self._owner._has_relative_energy = True

        return
        # }}}


import csv
import logging
import os

import kynetix.descriptors.descriptors as dc
from kynetix.parsers.parser_base import ParserBase
from kynetix.errors.error import *
from kynetix.functions import *


class AbsoluteEnergyParser(ParserBase):
    def __init__(self, owner):
        """
        Kinetic Model parser to parse absolute energies.

        Parameter:
        ----------
        owner: The KineticModel object.
        """
        super(AbsoluteEnergyParser, self).__init__(owner)

        # Set tools logger as child of model's
        self.__logger = logging.getLogger('model.parser.AbsoluteEnergyParser')

    def parse_data(self, filename="./abs_energy.py"):
        """
        Read data in absolute energy file.

        Parameters:
        -----------
        filename: file name of energy data.
        """
        globs, locs = {}, {}
        execfile(filename, globs, locs)

        self._owenr._absolute_energy = locs["absolute_energies"]

        # Set flag.
        self._owner._has_absolute_energy = True

        # Get relative energies from absolute energies.
        relative_energies = self.get_relative_from_absolute()
        self._owner._relative_energies = relative_energies

        # Set flag.
        self._owner._has_relative_energy = True

        return


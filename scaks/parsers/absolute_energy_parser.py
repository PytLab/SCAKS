import csv
import logging
import os

from .parser_base import ParserBase
from ..errors.error import *
from ..functions import *


class AbsoluteEnergyParser(ParserBase):
    """ Kinetic Model parser to parse absolute energies.  

    :param owner: The kinetic model that own this parser
    :type owner: KineticModel
    """
    def __init__(self, owner):
        super(AbsoluteEnergyParser, self).__init__(owner)

        # Set tools logger as child of model's
        self.__logger = logging.getLogger('model.parser.AbsoluteEnergyParser')

    def parse_data(self, filename="./abs_energy.py"):
        """ Read data in absolute energy file.

        :param filename: file name of energy data
        :type filename: str

        :return: :obj:`None`
        """
        globs, locs = {}, {}
        exec(open(filename, "rb").read(), globs, locs)

        self._owner._absolute_energies = locs["absolute_energies"]

        # Set flag.
        self._owner._has_absolute_energy = True

        # Get relative energies from absolute energies.
        relative_energies = self._get_relative_from_absolute()
        self._owner._relative_energies = relative_energies

        # Set flag.
        self._owner._has_relative_energy = True

        return


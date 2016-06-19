import csv
import logging
import os

from kynetix.parsers.parser_base import ParserBase
from kynetix.errors.error import *
from kynetix.functions import *


class CsvParser(ParserBase):
    def __init__(self, owner):
        """
        Kinetic Model parser to parse data in csv file.

        Parameter:
        ----------
        owner: The KineticModel object.
        """
        super(CsvParser, self).__init__(owner)

        # Set tools logger as child of model's
        self.__logger = logging.getLogger('model.parser.CsvParser')

    def parse_data(self, filename="./energy.csv", relative=False):
        # {{{
        """
        Read data in csv data file and update the species definitions.

        Parameters:
        -----------
        filename: file name of energy data.

        relative: A useless parameter for compatibility
                  with other parser_data method,
                  so just IGNORE it.

        Returns:
        --------
        species_definitions: The updated species definition of model.
        """
        # NOTE: Get the REFERENCE of model's species definitions.
        attribute_name = mangled_name(self._owner, "species_definitions")
        species_definitions = getattr(self._owner, attribute_name)

        # Check file existance.
        if not os.path.exists(filename):
            msg = "'{}' is not found.".format(filename)
            raise FilesError(msg)

        # Open data file.
        csvfile = open(filename, 'rU')
        reader = csv.DictReader(csvfile)

        # Loop over all data in file to update species definitions.
        for line_dict in reader:
            species_type = line_dict["species_type"]
            species_name = line_dict["species_name"]
            DFT_energy = float(line_dict["DFT_energy"])
            formation_energy = float(line_dict["formation_energy"])
            information = line_dict["information"]

            # Add to species_definition
            if not species_name in species_definitions:
                species_definitions.setdefault(species_name, {})
                for key in ['DFT_energy', 'formation_energy', 'information']:
                    species_definitions[species_name].setdefault(key, eval(key))
            else:
                for key in ['DFT_energy', 'formation_energy', 'information']:
                    species_definitions[species_name][key] = eval(key)

        # Close file.
        csvfile.close()

        # Set flag.
        attribute_name = mangled_name(self._owner, "has_absolute_energy")
        setattr(self._owner, attribute_name, True)

        # Get relative energies from absolute energies.
        relative_energies = self.get_relative_from_absolute()
        attribute_name = mangled_name(self._owner, "relative_energies")
        setattr(self._owner, attribute_name, relative_energies)

        # Set flag.
        attribute_name = mangled_name(self._owner, "has_relative_energy")
        setattr(self._owner, attribute_name, True)

        return
        # }}}

    def has_relative_energy(self):
        """
        Query function for relative energy flag.
        """
        return self.__has_relative_energy

    def has_absolute_energy(self):
        """
        Query function for absolute energy flag.
        """
        return self.__has_absolute_energy


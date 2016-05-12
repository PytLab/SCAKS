import csv
import logging
import os

from kynetix.parsers.parser_base import ParserBase
from kynetix.errors.error import *


class CsvParser(ParserBase):
    def __init__(self, owner, filename="./energy.csv"):
        """
        Kinetic Model parser to parse data in csv file.

        Parameter:
        ----------
        owner: The KineticModel object.
        """
        ParserBase.__init__(self, owner)

        # Set filename.
        self.__filename = filename

        # Set tools logger as child of model's
        self.__logger = logging.getLogger('model.parsers.CsvParser')

    def parse_data(self, relative=False):
        """
        Read data in csv data file and update the species definitions.

        Parameters:
        -----------
        relative: A useless parameter for compatibility
                  with other parser_data method,
                  so just IGNORE it.

        Returns:
        --------
        species_definitions: The updated species definition of model.
        """
        # Get the COPY of model's species.
        species_definitions = self._owner.species_definitions()

        # Check file existance.
        if not os.path.exists(self.__filename):
            msg = "'{}' is not found.".format(self.__filename)
            raise FilesError(msg)

        # Open data file.
        csvfile = open(self.__filename, 'rU')
        reader = csv.DictReader(csvfile)

        # Loop over all data in file to update species definitions.
        for line_dict in reader:
            site_name = line_dict["site_name"]
            species_name = line_dict["species_name"]
            information = line_dict["information"]
            DFT_energy = float(line_dict["DFT_energy"])
            formation_energy = float(line_dict["formation_energy"])

            # Get full name.
            # Gas.
            if site_name == 'gas':
                full_name = species_name + '_g'
            # Liquid.
            elif site_name == 'liquid':
                full_name = species_name + '_l'
            # Species on site.
            else:
                # Get site symbol of current line e.g. 's' or 'ss'.
                for site_symbol in self._owner.site_names():
                    if species_definitions[site_symbol]['site_name'] == site_name:
                        break
                if species_name == 'slab':
                    full_name = site_symbol
                else:
                    full_name = species_name + '_' + site_symbol

            # Check fullname.
            if full_name not in species_definitions:
                msg = (("'{}' is not in model's species definition, " +
                        "it will be igonred").format(full_name))
                self.__logger.warning(msg)
                continue

            # Update species_definition
            for key in ['DFT_energy', 'formation_energy', 'information']:
                species_definitions[full_name].setdefault(key, eval(key))

        # Close file.
        csvfile.close()

        return species_definitions


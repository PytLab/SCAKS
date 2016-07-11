import csv
import logging
import os
import string

from kynetix import mpi_master
from kynetix.functions import string2symbols
from kynetix.parsers.rxn_parser import *
from kynetix.table_makers.table_maker_base import *


class CsvMaker(TableMakerBase):

    # Class variables.
    fieldnames = ["species_type", "species_name",
                  "DFT_energy", "formation_energy",
                  "frequencies", "information"]

    def __init__(self, owner, filename="energy.csv"):
        super(CsvMaker, self).__init__(owner)

        # Set energy filename.
        self.__filename = filename

        # Set logger.
        self.__logger = logging.getLogger('model.table_maker.CsvMaker')

    def init_table(self, filename=None):
        """
        Function to initialize energy table.
        """
        # Set filename.
        if filename is None:
            filename = "init_" + self.__filename

        with open(filename, "wb") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=CsvMaker.fieldnames)

            # Write header.
            writer.writeheader()

            # Gas species.
            for gas in self._owner.gas_names():
                row = {"species_type": "gas",
                       "species_name": gas,
                       "DFT_energy": "0.0",
                       "formation_energy": "0.0",
                       "frequencies": [],
                       "information": "None"}
                writer.writerow(row)

            # Liquid names.
            for liquid in self._owner.liquid_names():
                row = {"species_type": "liquid",
                       "species_name": liquid,
                       "DFT_energy": "0.0",
                       "formation_energy": "0.0",
                       "frequencies": [],
                       "information": "None"}
                writer.writerow(row)

            # Intermediates.
            for intermediate in self._owner.adsorbate_names():
                row = {"species_type": "intermediate",
                       "species_name": intermediate,
                       "DFT_energy": "0.0",
                       "formation_energy": "0.0",
                       "frequencies": [],
                       "information": "None"}
                writer.writerow(row)

            # Transition states.
            for ts in self._owner.transition_state_names():
                row = {"species_type": "transition state",
                       "species_name": ts,
                       "DFT_energy": "0.0",
                       "formation_energy": "0.0",
                       "frequencies": [],
                       "information": "None"}
                writer.writerow(row)

            # Slab.
            for slab in self._owner.site_names():
                row = {"species_type": "slab",
                       "species_name": slab,
                       "DFT_energy": "0.0",
                       "formation_energy": "0.0",
                       "frequencies": [],
                       "information": "None"}
                writer.writerow(row)

        if mpi_master:
            self.__logger.info("Initialize data table - '{}'".format(filename))

    def __get_formation_energy(self, species_name, raw_energy):
        """
        Private function to get generalized formation energy of a species.
        """
        ref_energies = self._owner.ref_energies()
        energy = raw_energy

        # Single site.
        if species_name in self._owner.site_names():
            energy -= ref_energies[species_name]
            return energy

        # Get formula object.
        formula = ChemFormula(species_name)
        element_dict = formula.get_elements_dict()
        site_dict = formula.get_sites_dict()

        # Elements.
        for element, n in element_dict.iteritems():
            energy -= float(n)*float(ref_energies[element])

        # Sites.
        for site, n in site_dict.iteritems():
            # Ignore gas and liquid.
            if site in ("g", "l"):
                continue

            energy -= float(n)*float(ref_energies[site])

        return energy

    def update_table(self, infile=None, outfile=None, remove=True):
        """
        Function to update initial table.

        Parameters:
        -----------
        infile: The name of file read in, str.

        outfile: The name of output file, str.

        remove: Remove the infile or not, bool.
        """
        # Check infile existance.
        if infile is None:
            infile = "init_" + self.__filename

        if not os.path.exists(infile):
            raise IOError("{} not found in current path.".format(infile))

        # Set outfile.
        if outfile is None:
            outfile = self.__filename

        # Reader.
        csvin = open(infile)
        reader = csv.DictReader(csvin)

        # Writer.
        csvout = open(outfile, "wb")
        writer = csv.DictWriter(csvout, fieldnames=CsvMaker.fieldnames)

        fieldnames = reader.fieldnames
        writer.writeheader()

        for row in reader:
            species = row["species_name"]
            raw_energy = float(row["DFT_energy"])
            row["formation_energy"] = self.__get_formation_energy(species, raw_energy)
            writer.writerow(row)

        # Close files.
        csvin.close()
        csvout.close()

        # Remove initial table.
        if remove:
            os.remove(infile)


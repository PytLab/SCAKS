import logging
import os
from math import exp
from operator import mul

import numpy as np

# KMCLibX.
try:
    from KMCLib import *
except ImportError:
    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    print("!!!                                                   !!!")
    print("!!!          WARNING: KMCLib is not installed         !!!")
    print("!!! Any kMC calculation using KMCLib will be disabled !!!")
    print("!!!                                                   !!!")
    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

from ..compatutil import reduce
from ..database.thermo_data import kB_eV
from ..errors.error import *
from ..database.lattice_data import grid_neighbor_offsets
from ..functions import mangled_name
from .rxn_parser import *
from .relative_energy_parser import RelativeEnergyParser
from ..utilities.check_utilities import *


class KMCParser(RelativeEnergyParser):
    """ Parser class for KMC simulation.

    :param owner: The kinetic model that own this parser
    :type owner: KineticModel
    """
    def __init__(self, owner):
        super(KMCParser, self).__init__(owner)

        # Set logger.
        self.__logger = logging.getLogger('model.parsers.KMCParser')

    def parse_data(self,
                   energy_file="./rel_energy.py",
                   processes_file=None,
                   configuration_file=None,
                   sitesmap_file=None):
        """
        Function to get data from kMC input files, processes, configuration & sitesmap.

        :param processes_file: The name of processes definition file,  default name is "kmc_processes.py".
        :type process_file: str

        :param configuration_file: The name of configuration definition file, default name is "kmc_processes.py".
        :type configuration_file: str

        sitesmap_file: The name of sitesmap definition file, default name is "kmc_processes.py".
        :type sitesmap_file: str
        """
        # Basic parsing.
        super(KMCParser, self).parse_data(energy_file)

        # kMC parsing.
        process_dicts = self.parse_processes(filename=processes_file)
        configuration = self.parse_configuration(filename=configuration_file)
        sitesmap = self.construct_sitesmap(filename=sitesmap_file)

        # Pass data to model.
        for var in ["process_dicts", "configuration", "sitesmap"]:
            model_var = mangled_name(self._owner, var)
            setattr(self._owner, model_var, locals()[var])

    def construct_sitesmap(self, filename=None):
        """
        Function to read kmc_site file and create KMCLibSitesMap objects.

        :param filename: The name of configuration file
        :type filename: str

        :return: A KMCSitesMap objects.
        :rtype: :obj:`KMCSitesMap`
        """
        # {{{
        # Load data.
        if filename is None:
            filename = "kmc_sites.py"

        possible_site_types = self._owner.possible_site_types

        # Get site number.
        repetitions = self._owner.repetitions
        basis_sites = self._owner.basis_sites
        nsite = reduce(mul, repetitions)*len(basis_sites)

        def init_default_types():
            """
            Nested function to get default site types according to
            lattice repetitions and basis sites.
            """
            default_type = possible_site_types[0]
            return [default_type]*nsite

        # Get site types.
        if not os.path.exists(filename):
            site_types = init_default_types()
        else:
            globs, locs = {}, {}
            exec(open(filename, "rb").read(), globs, locs)

            if "site_types" not in locs:
                site_types = init_default_types()
            else:
                site_types = locs["site_types"]

        # Check length of site types.
        if len(site_types) != nsite:
            msg = "'site_types' must have {} elements.".format(nsite)
            raise SetupError(msg)

        # Check element type.
        for site_type in site_types:
            if site_type not in possible_site_types:
                msg = "Element '{}' not in possible_site_types '{}'."
                msg = msg.format(site_type, possible_site_types)
                raise SetupError(msg)

        # Construct lattice.
        lattice = self.construct_lattice()

        # Construct sitemap.
        sitesmap = KMCSitesMap(lattice=lattice,
                               types=site_types,
                               possible_types=possible_site_types)

        return sitesmap
        # }}}

    def construct_lattice(self):
        """ Function to construct KMCLattice object.

        :return: A kMC lattice grid
        :rtype: :obj:`KMCLattice`
        """
        # Construct unitcell.
        cell_vectors = np.array(self._owner.cell_vectors)
        basis_sites = np.array(self._owner.basis_sites)
        unit_cell = KMCUnitCell(cell_vectors=cell_vectors, basis_points=basis_sites)

        # Construct lattice.
        repetitions = self._owner.repetitions
        periodic = self._owner.periodic
        lattice = KMCLattice(unit_cell=unit_cell,
                             repetitions=repetitions,
                             periodic=periodic)

        return lattice

    def parse_configuration(self, filename=None):
        """
        Function to read configuration file and create KMCLibConfiguration objects.

        :param filename: The name of configuration file
        :type filename: str

        :return: A kMC configuration
        :rtype: :obj:`KMCConfiguration`
        """
        # {{{
        # Inner function to initialize emtpy lattice.
        def init_empty_types():
            repetitions = self._owner.repetitions
            basis_sites = self._owner.basis_sites
            nsite = reduce(mul, repetitions)*len(basis_sites)
            types = [empty_type for i in range(nsite)]
            return types

        possible_element_types = self._owner.possible_element_types

        # Check.
        empty_type = self._owner.empty_type
        check_string(empty_type, possible_element_types, "empty_type")

        # Load types data in file.
        if filename is None:
            filename = "kmc_configuration.py"

        # Use data in file.
        if os.path.exists(filename):
            globs, locs = {}, {}
            exec(open(filename, "rb").read(), globs, locs)
            if "types" in locs:
                types = locs["types"]
            else:
                types = init_empty_types()
        # Initialize as emtpy lattice.
        else:
            types = init_empty_types()

        # Construct lattice.
        lattice = self.construct_lattice()

        # Instantialize KMCLattice object.
        configuration = KMCConfiguration(lattice=lattice,
                                         types=types,
                                         possible_types=possible_element_types)

        return configuration
        # }}}

    def parse_processes(self, filename=None):
        """
        Function to read processes file and get process dicts.

        :param filename: The name of processes file
        :type filename: str

        :return: Elementary processes for kMC
        :rtype: list of dict
        """
        if filename is None:
            filename = "kmc_processes.py"

        globs, locs = {}, {}
        exec(open(filename, "rb").read(), globs, locs)

        # Get all possible process objects.
        if self._owner.log_allowed:
            msg = "Total {} processes dicts read in.".format(len(locs["processes"]))
            self.__logger.info(msg)

        return locs["processes"]


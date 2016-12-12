import logging
import os
from math import exp
from operator import mul

import numpy as np

# KMCLibX.
try:
    from KMCLib import *
except ImportError:
    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    print "!!!                                                   !!!"
    print "!!!          WARNING: KMCLib is not installed         !!!"
    print "!!! Any kMC calculation using KMCLib will be disabled !!!"
    print "!!!                                                   !!!"
    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

from kynetix.database.thermo_data import kB_eV
from kynetix.errors.error import *
from kynetix.database.lattice_data import grid_neighbor_offsets
from kynetix.functions import mangled_name
from kynetix.parsers.rxn_parser import *
from kynetix.parsers.relative_energy_parser import RelativeEnergyParser
from kynetix.solvers.solver_base import SolverBase
from kynetix.utilities.check_utilities import *


class KMCParser(RelativeEnergyParser):
    def __init__(self, owner):
        """
        Parser class for KMC simulation.
        """
        super(KMCParser, self).__init__(owner)

        # Set logger.
        self.__logger = logging.getLogger('model.parsers.KMCParser')

    def parse_data(self,
                   relative=False,
                   energy_file="./rel_energy.py",
                   processes_file=None,
                   configuration_file=None,
                   sitesmap_file=None):
        """
        Function to get data from kMC input files, processes, configuration & sitesmap.

        Parameters:
        -----------
        processes_file: The name of processes definition file, str.
                        the default name is "kmc_processes.py".

        configuration_file: The name of configuration definition file, str.
                            the default name is "kmc_processes.py".

        sitesmap_file: The name of sitesmap definition file, str.
                       the default name is "kmc_processes.py".
        """
        # Basic parsing.
        super(KMCParser, self).parse_data(relative, energy_file)

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

        Parameters:
        -----------
        filename: The name of configuration file, str.

        Returns:
        --------
        A KMCLibSitesMap objects.
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
            execfile(filename, globs, locs)

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
        """
        Function to construct KMCLattice object.
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

        Parameters:
        -----------
        lattice: The lattice of the configurartion as a KMCLattice.

        filename: The name of configuration file, str.

        Returns:
        --------
        A KMCLibConfiguration objects.
        """
        # {{{
        # Inner function to initialize emtpy lattice.
        def init_empty_types():
            repetitions = self._owner.repetitions
            basis_sites = self._owner.basis_sites
            nsite = reduce(mul, repetitions)*len(basis_sites)
            types = [empty_type for i in xrange(nsite)]
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
            execfile(filename, globs, locs)
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

        Parameters:
        -----------
        filename: The name of processes file, str.

        Returns:
        --------
        A list of process dicts.
        """
        if filename is None:
            filename = "kmc_processes.py"

        globs, locs = {}, {}
        execfile(filename, globs, locs)

        # Get all possible process objects.
        if self._owner.log_allowed:
            msg = "Total {} processes dicts read in.".format(len(locs["processes"]))
            self.__logger.info(msg)

        return locs["processes"]


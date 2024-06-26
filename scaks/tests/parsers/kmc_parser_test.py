import logging
import os
import unittest

from KMCLib import *

from ...models.kmc_model import KMCModel
from ...functions import *
from ...parsers import *

from .. import *


class KMCParserTest(unittest.TestCase):

    def setUp(self):
        self.maxDiff = None
        self.setup_dict = dict(
            rxn_expressions = [
                'CO_g + *_t -> CO_t',
                'CO_g + *_b -> CO_b',
                'O2_g + 2*_b -> 2O_b',
                'CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b',
                'CO_b + *_t <-> CO_t + *_b -> CO_b + *_t',
            ],

            species_definitions = {
                'CO_g': {'pressure': 0.01},
                'O2_g': {'pressure': 0.2},
                'CO2_g': {'pressure': 0.01},
                '*_b': {'site_name': 'bridge', 'type': 'site', 'total': 0.5},
                '*_t': {'site_name': 'top', 'type': 'site', 'total': 0.5},
            },

            temperature = 298.,
            parser = "KMCParser",
            corrector = "ThermodynamicCorrector",
            cell_vectors = [[3.0, 0.0, 0.0],
                            [0.0, 3.0, 0.0],
                            [0.0, 0.0, 3.0]],
            basis_sites = [[0.0, 0.0, 0.0],
                           [0.5, 0.0, 0.0],
                           [0.0, 0.5, 0.0],
                           [0.5, 0.5, 0.0]],
            unitcell_area = 9.0e-20,
            active_ratio = 4./9,
            repetitions = (3, 3, 1),
            periodic = (True, True, False),
            possible_element_types = ["O", "V", "O_s", "C"],
            empty_type = "V",
            possible_site_types = ["P"],
            nstep = 50000,
            random_seed = 13996,
            random_generator = 'MT',
        )

    def test_kmc_parser_construction(self):
        " Test kmc parser can be constructed correctly. "
        # Construction.
        model = KMCModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        parser = model.parser

        # Check the parser class and base class type.
        self.assertTrue(isinstance(parser, KMCParser))
        self.assertEqual(parser.__class__.__base__.__name__, "RelativeEnergyParser")


    def test_parse_processes(self):
        " Make sure we can parse all processes in kmc_processes.py correctly. "
        # Construction.
        model = KMCModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        parser = model.parser
        p = parser.parse_processes(filename=kmc_processes)

        self.assertEqual(11, len(p))
        self.assertTrue(isinstance(p, list))
        self.assertTrue(isinstance(p[0], dict))

    def test_construct_lattice(self):
        " Test we can construct lattice object correctly. "
        # {{{
        model = KMCModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        parser = model.parser
        lattice = parser.construct_lattice()

        # Check.
        self.assertTrue(isinstance(lattice, KMCLattice))

        ref_basis = [[0.0, 0.0, 0.0],
                     [0.5, 0.0, 0.0],
                     [0.0, 0.5, 0.0],
                     [0.5, 0.5, 0.0]]
        ret_basis = lattice.basis().tolist()
        self.assertListEqual(ref_basis, ret_basis)

        ref_periodic = (True, True, False)
        ret_periodic = lattice.periodic()
        self.assertTupleEqual(ref_periodic, ret_periodic)

        ref_repetitions = (3, 3, 1)
        ret_repetitions = lattice.repetitions()
        self.assertTupleEqual(ref_repetitions, ret_repetitions)
        ref_sites=[[0.000000, 0.000000, 0.000000],
                   [0.500000, 0.000000, 0.000000],
                   [0.000000, 0.500000, 0.000000],
                   [0.500000, 0.500000, 0.000000],
                   [0.000000, 1.000000, 0.000000],
                   [0.500000, 1.000000, 0.000000],
                   [0.000000, 1.500000, 0.000000],
                   [0.500000, 1.500000, 0.000000],
                   [0.000000, 2.000000, 0.000000],
                   [0.500000, 2.000000, 0.000000],
                   [0.000000, 2.500000, 0.000000],
                   [0.500000, 2.500000, 0.000000],
                   [1.000000, 0.000000, 0.000000],
                   [1.500000, 0.000000, 0.000000],
                   [1.000000, 0.500000, 0.000000],
                   [1.500000, 0.500000, 0.000000],
                   [1.000000, 1.000000, 0.000000],
                   [1.500000, 1.000000, 0.000000],
                   [1.000000, 1.500000, 0.000000],
                   [1.500000, 1.500000, 0.000000],
                   [1.000000, 2.000000, 0.000000],
                   [1.500000, 2.000000, 0.000000],
                   [1.000000, 2.500000, 0.000000],
                   [1.500000, 2.500000, 0.000000],
                   [2.000000, 0.000000, 0.000000],
                   [2.500000, 0.000000, 0.000000],
                   [2.000000, 0.500000, 0.000000],
                   [2.500000, 0.500000, 0.000000],
                   [2.000000, 1.000000, 0.000000],
                   [2.500000, 1.000000, 0.000000],
                   [2.000000, 1.500000, 0.000000],
                   [2.500000, 1.500000, 0.000000],
                   [2.000000, 2.000000, 0.000000],
                   [2.500000, 2.000000, 0.000000],
                   [2.000000, 2.500000, 0.000000],
                   [2.500000, 2.500000, 0.000000]]
        ret_sites = lattice.sites().tolist()
        self.assertListEqual(ref_sites, ret_sites)

        unitcell = lattice.unitCell()
        self.assertTrue(isinstance(unitcell, KMCUnitCell))
        # }}}

    def test_parse_configuration(self):
        " Make sure we can parse the configuration correctly. "
        # {{{
        model = KMCModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        parser = model.parser
        config = parser.parse_configuration(filename=kmc_config)

        # Check types.
        ref_types = ["V"]*36
        ret_types = config.types()
        self.assertListEqual(ref_types, ret_types)

        # Check atom_id_coordinates.
        ref_atom_id_coords=[[0.000000, 0.000000, 0.000000],
                            [0.500000, 0.000000, 0.000000],
                            [0.000000, 0.500000, 0.000000],
                            [0.500000, 0.500000, 0.000000],
                            [0.000000, 1.000000, 0.000000],
                            [0.500000, 1.000000, 0.000000],
                            [0.000000, 1.500000, 0.000000],
                            [0.500000, 1.500000, 0.000000],
                            [0.000000, 2.000000, 0.000000],
                            [0.500000, 2.000000, 0.000000],
                            [0.000000, 2.500000, 0.000000],
                            [0.500000, 2.500000, 0.000000],
                            [1.000000, 0.000000, 0.000000],
                            [1.500000, 0.000000, 0.000000],
                            [1.000000, 0.500000, 0.000000],
                            [1.500000, 0.500000, 0.000000],
                            [1.000000, 1.000000, 0.000000],
                            [1.500000, 1.000000, 0.000000],
                            [1.000000, 1.500000, 0.000000],
                            [1.500000, 1.500000, 0.000000],
                            [1.000000, 2.000000, 0.000000],
                            [1.500000, 2.000000, 0.000000],
                            [1.000000, 2.500000, 0.000000],
                            [1.500000, 2.500000, 0.000000],
                            [2.000000, 0.000000, 0.000000],
                            [2.500000, 0.000000, 0.000000],
                            [2.000000, 0.500000, 0.000000],
                            [2.500000, 0.500000, 0.000000],
                            [2.000000, 1.000000, 0.000000],
                            [2.500000, 1.000000, 0.000000],
                            [2.000000, 1.500000, 0.000000],
                            [2.500000, 1.500000, 0.000000],
                            [2.000000, 2.000000, 0.000000],
                            [2.500000, 2.000000, 0.000000],
                            [2.000000, 2.500000, 0.000000],
                            [2.500000, 2.500000, 0.000000]]
        ret_atom_id_coords = config.atomIDCoordinates().tolist()
        self.assertListEqual(ref_atom_id_coords, ret_atom_id_coords)

        # Check atom_id_types, should be the same as types.
        ref_atom_id_types = tuple(["V"]*36)
        ret_atom_id_types = config.atomIDTypes()
        self.assertTupleEqual(ref_atom_id_types, ret_atom_id_types)
        # }}}

    def test_construct_sitesmap(self):
        " Make sure we can construct sitesmap correctly. "
        # {{{
        model = KMCModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        parser = model.parser
        sitesmap = parser.construct_sitesmap(filename=kmc_sites)

        # Check sites.
        ref_sites=[[0.000000, 0.000000, 0.000000],
                   [0.500000, 0.000000, 0.000000],
                   [0.000000, 0.500000, 0.000000],
                   [0.500000, 0.500000, 0.000000],
                   [0.000000, 1.000000, 0.000000],
                   [0.500000, 1.000000, 0.000000],
                   [0.000000, 1.500000, 0.000000],
                   [0.500000, 1.500000, 0.000000],
                   [0.000000, 2.000000, 0.000000],
                   [0.500000, 2.000000, 0.000000],
                   [0.000000, 2.500000, 0.000000],
                   [0.500000, 2.500000, 0.000000],
                   [1.000000, 0.000000, 0.000000],
                   [1.500000, 0.000000, 0.000000],
                   [1.000000, 0.500000, 0.000000],
                   [1.500000, 0.500000, 0.000000],
                   [1.000000, 1.000000, 0.000000],
                   [1.500000, 1.000000, 0.000000],
                   [1.000000, 1.500000, 0.000000],
                   [1.500000, 1.500000, 0.000000],
                   [1.000000, 2.000000, 0.000000],
                   [1.500000, 2.000000, 0.000000],
                   [1.000000, 2.500000, 0.000000],
                   [1.500000, 2.500000, 0.000000],
                   [2.000000, 0.000000, 0.000000],
                   [2.500000, 0.000000, 0.000000],
                   [2.000000, 0.500000, 0.000000],
                   [2.500000, 0.500000, 0.000000],
                   [2.000000, 1.000000, 0.000000],
                   [2.500000, 1.000000, 0.000000],
                   [2.000000, 1.500000, 0.000000],
                   [2.500000, 1.500000, 0.000000],
                   [2.000000, 2.000000, 0.000000],
                   [2.500000, 2.000000, 0.000000],
                   [2.000000, 2.500000, 0.000000],
                   [2.500000, 2.500000, 0.000000]]
        ret_sites = sitesmap.sites().tolist()
        self.assertListEqual(ref_sites, ret_sites)

        # Check repetitions.
        ref_repetitions = (3, 3, 1)
        ret_repetitions = sitesmap.cellRepetitions()
        self.assertTupleEqual(ref_repetitions, ret_repetitions)

        # Check site types.
        ref_types = ["P"]*36
        ret_types = sitesmap.types()
        self.assertListEqual(ref_types, ret_types)

        # Check possible site types.
        ref_possible_types = {"*": 0, "P": 1}
        ret_possible_types = sitesmap.possibleTypes()
        self.assertDictEqual(ref_possible_types, ret_possible_types)

        # Check site type mapping.
        ref_types = [1]*36
        ret_types = sitesmap.siteTypesMapping(sitesmap.types())
        self.assertListEqual(ref_types, ret_types)
        # }}}

    def test_parse_data(self):
        " Make sure the kMC data can be parsed in correctly. "
        # Construction.
        model = KMCModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        model.parser.parse_data(energy_file=kmc_energy,
                                processes_file=kmc_processes,
                                configuration_file=kmc_config,
                                sitesmap_file=kmc_sites)

        self.assertTrue(model.has_relative_energy)
        self.assertTrue(hasattr(model, "_KMCModel__process_dicts"))
        self.assertTrue(hasattr(model, "_KMCModel__configuration"))
        self.assertTrue(hasattr(model, "_KMCModel__sitesmap"))

    def tearDown(self):
        cleanup()

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(KMCParserTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


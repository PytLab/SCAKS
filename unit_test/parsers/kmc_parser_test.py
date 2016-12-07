import logging
import os
import unittest

from KMCLib import *

from kynetix.models.kmc_model import KMCModel
from kynetix.functions import *
from kynetix.parsers import *

from unit_test import *


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
                'b': {'site_name': 'bridge', 'type': 'site', 'total': 0.5},
                't': {'site_name': 'top', 'type': 'site', 'total': 0.5},
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
        model = KMCModel(setup_dict=self.setup_dict, verbosity=logging.WARNING)
        parser = model.parser

        # Check the parser class and base class type.
        self.assertTrue(isinstance(parser, KMCParser))
        self.assertEqual(parser.__class__.__base__.__name__, "RelativeEnergyParser")

    def test_get_relative_energies(self):
        " Make sure we can get correct relative energies. "
        # Construction.
        model = KMCModel(setup_dict=self.setup_dict, verbosity=logging.WARNING)
        parser = model.parser
        parser.parse_data(filename=kmc_energy, relative=True)

        ref_e = (0.0, 1.92, -1.92)
        ret_e = parser._KMCParser__get_relative_energies('CO_g + *_t -> CO_t')
        self.assertTupleEqual(ref_e, ret_e)

        ref_e = (0.0, 2.09, -2.09)
        ret_e = parser._KMCParser__get_relative_energies('CO_g + *_b -> CO_b')
        self.assertTupleEqual(ref_e, ret_e)

        ref_e = (0.0, 3.48, -3.48)
        ret_e = parser._KMCParser__get_relative_energies('O2_g + 2*_b -> 2O_b')
        self.assertTupleEqual(ref_e, ret_e)

        ref_e = (0.39, 0.8500000000000001, -0.46)
        ret_e = parser._KMCParser__get_relative_energies('CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b')
        self.assertTupleEqual(ref_e, ret_e)

    def test_get_rxn_rates(self):
        " Make sure we can get correct forward and reverse rates for a reaction. "
        # Construction.
        model = KMCModel(setup_dict=self.setup_dict, verbosity=logging.WARNING)
        parser = model.parser
        parser.parse_data(filename=kmc_energy, relative=True)

        ref_r = (1575287.974387463, 3.8789566422291146e-14)
        ret_r = parser._KMCParser__get_rxn_rates('CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b')
        self.assertTupleEqual(ref_r, ret_r)

        ref_r = (215.85343473385328, 1.7062993852898129e-44)
        ret_r = parser._KMCParser__get_rxn_rates('O2_g + 2*_b -> 2O_b')
        self.assertTupleEqual(ref_r, ret_r)

        ref_r = (11.535554738754854, 1.3130247359797898e-18)
        ret_r = parser._KMCParser__get_rxn_rates('CO_g + *_t -> CO_t')
        self.assertTupleEqual(ref_r, ret_r)

    def test_parse_single_process(self):
        " Make sure we can parse a process dict correctly. "
        # {{{
        # Construction.
        model = KMCModel(setup_dict=self.setup_dict, verbosity=logging.WARNING)
        parser = model.parser
        parser.parse_data(filename=kmc_energy, relative=True)

        process_dict = {"reaction": "CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b",
                        "description": "CO and O couple and desorption.",
                        "coordinates_group":[[[0.0, 0.0, 0.0], [0.5, 0.5, 0.0]],
                                             [[0.0, 0.0, 0.0], [0.5, -0.5, 0.0]],
                                             [[0.0, 0.0, 0.0], [-0.5, 0.5, 0.0]],
                                             [[0.0, 0.0, 0.0], [-0.5, -0.5, 0.0]]],
                        "elements_before": ["V", "V"],
                        "elements_after": ["O_s", "O_s"],
                        "basis_sites": [1, 2],
                        "fast": True}

        processes = parser._KMCParser__parse_single_process(process_dict)

        # Check processes number.
        self.assertEqual(16, len(processes))

        # Check a the first process object.
        p = processes[0]
        self.assertListEqual(p.basisSites(), [1])
        self.assertListEqual(p.elementsBefore(), ["V", "V"])
        self.assertListEqual(p.elementsAfter(), ["O_s", "O_s"])
        self.assertTrue(p.fast())

        # Check coordinates.
        ref_coords = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.0]]
        c = p.localConfigurations()[0]
        ret_coords = c.coordinates().tolist()
        self.assertListEqual(ref_coords, ret_coords)

        # Check a the second process object.
        p = processes[-1]
        self.assertListEqual(p.basisSites(), [2])
        self.assertListEqual(p.elementsBefore(), ["O_s", "O_s"])
        self.assertListEqual(p.elementsAfter(), ["V", "V"])
        self.assertTrue(p.fast())

        # Check coordinates.
        ref_coords = [[0.0, 0.0, 0.0], [-0.5, -0.5, 0.0]]
        c = p.localConfigurations()[0]
        ret_coords = c.coordinates().tolist()
        self.assertListEqual(ref_coords, ret_coords)
        # }}}

    def test_parse_processes(self):
        " Make sure we can parse all processes in kmc_processes.py correctly. "
        # Construction.
        model = KMCModel(setup_dict=self.setup_dict, verbosity=logging.WARNING)
        parser = model.parser
        parser.parse_data(filename=kmc_energy, relative=True)
        p = parser.parse_processes(filename=kmc_processes)

        self.assertEqual(37, len(p))

    def test_construct_lattice(self):
        " Test we can construct lattice object correctly. "
        # {{{
        model = KMCModel(setup_dict=self.setup_dict, verbosity=logging.WARNING)
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
        model = KMCModel(setup_dict=self.setup_dict, verbosity=logging.WARNING)
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
        model = KMCModel(setup_dict=self.setup_dict, verbosity=logging.WARNING)
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

    def test_process_reactions_mapping(self):
        " Test process mapping query function. "
        # {{{
        model = KMCModel(setup_dict=self.setup_dict, verbosity=logging.WARNING)
        parser = model.parser
        parser.parse_data(filename=kmc_energy, relative=True)
        p = parser.parse_processes(filename=kmc_processes)

        ref_mapping = ['CO_g + *_t -> CO_t(->)',
                       'CO_g + *_t -> CO_t(<-)',
                       'CO_g + *_b -> CO_b(->)',
                       'CO_g + *_b -> CO_b(<-)',
                       'CO_g + *_b -> CO_b(->)',
                       'CO_g + *_b -> CO_b(<-)',
                       'O2_g + 2*_b -> 2O_b(->)',
                       'O2_g + 2*_b -> 2O_b(<-)',
                       'O2_g + 2*_b -> 2O_b(->)',
                       'O2_g + 2*_b -> 2O_b(<-)',
                       'O2_g + 2*_b -> 2O_b(->)',
                       'O2_g + 2*_b -> 2O_b(<-)',
                       'O2_g + 2*_b -> 2O_b(->)',
                       'O2_g + 2*_b -> 2O_b(<-)',
                       'CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b(->)',
                       'CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b(<-)',
                       'CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b(->)',
                       'CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b(<-)',
                       'CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b(->)',
                       'CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b(<-)',
                       'CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b(->)',
                       'CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b(<-)',
                       'CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b(->)',
                       'CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b(<-)',
                       'CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b(->)',
                       'CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b(<-)',
                       'CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b(->)',
                       'CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b(<-)',
                       'CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b(->)',
                       'CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b(<-)',]
        ret_mapping = parser.process_mapping()

        self.assertListEqual(ref_mapping, ret_mapping)
        # }}}

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(KMCParserTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


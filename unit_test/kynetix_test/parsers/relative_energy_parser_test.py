import logging
import unittest

from kynetix.model import KineticModel
from kynetix.parsers import *


class RelativeEnergyParserTest(unittest.TestCase):

    def setUp(self):
        self.maxDiff = None

    def test_relative_energy_parser_construction(self):
        " Test relative energy parser can be constructed. "
        # Construction.
        model = KineticModel(setup_file="input_files/relative_energy_parser.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()

        # Check the parser class and base class type.
        self.assertTrue(isinstance(parser, RelativeEnergyParser))
        self.assertEqual(parser.__class__.__base__.__name__, "ParserBase")

        # Check query.
        self.assertFalse(parser.has_relative_energy())
        self.assertFalse(parser.has_absolute_energy())

        self.assertListEqual([0.0, 0.0, 1.25], parser.Ga())
        self.assertListEqual([-0.758, -2.64, 0.324], parser.dG())
        ref_relative_energies = {'Ga': [0.0, 0.0, 1.25], 'dG': [-0.758, -2.64, 0.324]}
        self.assertDictEqual(ref_relative_energies, parser.relative_energies())

    def test_unknown_species(self):
        " Make sure we can get unknown species correctly. "

        # Construction.
        model = KineticModel(setup_file="input_files/relative_energy_parser.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()

        ref_unknown_species = ['O2_g', 'CO_s', 'O_s', 'CO-O_2s']
        ret_unknown_species = parser._RelativeEnergyParser__get_unknown_species()

        # Check.
        self.assertListEqual(ref_unknown_species, ret_unknown_species)

        # Now add 'O2_g' to ref_species.
        model._KineticModel__ref_species.append('O2_g')

        ref_unknown_species = ['CO_s', 'O_s', 'CO-O_2s']
        ret_unknown_species = parser._RelativeEnergyParser__get_unknown_species()

        # Check.
        self.assertListEqual(ref_unknown_species, ret_unknown_species)

    def test_unknown_coeff_vector(self):
        " Make sure we can get unknown species vector and energy value. "
        # Construction.
        model = KineticModel(setup_file="input_files/relative_energy_parser.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()

        rxn_list = [['CO_s', 'O_s'], ['CO-O_2s'], ['CO2_g', '2*_s']]
        ref_vectors = [[0, -1, -1, 1], [0, -1, -1, 0]]
        ref_Ga = 1.25
        ref_dG = 0.324

        ret_vectors, [ret_Ga, ret_dG] = \
            parser._RelativeEnergyParser__get_unknown_coeff_vector(rxn_list)

        self.assertListEqual(ref_vectors, ret_vectors)
        self.assertAlmostEqual(ref_Ga, ret_Ga)
        self.assertAlmostEqual(ref_dG, ret_dG)

        # Reaction without TS.
        rxn_list = [['CO_g', '*_s'], ['CO_s']]
        ref_vectors = [[0, 1, 0, 0]]
        ref_dG = -0.758

        ret_vectors, [ret_dG] = \
            parser._RelativeEnergyParser__get_unknown_coeff_vector(rxn_list)

        self.assertListEqual(ref_vectors, ret_vectors)
        self.assertAlmostEqual(ref_dG, ret_dG)

    def test_data_conversion(self):
        " Test relative energy can be converted to absolute energy. "
        # Construction.
        model = KineticModel(setup_file="input_files/relative_energy_parser.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()

        # Check G dict before conversion.
        self.assertDictEqual({}, parser._RelativeEnergyParser__G_dict)

        # Check after conversion.
        parser._RelativeEnergyParser__convert_data()
        ref_G_dict = {'CO-O_2s': 0.9259999999999999,
                      'CO2_g': 0.0,
                      'CO_g': 0.0,
                      'CO_s': -0.758,
                      'O2_g': 3.508,
                      'O_s': 0.43399999999999994,
                      's': 0.0}
        ret_G_dict = parser._RelativeEnergyParser__G_dict
        self.assertDictEqual(ref_G_dict, ret_G_dict)

    def test_data_parser(self):
        " Test data in relative energy file can be parsed correctly. "

        # Construction.
        model = KineticModel(setup_file="input_files/relative_energy_parser.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()

        # Check before parse.
        ref_species_definitions = {'CO-O_2s': {'elements': {'C': 1, 'O': 2},
                                    'site': 's',
                                    'site_number': 2,
                                    'type': 'transition_state'},
                                   'CO2_g': {'elements': {'C': 1, 'O': 2},
                                    'pressure': 0.0,
                                    'site': 'g',
                                    'site_number': 1,
                                    'type': 'gas'},
                                   'CO_g': {'elements': {'C': 1, 'O': 1},
                                    'pressure': 1.0,
                                    'site': 'g',
                                    'site_number': 1,
                                    'type': 'gas'},
                                   'CO_s': {'elements': {'C': 1, 'O': 1},
                                    'site': 's',
                                    'site_number': 1,
                                    'type': 'adsorbate'},
                                   'O2_g': {'elements': {'O': 2},
                                    'pressure': 0.3333333333333333,
                                    'site': 'g',
                                    'site_number': 1,
                                    'type': 'gas'},
                                   'O_s': {'elements': {'O': 1},
                                    'site': 's',
                                    'site_number': 1,
                                    'type': 'adsorbate'},
                                   's': {'site_name': '111', 'total': 1.0, 'type': 'site'}}
        self.assertDictEqual(ref_species_definitions, model.species_definitions())

        self.assertFalse(parser.has_absolute_energy())
        self.assertFalse(parser.has_relative_energy())

        self.assertTrue(hasattr(parser, "_RelativeEnergyParser__Ga"))
        self.assertTrue(hasattr(parser, "_RelativeEnergyParser__dG"))
        self.assertFalse(hasattr(parser, "_RelativeEnergyParser__Ea"))
        self.assertFalse(hasattr(parser, "_RelativeEnergyParser__dE"))

        # Parse absolute data.
        ret_species_definitions = parser.parse_data(relative=False)
        ref_species_definitions = {'CO-O_2s': {'elements': {'C': 1, 'O': 2},
                                    'formation_energy': 0.9259999999999999,
                                    'site': 's',
                                    'site_number': 2,
                                    'type': 'transition_state'},
                                   'CO2_g': {'elements': {'C': 1, 'O': 2},
                                    'formation_energy': 0.0,
                                    'pressure': 0.0,
                                    'site': 'g',
                                    'site_number': 1,
                                    'type': 'gas'},
                                   'CO_g': {'elements': {'C': 1, 'O': 1},
                                    'formation_energy': 0.0,
                                    'pressure': 1.0,
                                    'site': 'g',
                                    'site_number': 1,
                                    'type': 'gas'},
                                   'CO_s': {'elements': {'C': 1, 'O': 1},
                                    'formation_energy': -0.758,
                                    'site': 's',
                                    'site_number': 1,
                                    'type': 'adsorbate'},
                                   'O2_g': {'elements': {'O': 2},
                                    'formation_energy': 3.508,
                                    'pressure': 0.3333333333333333,
                                    'site': 'g',
                                    'site_number': 1,
                                    'type': 'gas'},
                                   'O_s': {'elements': {'O': 1},
                                    'formation_energy': 0.43399999999999994,
                                    'site': 's',
                                    'site_number': 1,
                                    'type': 'adsorbate'},
                                   's': {'formation_energy': 0.0,
                                    'site_name': '111',
                                    'total': 1.0,
                                    'type': 'site'}}
        self.assertDictEqual(ref_species_definitions, ret_species_definitions)
        self.assertTrue(parser.has_absolute_energy())
        self.assertFalse(parser.has_relative_energy())

        # Check if relative == true.
        # Construct again.
        model = KineticModel(setup_file="input_files/relative_energy_parser.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()

        # Check.
        ret_relative_energies = parser.parse_data(relative=True)
        ref_relative_energies = {'Ga': [0.0, 0.0, 1.25], 'dG': [-0.758, -2.64, 0.324]}
        self.assertDictEqual(ret_relative_energies, ref_relative_energies)
        self.assertFalse(parser.has_absolute_energy())
        self.assertTrue(parser.has_relative_energy())

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(RelativeEnergyParserTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


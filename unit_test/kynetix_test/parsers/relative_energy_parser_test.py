import logging
import os
import unittest

from kynetix.model import KineticModel
from kynetix.functions import *
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

        # Read relative energy data file.
        filename = "input_files/rel_energy.py"
        if os.path.exists(filename):
            globs, locs = {}, {}
            execfile(filename, globs, locs)

            # Set variables in data file as attr of parser
            for key in locs:
                setattr(parser, "_" + key, locs[key])
        else:
            raise IOError("{} is not found.".format(filename))

        rxn_expression = 'CO_s + O_s <-> CO-O_2s -> CO2_g + 2*_s'
        ref_vectors = [[0, -1, -1, 1], [0, -1, -1, 0]]
        ref_Ga = 1.25
        ref_dG = 0.324

        ret_vectors, [ret_Ga, ret_dG] = \
            parser._RelativeEnergyParser__get_unknown_coeff_vector(rxn_expression)

        self.assertListEqual(ref_vectors, ret_vectors)
        self.assertAlmostEqual(ref_Ga, ret_Ga)
        self.assertAlmostEqual(ref_dG, ret_dG)

        # Reaction without TS.
        rxn_expression = 'CO_g + *_s -> CO_s'
        ref_vectors = [[0, 1, 0, 0]]
        ref_dG = -0.758

        ret_vectors, [ret_dG] = \
            parser._RelativeEnergyParser__get_unknown_coeff_vector(rxn_expression)

        self.assertListEqual(ref_vectors, ret_vectors)
        self.assertAlmostEqual(ref_dG, ret_dG)

        # Another one.
        rxn_expression = 'O2_g + 2*_s -> 2O_s'
        ref_vectors = [[-1, 0, 2, 0]]
        ref_dG = -2.64

        ret_vectors, [ret_dG] = \
            parser._RelativeEnergyParser__get_unknown_coeff_vector(rxn_expression)

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

        # Read relative energy data file.
        filename = "input_files/rel_energy.py"
        if os.path.exists(filename):
            globs, locs = {}, {}
            execfile(filename, globs, locs)

            # Set variables in data file as attr of parser
            for key in locs:
                setattr(parser, "_" + key, locs[key])
        else:
            raise IOError("{} is not found.".format(filename))

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

    def test_data_parse(self):
        " Test data in relative energy file can be parsed correctly. "

        # Construction.
        model = KineticModel(setup_file="input_files/relative_energy_parser.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()

        # Check before parse.
        ref_species_definitions = {'CO2_g': {'pressure': 0.0},
                                   'CO_g': {'pressure': 1.0},
                                   'O2_g': {'pressure': 0.3333333333333333},
                                   's': {'site_name': '111', 'total': 1.0, 'type': 'site'}}
        self.assertDictEqual(ref_species_definitions, model.species_definitions())

        self.assertFalse(model.has_absolute_energy())
        self.assertFalse(model.has_relative_energy())

        # Parse absolute data.
        parser.parse_data(relative=False, filename="input_files/rel_energy.py")
        ref_species_definitions = {'CO-O_2s': {'formation_energy': 0.9259999999999999},
                                   'CO2_g': {'formation_energy': 0.0, 'pressure': 0.0},
                                   'CO_g': {'formation_energy': 0.0, 'pressure': 1.0},
                                   'CO_s': {'formation_energy': -0.758},
                                   'O2_g': {'formation_energy': 3.508, 'pressure': 0.3333333333333333},
                                   'O_s': {'formation_energy': 0.43399999999999994},
                                   's': {'formation_energy': 0.0,
                                         'site_name': '111',
                                         'total': 1.0,
                                         'type': 'site'}}

        self.assertTrue(hasattr(parser, "_Ga"))
        self.assertTrue(hasattr(parser, "_dG"))

        self.assertDictEqual(ref_species_definitions, model.species_definitions())
        self.assertTrue(model.has_absolute_energy())
        self.assertTrue(model.has_relative_energy())

        # Check if relative == true.
        # Construct again.
        model = KineticModel(setup_file="input_files/relative_energy_parser.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()

        # Check.
        parser.parse_data(relative=True, filename="input_files/rel_energy.py")
        ref_relative_energies = {'Gaf': [0.0, 0.0, 1.25],
                                 'Gar': [0.758, 2.64, 0.9259999999999999],
                                 'dG': [-0.758, -2.64, 0.324]}
        self.assertDictEqual(ref_relative_energies, model.relative_energies())
        self.assertFalse(model.has_absolute_energy())
        self.assertTrue(model.has_relative_energy())

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(RelativeEnergyParserTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


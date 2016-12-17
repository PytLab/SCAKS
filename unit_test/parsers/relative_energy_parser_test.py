import logging
import os
import unittest

from kynetix.models.micro_kinetic_model import MicroKineticModel
from kynetix.functions import *
from kynetix.parsers import *

from unit_test import *


class RelativeEnergyParserTest(unittest.TestCase):

    def setUp(self):
        self.maxDiff = None
        self.setup_dict = dict(
            rxn_expressions = [
                'CO_g + *_s -> CO_s',
                'O2_g + 2*_s -> 2O_s',
                'CO_s + O_s <-> CO-O_2s -> CO2_g + 2*_s',
            ],

            species_definitions = {
                'CO_g': {'pressure': 1.0},
                'O2_g': {'pressure': 1./3.},
                'CO2_g': {'pressure': 0.00},
                's': {'site_name': '111', 'type': 'site', 'total': 1.0},
            },

            temperature = 450.0,
            parser = "RelativeEnergyParser",
            solver = "SteadyStateSolver",
            corrector = "ThermodynamicCorrector",
            plotter = "EnergyProfilePlotter",
            rootfinding = 'ConstrainedNewton',
            tolerance = 1e-20,
            max_rootfinding_iterations = 100,
        )

    def test_relative_energy_parser_construction(self):
        " Test relative energy parser can be constructed. "
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        parser = model.parser

        # Check the parser class and base class type.
        self.assertTrue(isinstance(parser, RelativeEnergyParser))
        self.assertEqual(parser.__class__.__base__.__name__, "ParserBase")

    def test_data_parse(self):
        " Test data in relative energy file can be parsed correctly. "

        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        parser = model.parser

        # Check before parse.
        ref_species_definitions = {'CO2_g': {'pressure': 0.0},
                                   'CO_g': {'pressure': 1.0},
                                   'O2_g': {'pressure': 0.3333333333333333},
                                   's': {'site_name': '111', 'total': 1.0, 'type': 'site'}}
        self.assertDictEqual(ref_species_definitions, model.species_definitions)

        self.assertFalse(model.has_absolute_energy)
        self.assertFalse(model.has_relative_energy)

        # Parse absolute data.
        parser.parse_data(filename=mkm_energy)

        # Check.
        ref_relative_energies = {'Gaf': [0.0, 0.0, 1.25],
                                 'Gar': [0.758, 2.64, 0.9259999999999999],
                                 'dG': [-0.758, -2.64, 0.324]}
        self.assertDictEqual(ref_relative_energies, model.relative_energies)
        self.assertFalse(model.has_absolute_energy)
        self.assertTrue(model.has_relative_energy)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(RelativeEnergyParserTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


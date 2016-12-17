import logging
import unittest

from kynetix.models.micro_kinetic_model import MicroKineticModel
from kynetix.parsers import *

from unit_test import *


class MicroKineticModelTest(unittest.TestCase):

    def setUp(self):
        # Test case setting.
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
            ref_species = ['CO_g', 'CO2_g', 's']
        )

    def test_mkm_construction_query(self):
        " Test micro kinetic model can be constructed with parser. "
        # Test construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)

        # Load data in setup file.
        self.assertEqual(model.corrector.__class__.__name__, self.setup_dict["corrector"])
        self.assertEqual(model.plotter.__class__.__name__, self.setup_dict["plotter"])
        self.assertListEqual(model.rxn_expressions, self.setup_dict["rxn_expressions"])
        self.assertEqual(model.temperature, self.setup_dict["temperature"])
        self.assertListEqual(model.ref_species, self.setup_dict["ref_species"])
        self.assertEqual(model.logger_level, logging.WARNING)
        self.assertEqual(model.decimal_precision, 100)

        self.assertTrue(isinstance(model.parser, RelativeEnergyParser))

    def test_run(self):
        " Test micro kinetic model can run correctly. "
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        init_cvgs = [0.9, 0.1]
        model.run(data_file=mkm_energy, init_cvgs=init_cvgs)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(MicroKineticModelTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


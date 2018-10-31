import logging
import re
import unittest

import numpy as np
from mpmath import mpf

from catynetics.models.micro_kinetic_model import MicroKineticModel
from catynetics.plotters import *

from unit_test import *


class EnergyProfilePlotterTest(unittest.TestCase):

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
                '*_s': {'site_name': '111', 'type': 'site', 'total': 1.0},
            },

            temperature = 450.0,
            parser = "RelativeEnergyParser",
            solver = "SteadyStateSolver",
            corrector = "ThermodynamicCorrector",
            plotter = "EnergyProfilePlotter",
            rootfinding = 'ConstrainedNewton',
            decimal_precision = 100,
            tolerance = 1e-20,
            max_rootfinding_iterations = 100,
        )

    def test_construction_and_query(self):
        " Test plotter construction and query. "
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        parser = model.parser
        plotter = model.plotter

        self.assertTrue(isinstance(plotter, EnergyProfilePlotter))

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(EnergyProfilePlotterTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


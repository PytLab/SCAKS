import logging
import re
import unittest

import numpy as np
from mpmath import mpf

from kynetix.model import KineticModel
from kynetix.plotters import *


class EnergyProfilePlotterTest(unittest.TestCase):

    def setUp(self):
        # Test case setting.
        self.maxDiff = None

    def test_construction_and_query(self):
        " Test plotter construction and query. "
        # Construction.
        model = KineticModel(setup_file="input_files/energy_profile_plotter.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        plotter = model.plotter()

        self.assertTrue(isinstance(plotter, EnergyProfilePlotter))

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(EnergyProfilePlotterTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


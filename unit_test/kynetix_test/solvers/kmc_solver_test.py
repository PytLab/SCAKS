import logging
import re
import unittest

import numpy as np

from kynetix.model import KineticModel
from kynetix.solvers import *


class KMCSolverTest(unittest.TestCase):

    def setUp(self):
        # Test case setting.
        self.maxDiff = None

    def test_construction(self):
        " Make sure KMCSolver object can be constructed correctly. "
        model = KineticModel(setup_file="kmc_inputs/kmc_solver.mkm",
                             verbosity=logging.WARNING)
        solver = model.solver()

        self.assertTrue(isinstance(solver, KMCSolver))

    def test_get_control_parameter(self):
        " Make sure we can get KMCControlParameter object. "
        model = KineticModel(setup_file="kmc_inputs/kmc_solver.mkm",
                             verbosity=logging.WARNING)
        solver = model.solver()

        control_parameters = solver.get_control_parameters()

        # Check.
        self.assertEqual(10, control_parameters.dumpInterval())
        self.assertEqual(5, control_parameters.analysisInterval())
        self.assertEqual(50000, control_parameters.numberOfSteps())
        self.assertEqual(0, control_parameters.rngType())
        self.assertEqual(13996, control_parameters.seed())
        self.assertEqual(False, control_parameters.timeSeed())

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(KMCSolverTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


import logging
import re
import unittest

import numpy as np

from kynetix.model import KineticModel
from kynetix.solvers import *

from unit_test import *


class KMCSolverTest(unittest.TestCase):

    def setUp(self):
        # Test case setting.
        self.maxDiff = None
        self.setup = kmc_path + "/kmc_solver.mkm"

    def test_construction(self):
        " Make sure KMCSolver object can be constructed correctly. "
        model = KineticModel(setup_file=self.setup, verbosity=logging.WARNING)
        solver = model.solver()

        self.assertTrue(isinstance(solver, KMCSolver))

    def test_get_control_parameter(self):
        " Make sure we can get KMCControlParameter object. "
        model = KineticModel(setup_file=self.setup, verbosity=logging.WARNING)
        solver = model.solver()

        control_parameters = solver.get_control_parameters()

        # Check.
        self.assertEqual(10, control_parameters.dumpInterval())
        self.assertEqual(5, control_parameters.analysisInterval())
        self.assertEqual(50, control_parameters.numberOfSteps())
        self.assertEqual(0, control_parameters.rngType())
        self.assertEqual(13996, control_parameters.seed())
        self.assertEqual(False, control_parameters.timeSeed())

    def test_run(self):
        " Test the we can run the kmc model correctly. "
        model = KineticModel(setup_file=self.setup, verbosity=logging.WARNING)
        parser = model.parser()
        parser.parse_data(filename=kmc_energy, relative=True)
        solver = model.solver()
        model._KineticModel__processes = parser.parse_processes(filename=kmc_processes)
        model._KineticModel__configuration = parser.parse_configuration(filename=kmc_config)
        model._KineticModel__sitesmap = parser.construct_sitesmap(filename=kmc_sites)

        solver.run()

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(KMCSolverTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


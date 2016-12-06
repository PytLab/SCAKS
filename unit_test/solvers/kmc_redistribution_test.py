#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import unittest

from kynetix.model import KineticModel
from kynetix.solvers import *

from unit_test import *


class KMCRedistributionTest(unittest.TestCase):

    def setUp(self):
        # Test case setting.
        self.maxDiff = None

    def test_run_with_split_redistribution(self):
        " Make sure the model can run with split redistribution operation. "
        setup_file = kmc_path + "/kmc_split_redistribution.mkm"
        model = KineticModel(setup_file=setup_file, verbosity=logging.WARNING)
        parser = model.parser
        parser.parse_data(filename=kmc_energy, relative=True)

        # Run the model with redistribution.
        model.run_kmc(processes_file=kmc_processes,
                      configuration_file=kmc_config,
                      sitesmap_file=kmc_sites)

    def test_run_with_process_redistribution(self):
        " Make sure the model can run with process redistribution operation. "
        setup_file = kmc_path + "/kmc_process_redistribution.mkm"
        model = KineticModel(setup_file=setup_file, verbosity=logging.WARNING)
        parser = model.parser
        parser.parse_data(filename=kmc_energy, relative=True)

        # Run the model with redistribution.
        model.run_kmc(processes_file=kmc_processes,
                      configuration_file=kmc_config,
                      sitesmap_file=kmc_sites)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(KMCRedistributionTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


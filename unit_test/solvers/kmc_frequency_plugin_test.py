import logging
import re
import unittest

import numpy as np

from kynetix.models.kmc_model import KMCModel
from kynetix.solvers import *

from unit_test import *


class KMCFrequencyPluginTest(unittest.TestCase):

    def setUp(self):
        # Test case setting.
        self.maxDiff = None
        self.setup = kmc_path + "/kmc_frequency_plugin.mkm"

    def test_run_with_frequency(self):
        " Make sure KMCSolver object can be constructed correctly. "
        model = KMCModel(setup_file=self.setup, verbosity=logging.WARNING)
        parser = model.parser
        parser.parse_data(filename=kmc_energy, relative=True)
        
        # Run the model with analysis.
        model.run(processes_file=kmc_processes,
                  configuration_file=kmc_config,
                  sitesmap_file=kmc_sites)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(KMCFrequencyPluginTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


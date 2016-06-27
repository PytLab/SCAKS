import logging
import re
import unittest

import numpy as np

from kynetix.model import KineticModel
from kynetix.solvers import *


class KMCFrequencyPluginTest(unittest.TestCase):

    def setUp(self):
        # Test case setting.
        self.maxDiff = None

    def test_run_with_frequency(self):
        " Make sure KMCSolver object can be constructed correctly. "
        model = KineticModel(setup_file="kmc_inputs/kmc_frequency_plugin.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        parser.parse_data(filename="kmc_inputs/rel_energy.py", relative=True)
        
        # Run the model with analysis.
        model.run_kmc(processes_file="kmc_inputs/kmc_processes.py",
                      configuration_file="kmc_inputs/kmc_configuration.py",
                      sitesmap_file="kmc_inputs/kmc_sites.py")

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(KMCFrequencyPluginTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


import logging
import re
import unittest

import numpy as np

from mikiac.models.kmc_model import KMCModel
from mikiac.solvers import *

from unit_test import *


class KMCFrequencyPluginTest(unittest.TestCase):

    def setUp(self):
        # Test case setting.
        self.maxDiff = None
        self.setup_dict = dict(
            rxn_expressions = [
                'CO_g + *_t -> CO_t',
                'CO_g + *_b -> CO_b',
                'O2_g + 2*_b -> 2O_b',
                'CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b',
                'CO_b + *_t <-> CO_t + *_b -> CO_b + *_t',
            ],

            species_definitions = {
                'CO_g': {'pressure': 0.01},
                'O2_g': {'pressure': 0.2},
                'CO2_g': {'pressure': 0.01},
                '*_b': {'site_name': 'bridge', 'type': 'site', 'total': 0.5},
                '*_t': {'site_name': 'top', 'type': 'site', 'total': 0.5},
            },

            temperature = 298.,
            parser = "KMCParser",
            solver = "KMCSolver",
            corrector = "ThermodynamicCorrector",
            cell_vectors = [[3.0, 0.0, 0.0],
                            [0.0, 3.0, 0.0],
                            [0.0, 0.0, 3.0]],
            basis_sites = [[0.0, 0.0, 0.0],
                           [0.5, 0.0, 0.0],
                           [0.0, 0.5, 0.0],
                           [0.5, 0.5, 0.0]],
            unitcell_area = 9.0e-20,
            active_ratio = 4./9,
            repetitions = (3, 3, 1),
            periodic = (True, True, False),
            possible_element_types = ["O", "V", "O_s", "C"],
            empty_type = "V",
            possible_site_types = ["P"],
            nstep = 50,
            random_seed = 13996,
            random_generator = 'MT',
            trajectory_dump_interval = 10,
            analysis = ["FrequencyAnalysis"],
            analysis_interval = [1],
        )

    def test_run_with_frequency(self):
        " Make sure KMCSolver object can be constructed correctly. "
        model = KMCModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        model.parser.parse_data(energy_file=kmc_energy,
                                processes_file=kmc_processes,
                                configuration_file=kmc_config,
                                sitesmap_file=kmc_sites)
        
        # Run the model with analysis.
        model.run()

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(KMCFrequencyPluginTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


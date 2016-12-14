#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import unittest

from kynetix.models.kmc_model import KMCModel
from kynetix.solvers import *

from unit_test import *


class KMCRedistributionTest(unittest.TestCase):

    def setUp(self):
        # Test case setting.
        self.maxDiff = None

    def test_run_with_split_redistribution(self):
        " Make sure the model can run with split redistribution operation. "
        setup_dict = dict(
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
                'b': {'site_name': 'bridge', 'type': 'site', 'total': 0.5},
                't': {'site_name': 'top', 'type': 'site', 'total': 0.5},
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
            repetitions = (4, 4, 1),
            periodic = (True, True, False),
            possible_element_types = ["O", "V", "O_s", "C"],
            empty_type = "V",
            possible_site_types = ["P"],
            nstep = 50,
            random_seed = 13996,
            random_generator = 'MT',
            trajectory_dump_interval = 10,
            do_redistribution = True,
            redistribution_interval = 5,
            fast_species = ["V"],
            nsplits = (2, 2, 1),
        )
        model = KMCModel(setup_dict=setup_dict, logger_level=logging.WARNING)
        model.parser.parse_data(relative=True,
                                energy_file=kmc_energy,
                                processes_file=kmc_processes,
                                configuration_file=kmc_config,
                                sitesmap_file=kmc_sites)

        # Run the model with redistribution.
        model.run()

    def test_run_with_process_redistribution(self):
        " Make sure the model can run with process redistribution operation. "
        setup_dict = dict(
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
                'b': {'site_name': 'bridge', 'type': 'site', 'total': 0.5},
                't': {'site_name': 'top', 'type': 'site', 'total': 0.5},
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
            repetitions = (4, 4, 1),
            periodic = (True, True, False),
            possible_element_types = ["O", "V", "O_s", "C"],
            empty_type = "V",
            possible_site_types = ["P"],
            nstep = 50,
            random_seed = 13996,
            random_generator = 'MT',
            trajectory_dump_interval = 10,
            do_redistribution = True,
            redistribution_interval = 5,
            distributor_type = "ProcessRandomDistributor",
        )
        model = KMCModel(setup_dict=setup_dict, logger_level=logging.WARNING)
        model.parser.parse_data(relative=True,
                                energy_file=kmc_energy,
                                processes_file=kmc_processes,
                                configuration_file=kmc_config,
                                sitesmap_file=kmc_sites)

        # Run the model with redistribution.
        model.run()

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(KMCRedistributionTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


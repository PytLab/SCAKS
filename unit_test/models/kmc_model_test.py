import logging
import sys
import unittest

sys.path.append("..")

from mikiac.models.kmc_model import KMCModel
from unit_test import *


class KMCModelTest(unittest.TestCase):

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
            time_seed = True,
            random_generator = 'MT',
            analysis = [],
            analysis_interval = [],
            trajectory_dump_interval = 10,
            distributor_type = "ProcessRandomDistributor",
        )

    def test_kmc_model_construction_query(self):
        " Test kmc model can be constructed with relative energy parser. "
        # Test construction.
        model = KMCModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)

        # Test parameters in setup file have been parsed.
        self.assertTrue(hasattr(model, "_KMCModel__cell_vectors"))
        self.assertTrue(hasattr(model, "_KMCModel__basis_sites"))
        self.assertTrue(hasattr(model, "_KMCModel__unitcell_area"))
        self.assertTrue(hasattr(model, "_KMCModel__active_ratio"))
        self.assertTrue(hasattr(model, "_KMCModel__repetitions"))
        self.assertTrue(hasattr(model, "_KMCModel__periodic"))
        self.assertTrue(hasattr(model, "_KMCModel__nstep"))
        self.assertTrue(hasattr(model, "_KMCModel__random_seed"))
        self.assertTrue(hasattr(model, "_KMCModel__time_seed"))
        self.assertTrue(hasattr(model, "_KMCModel__random_generator"))
        self.assertTrue(hasattr(model, "_KMCModel__analysis"))
        self.assertTrue(hasattr(model, "_KMCModel__analysis_interval"))
        self.assertTrue(hasattr(model, "_KMCModel__possible_element_types"))
        self.assertTrue(hasattr(model, "_KMCModel__possible_site_types"))
        self.assertTrue(hasattr(model, "_KMCModel__empty_type"))
        self.assertTrue(hasattr(model, "_KMCModel__distributor_type"))

    def test_kmc_model_run(self):
        " Make sure KMC model can run properly. "
        # Construction.
        model = KMCModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)

        # Parse data.
        model.parser.parse_data(energy_file=kmc_energy,
                                processes_file=kmc_processes,
                                configuration_file=kmc_config,
                                sitesmap_file=kmc_sites)

        # Run the model.
        model.run()

        # Run the model with default sites types.
        model.run()

        # Run with default configuration.
        model.run()

        # Run with default configuration and sitemap.
        model.run()

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(KMCModelTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


import logging
import sys
import unittest

sys.path.append("..")

from kynetix.model import KineticModel
from unit_test import *


class KMCModelTest(unittest.TestCase):

    def setUp(self):
        # Test case setting.
        self.maxDiff = None

    def test_kmc_model_construction_query(self):
        " Test kmc model can be constructed with relative energy parser. "
        # Test construction.
        kmc_setup = kmc_path + "/kmc_model.mkm"
        model = KineticModel(setup_file=kmc_setup, verbosity=logging.WARNING)

        # Test parameters in setup file have been parsed.
        self.assertTrue(hasattr(model, "_KineticModel__cell_vectors"))
        self.assertTrue(hasattr(model, "_KineticModel__basis_sites"))
        self.assertTrue(hasattr(model, "_KineticModel__unitcell_area"))
        self.assertTrue(hasattr(model, "_KineticModel__active_ratio"))
        self.assertTrue(hasattr(model, "_KineticModel__repetitions"))
        self.assertTrue(hasattr(model, "_KineticModel__periodic"))
        self.assertTrue(hasattr(model, "_KineticModel__nstep"))
        self.assertTrue(hasattr(model, "_KineticModel__seed"))
        self.assertTrue(hasattr(model, "_KineticModel__random_generator"))
        self.assertTrue(hasattr(model, "_KineticModel__analysis"))
        self.assertTrue(hasattr(model, "_KineticModel__analysis_interval"))
        self.assertTrue(hasattr(model, "_KineticModel__possible_element_types"))
        self.assertTrue(hasattr(model, "_KineticModel__possible_site_types"))
        self.assertTrue(hasattr(model, "_KineticModel__empty_type"))

    def test_kmc_model_run(self):
        " Make sure KMC model can run properly. "
        # Construction.
        kmc_setup = kmc_path + "/kmc_model.mkm"
        model = KineticModel(setup_file=kmc_setup, verbosity=logging.WARNING)
        parser = model.parser()
        solver = model.solver()

        # Parse data.
        parser.parse_data(filename=kmc_energy, relative=True)

        # Run the model.
        model.run_kmc(processes_file=kmc_processes,
                      configuration_file=kmc_config,
                      sitesmap_file=kmc_sites)

        # Run the model with default sites types.
        model.run_kmc(processes_file=kmc_processes,
                      configuration_file=kmc_config)

        # Run with default configuration.
        model.run_kmc(processes_file=kmc_processes,
                      sitesmap_file=kmc_sites)

        # Run with default configuration and sitemap.
        model.run_kmc(processes_file=kmc_processes)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(KMCModelTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


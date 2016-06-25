import unittest
import logging

from kynetix.model import KineticModel


class KMCModelTest(unittest.TestCase):

    def setUp(self):
        # Test case setting.
        self.maxDiff = None

    def test_kmc_model_construction_query(self):
        " Test kmc model can be constructed with relative energy parser. "
        # Test construction.
        model = KineticModel(setup_file="kmc_inputs/kmc_model.mkm",
                             verbosity=logging.WARNING)

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
        self.assertTrue(hasattr(model, "_KineticModel__color_dict"))
        self.assertTrue(hasattr(model, "_KineticModel__circle_attrs"))

    def test_kmc_model_run(self):
        " Make sure KMC model can run properly. "
        # Construction.
        model = KineticModel(setup_file="kmc_inputs/kmc_model.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        solver = model.solver()

        # Parse data.
        parser.parse_data(filename="kmc_inputs/rel_energy.py", relative=True)

        # Run the model.
        model.run_kmc(processes_file="kmc_inputs/kmc_processes.py",
                      configuration_file="kmc_inputs/kmc_configuration.py",
                      sitesmap_file="kmc_inputs/kmc_sites.py")

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(KMCModelTest)
    unittest.TextTestRunner(verbosity=2).run(suite)

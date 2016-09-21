import logging
import unittest

from kynetix.model import KineticModel
from kynetix.parsers import *

from unit_test import *


class MicroKineticModelTest(unittest.TestCase):

    def setUp(self):
        # Test case setting.
        self.maxDiff = None
        self.setup_file = mkm_path + "/mkm_model.mkm"

    def test_mkm_construction_query(self):
        " Test micro kinetic model can be constructed with parser. "
        # Test construction.
        model = KineticModel(setup_file=self.setup_file, verbosity=logging.WARNING)

        # Test member data query.
        self.assertEqual(model.h(), 4.135667662e-15)
        self.assertEqual(model.kB(), 8.6173324e-5)

        # Load data in setup file.
        glob, loc = {}, {}
        execfile(self.setup_file, glob, loc)
        self.assertEqual(model.tools(), loc["tools"])
        self.assertEqual(model.corrector(), loc["corrector"])
        self.assertEqual(model.plotter(), loc["plotter"])
        self.assertEqual(model.rxn_expressions(), loc["rxn_expressions"])
        self.assertEqual(model.temperature(), loc["temperature"])
        self.assertEqual(model.ref_species(), loc["ref_species"])
        self.assertEqual(model.surface_name(), loc["surface_name"])
        self.assertEqual(model.verbosity(), logging.WARNING)
        self.assertEqual(model.decimal_precision(), 100)

        self.assertTrue(isinstance(model.parser(), RelativeEnergyParser))

    def test_check_inputs(self):
        " Test the __check_inputs helper functions. "
        # Test construction.
        model = KineticModel(setup_file=self.setup_file, verbosity=logging.WARNING)

        # Check a valid input dict.
        inputs_dict = {"parser": "CsvParser",
                       "tools": ["parser", "solver"]}
        ret_inputs_dict = model._KineticModel__check_inputs(inputs_dict)

        self.assertIs(inputs_dict, ret_inputs_dict)

        # Check input dict with invalid parameter.
        inputs_dict = {"parser": "CsvParser",
                       "tools": ["parser", "solver"],
                       "whatname": (1, 2, 3)}
        ref_inputs_dict = {"parser": "CsvParser",
                           "tools": ["parser", "solver"]}
        ret_inputs_dict = model._KineticModel__check_inputs(inputs_dict)

        self.assertDictEqual(ret_inputs_dict, ref_inputs_dict)

    def test_run_mkm(self):
        " Test micro kinetic model can run correctly. "
        self.setup_file = mkm_path + "/mkm_model.mkm"
        model = KineticModel(setup_file=self.setup_file, verbosity=logging.WARNING)
        model.run_mkm(data_file=mkm_energy)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(MicroKineticModelTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


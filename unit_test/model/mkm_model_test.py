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

        # Load data in setup file.
        glob, loc = {}, {}
        execfile(self.setup_file, glob, loc)
        self.assertEqual(model.corrector.__class__.__name__, loc["corrector"])
        self.assertEqual(model.plotter.__class__.__name__, loc["plotter"])
        self.assertListEqual(model.rxn_expressions, loc["rxn_expressions"])
        self.assertEqual(model.temperature, loc["temperature"])
        self.assertListEqual(model.ref_species, loc["ref_species"])
        self.assertEqual(model.verbosity, logging.WARNING)
        self.assertEqual(model.decimal_precision, 100)

        self.assertTrue(isinstance(model.parser, RelativeEnergyParser))

    def test_run_mkm(self):
        " Test micro kinetic model can run correctly. "
        self.setup_file = mkm_path + "/mkm_model.mkm"
        model = KineticModel(setup_file=self.setup_file, verbosity=logging.WARNING)
        model.run_mkm(data_file=mkm_energy)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(MicroKineticModelTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


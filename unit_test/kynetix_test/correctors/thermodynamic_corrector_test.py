import logging
import re
import unittest

import numpy as np
from mpmath import mpf

from kynetix.model import KineticModel
from kynetix.correctors import *


class ThermodynamicCorrectorTest(unittest.TestCase):

    def setUp(self):
        # Test case setting.
        self.maxDiff = None

    def test_construction_and_query(self):
        " Test plotter construction and query. "
        # Construction.
        model = KineticModel(setup_file="input_files/thermodynamic_corrector.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        corrector = model.corrector()

        self.assertTrue(isinstance(corrector, ThermodynamicCorrector))

    def test_shomate_correction(self):
        " Test we can get correct energy correction by Shomate equation. "
        # Construction.
        model = KineticModel(setup_file="input_files/thermodynamic_corrector.mkm",
                             verbosity=logging.WARNING)
        corrector = model.corrector()

        # Check.
        gas = "CO_g"
        ret_delta = corrector.shomate_correction(gas)
        ref_delta = (-0.84190990132027166, 0.13611317171618073, -0.97802307303645242)
        self.assertTupleEqual(ret_delta, ref_delta)

        gas = "O2_g"
        ret_delta = corrector.shomate_correction(gas)
        ref_delta = (-0.87627116516400394, 0.13787890338175926, -1.0141500685457632)
        self.assertTupleEqual(ret_delta, ref_delta)

    def test_entropy_correction(self):
        " Make sure we can get correct entropy correction. "
        # Construction.
        model = KineticModel(setup_file="input_files/thermodynamic_corrector.mkm",
                             verbosity=logging.WARNING)
        corrector = model.corrector()

        # Check.
        gas = "CO_g"
        ret_delta = corrector.entropy_correction(gas)
        ref_delta = -1.1538116935108251
        self.assertEqual(ref_delta, ret_delta)

        gas = "O2_g"
        ret_delta = corrector.entropy_correction(gas)
        ref_delta = -1.2250150716175705
        self.assertEqual(ref_delta, ret_delta)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(ThermodynamicCorrectorTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


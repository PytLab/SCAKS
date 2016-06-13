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
        " Test we can get correct energy correction dict by Shomate equation. "
        # Construction.
        model = KineticModel(setup_file="input_files/thermodynamic_corrector.mkm",
                             verbosity=logging.WARNING)
        corrector = model.corrector()

        # Check.
        ref_dict = {'CO2_g': mpf('-0.89590348800350272373549387339153327047824859619140625'),
                    'CO_g': mpf('-0.84190990132027165859796014046878553926944732666015625'),
                    'O2_g': mpf('-0.87627116516400394008456942174234427511692047119140625')}
        ret_dict = corrector.shomate_correction()

        self.assertDictEqual(ref_dict, ret_dict)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(ThermodynamicCorrectorTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


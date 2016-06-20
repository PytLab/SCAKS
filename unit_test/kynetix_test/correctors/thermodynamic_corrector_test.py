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
        ref_delta = -0.84190990132027166
        self.assertEqual(ret_delta, ref_delta)

        gas = "O2_g"
        ret_delta = corrector.shomate_correction(gas)
        ref_delta = -0.87627116516400394
        self.assertEqual(ret_delta, ref_delta)

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

    def test_solvers_correction_energy(self):
        " Test solver's correction energy function. "
        model = KineticModel(setup_file="input_files/thermodynamic_corrector.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        parser.parse_data(filename="input_files/rel_energy.py")
        solver = model.solver()
        solver.get_data()

        ref_e = {'*_s': mpf('0.0'),
                 'CO-O_2s': mpf('0.9259999999999999342747969421907328069210052490234375'),
                 'CO2_g': mpf('0.0'),
                 'CO_g': mpf('0.0'),
                 'CO_s': mpf('-0.75800000000000000710542735760100185871124267578125'),
                 'O2_g': mpf('3.50800000000000000710542735760100185871124267578125'),
                 'O_s': mpf('0.4339999999999999413802242997917346656322479248046875')}
        ret_e = solver._G

        self.assertDictEqual(ref_e, ret_e)

        # Correction.
        solver.correct_energies()
        ref_e = {'*_s': mpf('0.0'),
                 'CO-O_2s': mpf('0.9259999999999999342747969421907328069210052490234375'),
                 'CO2_g': mpf('-0.89590348800350272373549387339153327047824859619140625'),
                 'CO_g': mpf('-0.84190990132027165859796014046878553926944732666015625'),
                 'CO_s': mpf('-0.75800000000000000710542735760100185871124267578125'),
                 'O2_g': mpf('2.63172883483599606702085793585865758359432220458984375'),
                 'O_s': mpf('0.4339999999999999413802242997917346656322479248046875')}
        ret_e = solver._G

        self.assertDictEqual(ref_e, ret_e)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(ThermodynamicCorrectorTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


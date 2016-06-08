import logging
import re
import unittest

import numpy as np

from kynetix.parsers.rxn_parser import *


class RxnEquationTest(unittest.TestCase):

    def setUp(self):
        # Test case setting.
        self.maxDiff = None

    def test_check_conservation(self):
        " Test conservation checking function. "
        # Construction.
        equation = RxnEquation("CO_s + O_s <-> CO-O_2s -> CO2_g + 2*_s")
        self.assertTrue(equation.check_conservation())

        # Type Exception.
        equation = RxnEquation("CO_s + K_s <-> CO-O_2s -> CO2_g + 2*_s")
        self.assertRaisesRegexp(ValueError, r"^Mass", equation.check_conservation)

        equation = RxnEquation("CO_s + O2_s <-> CO-O_2s -> CO2_g + 2*_s")
        self.assertRaisesRegexp(ValueError, r"^Mass", equation.check_conservation)

        # Site Exception.
        equation = RxnEquation("CO_s + O_2s <-> CO-O_2s -> CO2_g + 2*_s")
        self.assertRaisesRegexp(ValueError, r"^Site", equation.check_conservation)

        equation = RxnEquation("CO_s + O_t <-> CO-O_2s -> CO2_g + 2*_s")
        self.assertRaisesRegexp(ValueError, r"^Site", equation.check_conservation)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(RxnEquationTest)
    unittest.TextTestRunner(verbosity=2).run(suite)

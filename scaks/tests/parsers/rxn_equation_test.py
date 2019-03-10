import logging
import re
import unittest

import numpy as np

from ...parsers.rxn_parser import *
from .. import cleanup

class RxnEquationTest(unittest.TestCase):

    def setUp(self):
        # Test case setting.
        self.maxDiff = None

    def test_to_formula_list(self):
        " Make sure we can get correct formula list for a rxn equation. "
        # Construction.
        equation = RxnEquation("CO_s + O_s <-> CO-O_2s -> CO2_g + 2*_s")

        ref_formula_list = [['CO_s', 'O_s'], ['CO-O_2s'], ['CO2_g', '2*_s']]
        ret_formula_list = equation.to_formula_list()

        # Check.
        for state, state_str in zip(ret_formula_list, ref_formula_list):
            for formula, formula_str in zip(state, state_str):
                self.assertEqual(formula.formula(), formula_str)

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

    def test_revert(self):
        " Test the reaction can revert correctly. "
        equation = RxnEquation("CO_s + O_s <-> CO-O_2s -> CO2_g + 2*_s")
        ref_reverse = "CO2_g + 2*_s <-> CO-O_2s -> CO_s + O_s"
        ret_reverse = equation.revert().rxn_equation()

        equation = RxnEquation("CO_s + O_s -> CO2_g + 2*_s")
        ref_reverse = "CO2_g + 2*_s -> CO_s + O_s"
        ret_reverse = equation.revert().rxn_equation()

    def test_adsorption_gases(self):
        " Make sure we can get the correct gas species in IS. "
        equation = RxnEquation("CO2_g + 2*_s -> CO_s + O_s")
        gases = equation.adsorption_gases()
        ref_gas_names = ["CO2_g"]
        ret_gas_names = [gas.formula() for gas in gases]

        self.assertListEqual(ref_gas_names, ret_gas_names)

    def test_desorption_gases(self):
        " Make sure we can get the correct gas species in FS. "
        equation = RxnEquation("CO_s + O_s <-> CO-O_2s -> CO2_g + 2*_s")
        gases = equation.desorption_gases()
        ref_gas_names = ["CO2_g"]
        ret_gas_names = [gas.formula() for gas in gases]

        self.assertListEqual(ref_gas_names, ret_gas_names)

    def test_adsorbate_search(self):
        " Make sure adsorbate can be searched correctly "
        adsorbate = 'CO_s'
        equation = RxnEquation("CO_s + O_s <-> CO-O_2s -> CO2_g + 2*_s")
        res = equation.search_intermediate(adsorbate)
        self.assertEqual(res['state'], 'IS')
        self.assertEqual(res['intermediate'].formula, adsorbate)
        self.assertEqual(res['intermediate'].stoichiometry, 1)

        res = equation.search_intermediate('H_s')
        self.assertIsNone(res)

    def tearDown(self):
        cleanup()

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(RxnEquationTest)
    unittest.TextTestRunner(verbosity=2).run(suite)

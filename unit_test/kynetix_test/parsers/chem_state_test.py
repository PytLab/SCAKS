import logging
import re
import unittest

import numpy as np

from kynetix.parsers.rxn_parser import *


class ChemStateTest(unittest.TestCase):

    def setUp(self):
        # Test case setting.
        self.maxDiff = None

    def test_get_elements_dict(self):
        " Test we can get correct elements dictionary. "
        # Construction.
        state_str = "CO2_g + 2*_s"
        state = ChemState(state_str)

        # Test get_elements_dict.
        ref_elements_dict = {'C': 1, 'O': 2}
        ret_elements_dict = state.get_elements_dict()
        self.assertDictEqual(ref_elements_dict, ret_elements_dict)

        # More stoichmetries.
        state_str = "2CO2_g + 2*_s"
        state = ChemState(state_str)

        # Test get_elements_dict.
        ref_elements_dict = {'C': 2, 'O': 4}
        ret_elements_dict = state.get_elements_dict()
        self.assertDictEqual(ref_elements_dict, ret_elements_dict)

    def test_get_sites_dict(self):
        " Test we can get sites dict correctly. "
        # Construction.
        state_str = "CO2_g + 2*_s"
        state = ChemState(state_str)

        # Test get_elements_dict.
        ref_sites_dict = {'s': 2}
        ret_sites_dict = state.get_sites_dict()
        self.assertDictEqual(ref_sites_dict, ret_sites_dict)

    def test_conserve(self):
        " Test conservation checking function. "
        # Construction.
        state_str = "CO2_s + O_s"
        state_1 = ChemState(state_str)

        state_str = "CO-O_s + O_s"
        state_2 = ChemState(state_str)

        self.assertTrue(state_2.conserve(state_1))
        self.assertTrue(state_1.conserve(state_2))

        # Once again.
        state_str = "CO2_g + 2*_s"
        state_1 = ChemState(state_str)

        state_str = "CO-O_2s"
        state_2 = ChemState(state_str)

        self.assertTrue(state_2.conserve(state_1))
        self.assertTrue(state_1.conserve(state_2))

        # Exception check.
        state_str = "CO_g + 2*_s"
        state_1 = ChemState(state_str)

        state_str = "CO-O_2s"
        state_2 = ChemState(state_str)

        self.assertRaisesRegexp(ValueError, r"^Mass", state_1.conserve, state_2)
        self.assertRaisesRegexp(ValueError, r"^Mass", state_2.conserve, state_1)

        # Exception check.
        state_str = "CO2_g + *_s"
        state_1 = ChemState(state_str)

        state_str = "CO-O_2s"
        state_2 = ChemState(state_str)

        self.assertRaisesRegexp(ValueError, r"^Site", state_1.conserve, state_2)
        self.assertRaisesRegexp(ValueError, r"^Site", state_2.conserve, state_1)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(ChemStateTest)
    unittest.TextTestRunner(verbosity=2).run(suite)

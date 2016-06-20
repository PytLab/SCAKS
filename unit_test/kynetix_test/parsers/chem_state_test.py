import logging
import re
import unittest

import numpy as np

from kynetix.parsers.rxn_parser import *


class ChemStateTest(unittest.TestCase):

    def setUp(self):
        # Test case setting.
        self.maxDiff = None

    def test_split(self):
        " Test we can get formula string list correctly. "
        # Construction.
        state_str = "CO2_g + 2*_s"
        state = ChemState(state_str)

        ref_formula_list = ['CO2_g', '2*_s']
        ret_formula_list = state.split()

        self.assertListEqual(ref_formula_list, ret_formula_list)

    def test_tolist(self):
        " Test we can get correct ChemFormula list. "
        # Construction.
        state_str = "CO2_g + 2*_s"
        state = ChemState(state_str)

        formula_list = state.tolist()

        # Check length.
        self.assertEqual(len(formula_list), 2)

        # Check formula.
        formula_str_list = ['CO2_g', '2*_s']
        for formula, formula_str in zip(formula_list, formula_str_list):
            self.assertTrue(isinstance(formula, ChemFormula))
            self.assertEqual(formula.formula(), formula_str)

    def test_get_species_site_list(self):
        " Make sure we can get species_site list correctly. "
        # Construction.
        state_str = "CO2_g + 2*_s"
        state = ChemState(state_str)

        ref_species_site_list = ['CO2_g', '*_s']
        ret_species_site_list = state.get_species_site_list()

        self.assertListEqual(ref_species_site_list, ret_species_site_list)

    def test_get_species_site_dict(self):
        " Make sure we can get species_site dict correctly. "
        # Construction.
        state_str = "CO2_g + 2*_s"
        state = ChemState(state_str)

        ref_species_site_dict = {'*_s': 2, 'CO2_g': 1}
        ret_species_site_dict = state.get_species_site_dict()

        self.assertDictEqual(ref_species_site_dict, ret_species_site_dict)

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

    def test_texen(self):
        " Make sure we can get correct tex string for a state. "
        state_str = "HCOOH_g + 2*_s"
        state = ChemState(state_str)

        ref_tex = r'HCOOH(g) + 2*(s)'
        ret_tex = state.texen()

        self.assertEqual(ref_tex, ret_tex)

        state_str = "HCOO-H_s + *_s"
        state = ChemState(state_str)

        ref_tex = r'HCOO-H(s) + *(s)'
        ret_tex = state.texen()

        self.assertEqual(ref_tex, ret_tex)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(ChemStateTest)
    unittest.TextTestRunner(verbosity=2).run(suite)

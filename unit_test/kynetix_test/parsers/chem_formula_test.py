import logging
import re
import unittest

import numpy as np

from kynetix.parsers.rxn_parser import *


class ChemFormulaTest(unittest.TestCase):

    def setUp(self):
        # Test case setting.
        self.maxDiff = None

    def test_type(self):
        " Make sure formula can return correct species type. "
        # Gas.
        formula_str = "CO2_g"
        formula = ChemFormula(formula_str)
        self.assertEqual("gas", formula.type())

        # Liquid.
        formula_str = "CO2_l"
        formula = ChemFormula(formula_str)
        self.assertEqual("liquid", formula.type())

        # Site.
        formula_str = "3*_s"
        formula = ChemFormula(formula_str)
        self.assertEqual("site", formula.type())

        # Adsorbate.
        formula_str = "CO2_s"
        formula = ChemFormula(formula_str)
        self.assertEqual("adsorbate", formula.type())

    def test_get_elements_dict(self):
        " Test we can get correct elements dictionary. "
        # Construction.
        formula_str = "CO2_g"
        formula = ChemFormula(formula_str)

        # Test get_elements_dict.
        ref_elements_dict = {'C': 1, 'O': 2}
        ret_elements_dict = formula.get_elements_dict()
        self.assertDictEqual(ref_elements_dict, ret_elements_dict)

        # More than one molecules.
        formula_str = "3CO2_g"
        formula = ChemFormula(formula_str)

        # Test get_elements_dict.
        ref_elements_dict = {'C': 3, 'O': 6}
        ret_elements_dict = formula.get_elements_dict()
        self.assertDictEqual(ref_elements_dict, ret_elements_dict)

        # Construction.
        formula_str = "3*_s"
        formula = ChemFormula(formula_str)

        # Test get_elements_dict.
        ref_elements_dict = {}
        ret_elements_dict = formula.get_elements_dict()
        self.assertDictEqual(ref_elements_dict, ret_elements_dict)

    def test_get_sites_dict(self):
        " Make sure we can get correct site dictionary. "
        # Construction.
        formula_str = "CO2_g"
        formula = ChemFormula(formula_str)

        ref_sites_dict = {'g': 1}
        ret_sites_dict = formula.get_sites_dict()
        self.assertDictEqual(ref_sites_dict, ret_sites_dict)

        # Construction.
        formula_str = "3*_s"
        formula = ChemFormula(formula_str)

        ref_sites_dict = {'s': 3}
        ret_sites_dict = formula.get_sites_dict()
        self.assertDictEqual(ref_sites_dict, ret_sites_dict)

    def test_conserve(self):
        " Test conservation checking function. "
        # Construction.
        formula_str = "CO2_s"
        formula_1 = ChemFormula(formula_str)

        formula_str = "CO-O_s"
        formula_2 = ChemFormula(formula_str)

        self.assertTrue(formula_1.conserve(formula_2))
        self.assertTrue(formula_2.conserve(formula_1))

        # Test elements nonconservation.
        formula_str = "CO_s"
        formula_1 = ChemFormula(formula_str)

        formula_str = "CO-O_s"
        formula_2 = ChemFormula(formula_str)

        self.assertRaisesRegexp(ValueError, r"^Mass", formula_2.conserve, formula_1)

        # Test sites nonconservation.
        formula_str = "CO2_2s"
        formula_1 = ChemFormula(formula_str)

        formula_str = "CO-O_s"
        formula_2 = ChemFormula(formula_str)

        self.assertRaisesRegexp(ValueError, r"^Site", formula_2.conserve, formula_1)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(ChemFormulaTest)
    unittest.TextTestRunner(verbosity=2).run(suite)

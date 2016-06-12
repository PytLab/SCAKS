import logging
import re
import unittest

import numpy as np

from kynetix.model import KineticModel
from kynetix.parsers import *


class ParserBaseTest(unittest.TestCase):

    def setUp(self):
        # Test case setting.
        self.maxDiff = None

    def test_parser_construction(self):
        " Test parser can be constructed in kinetic model. "
        # Construction.
        model = KineticModel(setup_file="input_files/parser_base.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()

        # Check the parser class and base class type.
        self.assertTrue(isinstance(parser, RelativeEnergyParser))
        self.assertEqual(parser.__class__.__base__.__name__, "ParserBase")

#    def test_parser_base_query(self):
#        " Test parser base query functions. "
#        # Construction.
#        model = KineticModel(setup_file="input_files/parser_base.mkm",
#                             verbosity=logging.WARNING)
#        parser = model.parser()
#
#        # Test species definitions.
#        ref_species_definitions = {'CO-O_2s': {'elements': {'C': 1, 'O': 2},
#                                    'site': 's',
#                                    'site_number': 2,
#                                    'type': 'transition_state'},
#                                   'CO2_g': {'elements': {'C': 1, 'O': 2},
#                                    'pressure': 0.0,
#                                    'site': 'g',
#                                    'site_number': 1,
#                                    'type': 'gas'},
#                                   'CO_g': {'elements': {'C': 1, 'O': 1},
#                                    'pressure': 1.0,
#                                    'site': 'g',
#                                    'site_number': 1,
#                                    'type': 'gas'},
#                                   'CO_s': {'elements': {'C': 1, 'O': 1},
#                                    'site': 's',
#                                    'site_number': 1,
#                                    'type': 'adsorbate'},
#                                   'O2_g': {'elements': {'O': 2},
#                                    'pressure': 0.3333333333333333,
#                                    'site': 'g',
#                                    'site_number': 1,
#                                    'type': 'gas'},
#                                   'O_s': {'elements': {'O': 1},
#                                    'site': 's',
#                                    'site_number': 1,
#                                    'type': 'adsorbate'},
#                                   's': {'site_name': '111', 'total': 1.0, 'type': 'site'}}
#        ret_species_definitions = parser.species_definitions()
#
#        self.assertDictEqual(ret_species_definitions, ref_species_definitions)

    def test_stoichiometry_matrices(self):
        " Make sure we can get the reactant product matrix and intermediate matrix correctly."
        # Construction.
        model = KineticModel(setup_file="input_files/parser_base.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()

        ref_reapro_matrix = np.matrix([[1.0, -1.0, 0.0],
                                       [2.0, 0.0, -2.0],
                                       [-2.0, 1.0, 1.0]])
        ref_site_matrix = np.matrix([[0.0, 1.0, 0.0],
                                    [0.0, 0.0, 1.0],
                                    [-1.0, 0.0, 0.0]])
        ret_reapro_matrix, ret_site_matrix = parser.get_stoichiometry_matrices()

        self.assertTrue(np.allclose(ref_reapro_matrix, ret_reapro_matrix))
        self.assertTrue(np.allclose(ref_site_matrix, ret_site_matrix))

    def test_elemtary_rxns_parse(self):
        " Test all elementary reaction equations can be parsed correctly. "

        # Construction.
        model = KineticModel(setup_file="input_files/parser_base.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()

        elementary_rxns = ['CO_g + *_s -> CO_s',
                           'O2_g + 2*_s -> 2O_s',
                           'CO_s + O_s <-> CO-O_2s -> CO2_g + 2*_s']

        ref_adsorbate_names = ('CO_s', 'O_s')
        ref_gas_names = ('CO2_g', 'CO_g', 'O2_g')
        ref_liquid_names = ()
        ref_site_names = ('s', )
        ref_transition_state_names = ('CO-O_2s', )

        (ret_adsorbate_names,
         ret_gas_names,
         ret_liquid_names,
         ret_site_names,
         ret_transition_state_names,
         ret_elementary_rxns_list) = parser.parse_elementary_rxns(elementary_rxns)

        # Check.
        self.assertTupleEqual(ref_adsorbate_names, ret_adsorbate_names)
        self.assertTupleEqual(ref_gas_names, ret_gas_names)
        self.assertTupleEqual(ref_liquid_names, ret_liquid_names)
        self.assertTupleEqual(ref_site_names, ret_site_names)
        self.assertTupleEqual(ref_transition_state_names, ret_transition_state_names)

    def test_total_rxn_equation(self):
        " Test we can get the total reaction equation from elementary reactions. "

        # Construction.
        model = KineticModel(setup_file="input_files/parser_base.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()

        ref_total_rxn_equation = "2CO_g + O2_g -> 2CO2_g"
        ret_total_rxn_equation = parser.get_total_rxn_equation()

        # Check.
        self.assertEqual(ref_total_rxn_equation, ret_total_rxn_equation)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(ParserBaseTest)
    unittest.TextTestRunner(verbosity=2).run(suite)

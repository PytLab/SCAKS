import logging
import re
import unittest

import numpy as np

from kynetix.model import KineticModel
from kynetix.parsers import *


class TestParserBase(unittest.TestCase):

    def setUp(self):
        # Test case setting.
        self.maxDiff = None

    def test_parser_construction(self):
        " Test parser can be constructed in kinetic model. "
        # Construction.
        model = KineticModel(setup_file="input_files/setup.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()

        # Check the parser class and base class type.
        self.assertTrue(isinstance(parser, RelativeEnergyParser))
        self.assertEqual(parser.__class__.__base__.__name__, "ParserBase")

    def test_parser_base_query(self):
        " Test parser base query functions. "
        # Construction.
        model = KineticModel(setup_file="input_files/setup.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()

        # Test regex_dict().
        ref_regex_dict = {}

        states_regex = re.compile(r'([^\<\>]*)(?:\<?\-\>)' +
                                  r'(?:([^\<\>]*)(?:\<?\-\>))?([^\<\>]*)')
        ref_regex_dict['IS_TS_FS'] = [states_regex, ['IS', 'TS', 'FS']]

        species_regex = re.compile(r'(\d*)([^\_\+\*\<\>]+)_(\d*)(\w+)')
        ref_regex_dict['species'] = \
            [species_regex, ['stoichiometry', 'name', 'site_number', 'site']]

        site_regex = re.compile(r'(\d*)(?:\*\_)(\w+)')
        ref_regex_dict['empty_site'] = [site_regex, ['stoichiometry', 'site']]

        self.assertDictEqual(parser.regex_dict(), ref_regex_dict)

        # Test species definitions.
        ref_species_definitions = {'CO-O_2s': {'elements': {'C': 1, 'O': 2},
                                    'site': 's',
                                    'site_number': 2,
                                    'type': 'transition_state'},
                                   'CO2_g': {'elements': {'C': 1, 'O': 2},
                                    'pressure': 0.0,
                                    'site': 'g',
                                    'site_number': 1,
                                    'type': 'gas'},
                                   'CO_g': {'elements': {'C': 1, 'O': 1},
                                    'pressure': 1.0,
                                    'site': 'g',
                                    'site_number': 1,
                                    'type': 'gas'},
                                   'CO_s': {'elements': {'C': 1, 'O': 1},
                                    'site': 's',
                                    'site_number': 1,
                                    'type': 'adsorbate'},
                                   'O2_g': {'elements': {'O': 2},
                                    'pressure': 0.3333333333333333,
                                    'site': 'g',
                                    'site_number': 1,
                                    'type': 'gas'},
                                   'O_s': {'elements': {'O': 1},
                                    'site': 's',
                                    'site_number': 1,
                                    'type': 'adsorbate'},
                                   's': {'site_name': '111', 'total': 1.0, 'type': 'site'}}
        ret_species_definitions = parser.species_definitions()

        self.assertDictEqual(ret_species_definitions, ref_species_definitions)

    def test_state_expression_parse(self):
        " Test state expression can be parsed correctly. "

        # Construction.
        model = KineticModel(setup_file="input_files/setup.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()

        # Test a normal state.
        state_expression = "CO_g + *_s"
        ref_state_dict = {'CO_g': {'elements': {'C': 1, 'O': 1},
                                   'number': 1,
                                   'site': 'g',
                                   'site_number': 1}}

        ref_site_dict = {'s': {'number': 1, 'type': 's'}}
        ref_species_list = ['CO_g', '*_s']
        ret_state_dict, ret_site_dict, ret_species_list = \
            parser._ParserBase__parse_state_expression(state_expression)

        self.assertDictEqual(ref_state_dict, ret_state_dict)
        self.assertDictEqual(ref_site_dict, ret_site_dict)
        self.assertListEqual(ref_species_list, ref_species_list)

        # Test a 2-sites state.
        state_expression = "CO-O_2s"
        ref_state_dict = {'CO-O_2s': {'elements': {'C': 1, 'O': 2},
                                      'number': 1,
                                      'site': 's',
                                      'site_number': 2}}
        ref_site_dict = {}
        ref_species_list = ['CO-O_2s']

        ret_state_dict, ret_site_dict, ret_species_list = \
            parser._ParserBase__parse_state_expression(state_expression)

        self.assertDictEqual(ref_state_dict, ret_state_dict)
        self.assertDictEqual(ref_site_dict, ret_site_dict)
        self.assertListEqual(ref_species_list, ref_species_list)

    def test_species_expression_parse(self):
        " Test species expression can be parsed correctly. "
        # Construction.
        model = KineticModel(setup_file="input_files/setup.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()

        # Gas.
        species = "CO_g"
        ref_sp_dict = {'CO_g': {'elements': {'C': 1, 'O': 1},
                                'number': 1,
                                'site': 'g',
                                'site_number': 1}}
        ret_sp_dict = parser._ParserBase__parse_species_expression(species)
        self.assertDictEqual(ref_sp_dict, ret_sp_dict)

        # Intermediate.
        species = "C2H4_s"
        ref_sp_dict = {'C2H4_s': {'elements': {'C': 2, 'H': 4},
                                  'number': 1,
                                  'site': 's',
                                  'site_number': 1}}
        ret_sp_dict = parser._ParserBase__parse_species_expression(species)
        self.assertDictEqual(ref_sp_dict, ret_sp_dict)

        # Transition state.
        species = "O-O_2s"
        ref_sp_dict = {'O-O_2s': {'elements': {'O': 2},
                                  'number': 1,
                                  'site': 's',
                                  'site_number': 2}}
        ret_sp_dict = parser._ParserBase__parse_species_expression(species)
        self.assertDictEqual(ref_sp_dict, ret_sp_dict)

    def test_site_expression_parse(self):
        " Test site expression can be parsed correctly. "
        # Construction.
        model = KineticModel(setup_file="input_files/setup.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()

        site = "2*_s"
        ref_site_dict = {'s': {'number': 2, 'type': 's'}}
        ret_site_dict = parser._ParserBase__parse_site_expression(site)
        self.assertDictEqual(ref_site_dict, ret_site_dict)

        site = "*_t"
        ref_site_dict = {'t': {'number': 1, 'type': 't'}}
        ret_site_dict = parser._ParserBase__parse_site_expression(site)
        self.assertDictEqual(ref_site_dict, ret_site_dict)

    def test_stoichiometry_matrices(self):
        " Make sure we can get the reactant product matrix and intermediate matrix correctly."
        # Construction.
        model = KineticModel(setup_file="input_files/setup.mkm",
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

    def test_update_species_definitions(self):
        " Test species definitions of parser can be updated correctly. "

        # Construction.
        model = KineticModel(setup_file="input_files/setup.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()

        # Test species original definitions.
        ref_species_definitions = {'CO2_g': {'pressure': 0.0},
                                  'CO_g': {'pressure': 1.0},
                                  'O2_g': {'pressure': 0.3333333333333333},
                                  's': {'site_name': '111', 'total': 1.0, 'type': 'site'}}

        # Here we change the default value of parser's definitions.
        parser._ParserBase__species_definitions = ref_species_definitions

        # Em... this assertion maybe redundant.
        ret_species_definitions = parser.species_definitions()
        self.assertDictEqual(ret_species_definitions, ref_species_definitions)

        # Update with one species dict.
        sp_dict = {'CO_g': {'elements': {'C': 1, 'O': 1},
                            'number': 1,
                            'site': 'g',
                            'site_number': 1}}
        parser._ParserBase__update_species_definitions(sp_dict)

        ref_species_definitions = {'CO2_g': {'pressure': 0.0},
                                   'CO_g': {'elements': {'C': 1, 'O': 1},
                                            'pressure': 1.0,
                                            'site': 'g',
                                            'site_number': 1,
                                            'type': 'gas'},
                                   'O2_g': {'pressure': 0.3333333333333333},
                                   's': {'site_name': '111', 'total': 1.0, 'type': 'site'}}
        self.assertDictEqual(parser.species_definitions(), ref_species_definitions)

    def test_single_rxn_expression_parse(self):
        " Test a single reaction expression can be parsed correctly. "

        # Construction.
        model = KineticModel(setup_file="input_files/setup.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()

        rxn_expression = "CO_s + O_s <-> CO-O_2s -> CO2_g + 2*_s"
        ref_states_dict = {'FS': {'empty_sites_dict': {'s': {'number': 2, 'type': 's'}},
                            'species_dict': {'CO2_g': {'elements': {'C': 1, 'O': 2},
                              'number': 1,
                              'site': 'g',
                              'site_number': 1}},
                            'state_expression': 'CO2_g + 2*_s'},
                           'IS': {'empty_sites_dict': {},
                            'species_dict': {'CO_s': {'elements': {'C': 1, 'O': 1},
                              'number': 1,
                              'site': 's',
                              'site_number': 1},
                             'O_s': {'elements': {'O': 1}, 'number': 1, 'site': 's', 'site_number': 1}},
                            'state_expression': 'CO_s + O_s'},
                           'TS': {'empty_sites_dict': {},
                            'species_dict': {'CO-O_2s': {'elements': {'C': 1, 'O': 2},
                              'number': 1,
                              'site': 's',
                              'site_number': 2}},
                            'state_expression': 'CO-O_2s'}}
        ref_rxn_list = [['CO_s', 'O_s'], ['CO-O_2s'], ['CO2_g', '2*_s']]
        ret_states_dict, ret_rxn_list = parser.parse_single_elementary_rxn(rxn_expression)

        # Check.
        self.assertDictEqual(ref_states_dict, ret_states_dict)
        self.assertListEqual(ref_rxn_list, ret_rxn_list)

    def test_elemtary_rxns_parse(self):
        " Test all elementary reaction equations can be parsed correctly. "

        # Construction.
        model = KineticModel(setup_file="input_files/setup.mkm",
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
        ref_elementary_rxns_list = [[['CO_g', '*_s'], ['CO_s']],
                                    [['O2_g', '2*_s'], ['2O_s']],
                                    [['CO_s', 'O_s'], ['CO-O_2s'], ['CO2_g', '2*_s']]]
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
        self.assertListEqual(ref_elementary_rxns_list, ret_elementary_rxns_list)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestParserBase)
    unittest.TextTestRunner(verbosity=2).run(suite)


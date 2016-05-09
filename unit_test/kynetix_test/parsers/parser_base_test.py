import unittest
import logging
import re

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

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestParserBase)
    unittest.TextTestRunner(verbosity=2).run(suite)


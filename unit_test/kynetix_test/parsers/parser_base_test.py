import unittest
import logging
import re

from kynetix.model import KineticModel
from kynetix.parsers import *


class TestParserBase(unittest.TestCase):

    def setUp(self):
        # Test case setting.
        self.maxDiff = None

        # Construction.
        model = KineticModel(setup_file="input_files/setup.mkm",
                             verbosity=logging.WARNING)
        self.parser = model.parser()

    def test_parser_construction(self):
        " Test parser can be constructed in kinetic model. "
        # Check the parser class and base class type.
        self.assertTrue(isinstance(self.parser, RelativeEnergyParser))
        self.assertEqual(self.parser.__class__.__base__.__name__, "ParserBase")

    def test_parser_base_query(self):
        " Test parser base query functions. "

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

        self.assertDictEqual(self.parser.regex_dict(), ref_regex_dict)

    def test_state_expression_parse(self):
        " Test state expression can be parsed correctly. "

        # Test a normal state.
        state_expression = "CO_g + *_s"
        ref_state_dict = {'CO_g': {'elements': {'C': 1, 'O': 1},
                                   'number': 1,
                                   'site': 'g',
                                   'site_number': 1}
                         }

        ref_site_dict = {'s': {'number': 1, 'type': 's'}}
        ref_species_list = ['CO_g', '*_s']
        ret_state_dict, ret_site_dict, ret_species_list = \
            self.parser._ParserBase__parse_state_expression(state_expression)

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
            self.parser._ParserBase__parse_state_expression(state_expression)

        self.assertDictEqual(ref_state_dict, ret_state_dict)
        self.assertDictEqual(ref_site_dict, ret_site_dict)
        self.assertListEqual(ref_species_list, ref_species_list)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestParserBase)
    unittest.TextTestRunner(verbosity=2).run(suite)


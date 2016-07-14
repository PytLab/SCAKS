import logging
import re
import unittest

from kynetix.errors.error import *
from kynetix.model import KineticModel
from kynetix.utilities.check_utilities import *


class CheckUtilitiesTest(unittest.TestCase):

    def setUp(self):
        # Test case setting.
        self.maxDiff = None

    def test_check_list_tuple(self):
        " Test check_list_tuple function can work properly. "
        # Check valid string sequence.
        tools = ['parser', 'solver']
        check_list_tuple(tools, str, "tools")

        # Check invalid sequence with different type.
        tools = ['parser', 12]
        self.assertRaises(SetupError, check_list_tuple, tools,
                          entry_type=str, param_name="tools")

        # Not a sequence.
        tools = 'parser'
        self.assertRaises(SetupError, check_list_tuple, tools,
                          entry_type=str, param_name="tools")

    def test_check_string(self):
        "Test check_string function can work correctly."
        parsers_range = ("RelativeEneryParser", "CsvParser", "KMCParser")

        # Check a valid string.
        parser = "KMCParser"
        check_string(parser, string_range=parsers_range, param_name="parser")

        # Invalid string (not in range).
        parser = "AsdParser"
        self.assertRaisesRegexp(SetupError,
                                r".+ must be one of \(.+, .+, .+\)",
                                check_string,
                                parser,
                                string_range=parsers_range,
                                param_name="parser")

        # Invalid string (not a string).
        parser = 123
        self.assertRaisesRegexp(SetupError,
                                r".+ is not a string",
                                check_string,
                                parser,
                                string_range=parsers_range,
                                param_name="parser")

    def test_species_definitions(self):
        " Test species definitions can be check correctly. "

        # A valid species definitions.
        species_definitions = {}
        species_definitions['CO_g'] = {'pressure': 1.0}    # define the gas pressures
        species_definitions['O2_g'] = {'pressure': 1./3.}  # 0.094
        species_definitions['CO2_g'] = {'pressure': 0.00}
        species_definitions['s'] = {'site_name': '111', 'type': 'site', 'total': 1.0}

        # Check.
        check_species_definitions(species_definitions)

        # Invalid (not a dict).
        species_definitions = (1, 2, 3)
        self.assertRaisesRegexp(SetupError,
                                "species definitions is not a dict.",
                                check_species_definitions,
                                species_definitions)

        # Invalid (pressure not set).
        species_definitions = {}
        species_definitions['CO_g'] = {}    # define the gas pressures
        species_definitions['O2_g'] = {'pressure': 1./3.}  # 0.094
        species_definitions['CO2_g'] = {'pressure': 0.00}
        species_definitions['s'] = {'site_name': '111', 'type': 'site', 'total': 1.0}

        self.assertRaisesRegexp(SetupError,
                                "No pressure info for gas species CO_g.",
                                check_species_definitions,
                                species_definitions)

    def test_check_process_dict(self):
        " Test process dict can be checked correctly. "
        # {{{
        # Check keys.
        process = {"reaction": "CO_g + *_t -> CO_t",
                   "description": "CO adsorption at top site.",
                   "elements_before": ["V"],
                   "elements_after": ["C"],
                   "basis_sites": [0]}

        self.assertRaisesRegexp(SetupError,
                                r"^key '\w+' is not in process_dict",
                                check_process_dict,
                                process)

        # Check reaction string.
        process = {"reaction": ["CO_g + *_t -> CO_t"],
                   "description": "CO adsorption at top site.",
                   "coordinates_group": [[[0.0, 0.0, 0.0]]],
                   "elements_before": ["V"],
                   "elements_after": ["C"],
                   "basis_sites": [0]}

        self.assertRaisesRegexp(SetupError,
                                r"^reaction must be a string.$",
                                check_process_dict,
                                process)

        # Check coordinates.
        process = {"reaction": "CO_g + *_t -> CO_t",
                   "description": "CO adsorption at top site.",
                   "coordinates_group": [[[0.0, 0.0]]],
                   "elements_before": ["V"],
                   "elements_after": ["C"],
                   "basis_sites": [0]}

        self.assertRaisesRegexp(SetupError,
                                r"must have 3 entries.$",
                                check_process_dict,
                                process)

        # Check elements.
        process = {"reaction": "CO_g + *_t -> CO_t",
                   "description": "CO adsorption at top site.",
                   "coordinates_group": [[[0.0, 0.0, 0.0]]],
                   "elements_before": ["V", "V"],
                   "elements_after": ["C"],
                   "basis_sites": [0]}

        self.assertRaisesRegexp(SetupError, r"^Lengths of", check_process_dict, process)
        # }}}

    def test_check_process_coord(self):
        " Test process coordinates can be checked correctly. "

        c = [[0.0, 0.0, 0.0],
             [0.0, 0.5, 0.0],
             [0.5, 0.5, 0.0],
             [0.5, 0.0, 0.0],
             [0.5, -0.5, 0.0],
             [0.0, -0.5, 0.0],
             [-0.5, -0.5, 0.0],
             [-0.5, 0.0, 0.0],
             [-0.5, 0.5, 0.0]]

        check_process_coordinates(c)

        c = [[0.0, 0.0, 0.0],
             [0.0, 0.5, 0.0],
             [0.5, 0.5, 0.0],
             [0.5, 0.0, 0.0],
             [0.5, -0.5, 0.0],
             [0.0, -0.5, 0.0],
             [-0.5, -0.5, 0.0],
             [-0.5, 0.0, 0.0],
             [0.0, 0.0, 0.0]]

        self.assertRaisesRegexp(SetupError, "^Found", check_process_coordinates, c)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(CheckUtilitiesTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


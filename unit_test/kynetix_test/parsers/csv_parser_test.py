import logging
import unittest

from kynetix.model import KineticModel
from kynetix.parsers import *


class CsvParserTest(unittest.TestCase):

    def setUp(self):
        self.maxDiff = None

    def test_csv_parser_construction(self):
        " Test csv parser can be constructed. "
        # Construction.
        model = KineticModel(setup_file="input_files/csv_parser.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()

        # Check the parser class and base class type.
        self.assertTrue(isinstance(parser, CsvParser))
        self.assertEqual(parser.__class__.__base__.__name__, "ParserBase")

    def test_data_parse(self):
        " Test data in csv file can be read correctly. "
        # Construction.
        model = KineticModel(setup_file="input_files/csv_parser.mkm",
                             verbosity=logging.ERROR)
        parser = model.parser()

        ref_species_definitions = {'CO-O_s': {'DFT_energy': -115241.8617,
                                    'elements': {'C': 1, 'O': 2},
                                    'formation_energy': 4.2,
                                    'information': 'None',
                                    'site': 's',
                                    'site_number': 1,
                                    'type': 'transition_state'},
                                   'CO2_g': {'DFT_energy': -32.96253087,
                                    'elements': {'C': 1, 'O': 2},
                                    'formation_energy': 2.45,
                                    'information': 'None',
                                    'pressure': 0.0,
                                    'site': 'g',
                                    'site_number': 1,
                                    'type': 'gas'},
                                   'CO_g': {'DFT_energy': -626.6119705,
                                    'elements': {'C': 1, 'O': 1},
                                    'formation_energy': 2.74,
                                    'information': 'None',
                                    'pressure': 1.0,
                                    'site': 'g',
                                    'site_number': 1,
                                    'type': 'gas'},
                                   'CO_s': {'DFT_energy': -115390.4456,
                                    'elements': {'C': 1, 'O': 1},
                                    'formation_energy': 1.55,
                                    'information': 'None',
                                    'site': 's',
                                    'site_number': 1,
                                    'type': 'adsorbate'},
                                   'O2_g': {'DFT_energy': -496.4113942,
                                    'elements': {'O': 2},
                                    'formation_energy': 5.42,
                                    'information': 'None',
                                    'pressure': 0.3333333333333333,
                                    'site': 'g',
                                    'site_number': 1,
                                    'type': 'gas'},
                                   'O_s': {'DFT_energy': -115225.1065,
                                    'elements': {'O': 1},
                                    'formation_energy': 1.55,
                                    'information': 'None',
                                    'site': 's',
                                    'site_number': 1,
                                    'type': 'adsorbate'},
                                   's': {'DFT_energy': -114762.2548,
                                    'formation_energy': 0.0,
                                    'information': 'None',
                                    'site_name': '111',
                                    'total': 1.0,
                                    'type': 'site'}}

        parser.parse_data()
        self.assertDictEqual(ref_species_definitions, model.species_definitions())

    def test_get_barriers(self):
        " Test parsers can calculate reaction barriers correctly. "
        # Construction.
        model = KineticModel(setup_file="input_files/csv_parser.mkm",
                             verbosity=logging.ERROR)
        parser = model.parser()

        # Before get absolute data.
        rxn_list = model.elementary_rxns_list()[-1]
        self.assertRaises(AttributeError, parser.get_barrier, rxn_list)

        # Parse absolute data.
        parser.parse_data()
        ret_f_barrier, ret_r_barrier = parser.get_barrier(rxn_list)
        ref_f_barrier, ref_r_barrier = (1.1, 1.75)
        self.assertEqual(ref_f_barrier, ret_f_barrier)
        self.assertEqual(ref_r_barrier, ret_r_barrier)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(CsvParserTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


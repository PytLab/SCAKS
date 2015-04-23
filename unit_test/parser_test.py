import unittest
import sys

class TestParser(unittest.TestCase):

    def setUp(self):
        #create an instance of KineticModel
        sys.path.append('E:\BaiduYun\MyDocuments\ECUST+\Python\script')
        from kinetic import model
        self.km = model.KineticModel(setup_file='methanation.mkm')
        self.parser = self.km.parser

    def test_set_logger(self):
        "Make sure parser has a logger"
        self.assertTrue(hasattr(self.parser, 'logger'))

        #should raise an exception for an AttributeError
        self.assertRaises(AttributeError)

    def test_parse_site_expression(self):
        "Make sure site can be parsed successfully"
        site_expression = '3*_s'
        site_dict = self.parser.parse_site_expression(site_expression)
        target_dict = {'s': {'number': 3, 'type': 's'}}
        self.assertDictEqual(target_dict, site_dict)

    def test_parse_species_expression(self):
        "Make sure species expression e.g. '2CH3_s' can be parsed successfully"

        #test adsorbate
        adsorbate_expression = '2CH3_s'
        adsorbate_dict = \
                self.parser.parse_species_expression(adsorbate_expression)
        target_adsorbate_dict = {'CH3_s': {'number': 2, 
                                 'site': 's', 
                                 'elements': {'C': 1, 'H': 3}}}
        self.assertDictEqual(adsorbate_dict, target_adsorbate_dict)

        #test transition state
        ts_expression = '3CH2-H_s'
        ts_dict = self.parser.parse_species_expression(ts_expression)
        target_ts_dict = {'CH2-H_s': {'number': 3, 
                                      'site': 's', 
                                      'elements': {'C': 1, 'H': 3}}}
        self.assertDictEqual(ts_dict, target_ts_dict)

    def test_parse_state_expression(self):
        "Make sure state expression e.g. 'CH-H_s + *_s' can be parsed"

        #test IS or FS
        is_expression = 'CH-H_s + *_s'
        is_species_dict = self.parser.parse_state_expression(is_expression)[0]
        target_is_species_dict = \
           {'CH-H_s': {'elements': {'C': 1, 'H': 2}, 'number': 1, 'site': 's'}}
        self.assertDictEqual(is_species_dict, target_is_species_dict)

        is_sites_dict = self.parser.parse_state_expression(is_expression)[1]
        target_is_sites_dict = {'s': {'number': 1, 'type': 's'}}
        self.assertDictEqual(is_sites_dict, target_is_sites_dict)

        is_species_list = self.parser.parse_state_expression(is_expression)[2]
        target_is_species_list = ['CH-H_s', '*_s']
        self.assertListEqual(is_species_list, target_is_species_list)

        #test TS
        ts_expression = 'CH3-H_s + *_s'
        ts_species_dict = self.parser.parse_state_expression(ts_expression)[0]
        target_ts_species_dict = {'CH3-H_s': {'elements': {'C': 1, 'H': 4}, 
                                  'number': 1, 
                                  'site': 's'}}
        self.assertDictEqual(ts_species_dict, target_ts_species_dict)

        ts_sites_dict = self.parser.parse_state_expression(ts_expression)[1]
        target_ts_sites_dict = {'s': {'number': 1, 'type': 's'}}
        self.assertDictEqual(ts_sites_dict, target_ts_sites_dict)

        ts_species_list = self.parser.parse_state_expression(ts_expression)[2]
        target_ts_species_list = ['CH3-H_s', '*_s']
        self.assertListEqual(ts_species_list, target_ts_species_list)

    def test_parse_single_elementary_rxn(self):
        "Make sure one rxn equation can be parsed successfully"
        rxn_equation = 'O_s + H_s <-> H-O_s + *_s -> HO_s + 2*_s'
        target_states_dict = \
        {'FS': {'empty_sites_dict': {'s': {'number': 2, 'type': 's'}},
                'species_dict': {'HO_s': {'elements': {'H': 1, 'O': 1},
                                          'number': 1,
                                          'site': 's'}},
                'state_expression': 'HO_s + 2*_s'},
        'IS': {'empty_sites_dict': {},
               'species_dict': 
                   {'H_s': {'elements': {'H': 1}, 'number': 1, 'site': 's'},
                    'O_s': {'elements': {'O': 1}, 'number': 1, 'site': 's'}},
               'state_expression': 'O_s + H_s'},
        'TS': {'empty_sites_dict': {'s': {'number': 1, 'type': 's'}},
               'species_dict': {'H-O_s': {'elements': {'H': 1, 'O': 1},
                                          'number': 1,
                                          'site': 's'}},
               'state_expression': 'H-O_s + *_s'}}
        target_rxn_list = [['O_s', 'H_s'], ['H-O_s', '*_s'], ['HO_s', '2*_s']]
        state_dict, rxn_list = \
                        self.parser.parse_single_elementary_rxn(rxn_equation)

        self.assertDictEqual(target_states_dict, state_dict)
        self.assertListEqual(target_rxn_list, rxn_list)

    def test_parse_elementary_rxns(self):
        "Make sure all attrs of elementary rxns have been parsed into model"
        self.parser.parse_elementary_rxns(self.km.rxn_expressions)
        self.assertTrue(hasattr(self.km, 'adsorbate_names'))
        self.assertTrue(hasattr(self.km, 'gas_names'))
        self.assertTrue(hasattr(self.km, 'site_names'))
        self.assertTrue(hasattr(self.km, 'transition_state_names'))

    def test_parse_data(self):
        "Make sure data input by user is parsed into model"
        self.parser.parse_data()
        for species in self.parser._owner.species_definitions:
            for key in ['formation_energy', 'frequencies', 'information']:
                self.assertTrue(self.km.
                        species_definitions[species].has_key(key))

#if __name__ == '__main__':
#    unittest.main()
suite = unittest.TestLoader().loadTestsFromTestCase(TestParser)
unittest.TextTestRunner(verbosity=2).run(suite)
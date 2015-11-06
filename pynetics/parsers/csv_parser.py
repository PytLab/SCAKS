import logging

from parser_base import *


class CsvParser(ParserBase):
    def __init__(self, owner):
        ParserBase.__init__(self, owner)
        # set tools logger as child of model's
        self.logger = logging.getLogger('model.parsers.CsvParser')

    def parse_data(self):
        """
        Parse input file, add data
        e.g. formation energy, frequencies, information,
        to _owner.species_definition
        and set corresponding dict as model's attrs.
        """
        sp_dict = self._owner.species_definitions  # reference varibles
        file_obj = open('./energy.csv', 'rU')
        file_obj.readline()
        #formation_energy_dict, frequencies_dict = {}, {}
        for line_str in file_obj:
            line_list = line_str.strip('\n').split(',')
            site_name, species_name, information = \
                line_list[1], line_list[2], line_list[5]
            DFT_energy = float(line_list[3])
            formation_energy = float(line_list[4])

            #get full name
            #for gas
            if site_name == 'gas':
                full_name = species_name + '_g'
            elif site_name == 'liquid':
                full_name = species_name + '_l'
            #for any species on site
            else:
                #get site symbol of current line e.g. 's' or 'ss'
                for site_symbol in self._owner.site_names:
                    if sp_dict[site_symbol]['site_name'] == site_name:
                        break
                if species_name == 'slab':
                    full_name = site_symbol
                else:
                    full_name = species_name + '_' + site_symbol

            #update _owner.species_definition
            for key in ['DFT_energy', 'formation_energy', 'information']:
                #print eval(key)
                sp_dict[full_name].setdefault(key, eval(key))
        file_obj.close()

        setattr(self._owner, 'hasdata', True)

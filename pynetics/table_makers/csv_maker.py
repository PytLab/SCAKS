import string
import logging

from table_maker_base import *
from pynetics.functions import string2symbols


class CsvMaker(TableMakerBase):
    def __init__(self, owner):
        TableMakerBase.__init__(self, owner)
        self.logger = logging.getLogger('model.table_maker.CsvMaker')

    def create_initial_table(self):
        #create input file
        f = open('./energy.csv', 'w')
        f.write(self.add_lines())
        f.close()

    def add_lines(self):
        """
        Add lines to input files
        """
        gas_lines = ''
        liquid_lines = ''
        adsorbate_lines = ''
        transition_lines = ''
        site_lines = ''

        def update_line(line_tuple):
            row_list = [''] * 5
            #change row list
            for idx in range(5):
                row_list[idx] = line_tuple[idx]
            return ','.join(row_list) + '\n'

        #header
        header_str = ','.join(self.header) + '\n'
        #species part
        if self.species_definitions:
            for species in self.species_definitions:
                if self.species_definitions[species]['type'] == 'site':
                    surface_name = self.surface_name
                    site_name = self.species_definitions[species]['site_name']
                    species_name = 'slab'
                elif self.species_definitions[species]['type'] == 'gas':
                    surface_name = 'None'
                    site_name = 'gas'
                    species_name = self.species_definitions[species]['name']
                elif self.species_definitions[species]['type'] == 'liquid':
                    surface_name = 'None'
                    site_name = 'liquid'
                    species_name = self.species_definitions[species]['name']
                else:  # for adsorbate & transition state
                    surface_name = self.surface_name
                    site = self.species_definitions[species]['site']
                    site_name = self.species_definitions[site]['site_name']
                    species_name = self.species_definitions[species]['name']

                DFT_energy = '0.0'
                information = 'None'
                line_tuple = (surface_name, site_name, species_name,
                              DFT_energy, information)

                if self.species_definitions[species]['type'] == 'adsorbate':
                    adsorbate_lines += update_line(line_tuple)
                elif self.species_definitions[species]['type'] == 'transition_state':
                    transition_lines += update_line(line_tuple)
                elif self.species_definitions[species]['type'] == 'site':
                    site_lines += update_line(line_tuple)
                elif self.species_definitions[species]['type'] == 'gas':
                    gas_lines += update_line(line_tuple)
                elif self.species_definitions[species]['type'] == 'liquid':
                    liquid_lines += update_line(line_tuple)

        return header_str + gas_lines + liquid_lines + site_lines + \
            adsorbate_lines + transition_lines

    def get_new_row(self, row_str, mode):
        """
        Analyse each line of 'energy.csv',
        return a new line containing 'generalized formation energy'.
        """
        striped_str = string.whitespace
        row_list = row_str.strip(string.whitespace).split(',')
        site_name, species_name, energy = \
            row_list[1], row_list[2], float(row_list[3])
        if species_name == 'slab':
            element_list = []
        else:
            species_name = species_name.replace('-', '')
            element_list = string2symbols(species_name)

        if site_name == 'gas' or site_name == 'liquid':
            for element in element_list:
                energy -= self.ref_dict[element]
        else:  # for adsorbates, transition_state, slab
            for element in element_list + [site_name]:
                energy -= self.ref_dict[element]
        #generate a new line
        if mode == 'add':
            row_list.insert(4, str(round(energy, 3)))
        if mode == 'update':
            row_list[4] = str(round(energy, 3))
        new_row_str = ','.join(row_list)
        return new_row_str

    def create_new_table(self, mode):
        """
        Read initial input file, calculate generalized formation energy.
        Create a new input file containing
        a column of generalized formation energy.
        """
        #get old table content
        f = open('./energy.csv', 'rU')
        lines_list = f.readlines()
        f.close()

        content_str = ''

        #modify header
        header_str = lines_list[0].strip(string.whitespace)
        if mode == 'add':
            header_list = header_str.split(',')
            header_list.insert(4, 'formation_energy')
            new_header_str = ','.join(header_list) + '\n'
            content_str += new_header_str
        if mode == 'update':
            content_str += header_str + '\n'

        #modify species part
        for row_str in lines_list[1:]:
            new_row_str = self.get_new_row(row_str, mode=mode) + '\n'
            content_str += new_row_str

        #create new input file
        f = open('./energy.csv', 'w')
        f.write(content_str)
        f.close()

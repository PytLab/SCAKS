import cPickle

from functions import *


__version__ = '0.0.2'


class ModelShell(object):
    """
    A non-functional parent class to be inherited by
    other tools class of kinetic model.
    """
    def __init__(self, owner):
        self._owner = owner
        self.archived_data_dict = {}

    def split_species(self, species_str):
        "Split species_str to number(int) and species_name(str)"
        if not '*' in species_str:  # for species adsorbated on surface
            m = self._owner.regex_dict['species'][0].search(species_str)
            if not m.group(1):
                stoichiometry = 1
            else:
                stoichiometry = int(m.group(1))
            species_name = m.group(2) + '_' + m.group(3) + m.group(4)
            return stoichiometry, species_name
        else:  # for site
            m = self._owner.regex_dict['empty_site'][0].search(species_str)
            if not m.group(1):
                stoichiometry = 1
            else:
                stoichiometry = int(m.group(1))
            site_name = '*_' + m.group(2)
            return stoichiometry, site_name

    def update_defaults(self, defaults):
        """
        Update values in defaults dict,
        if there are custom parameters in setup file.
        """
        for parameter_name in defaults:
            if hasattr(self._owner, parameter_name):
                defaults[parameter_name] = \
                    getattr(self._owner, parameter_name)

        return defaults

    def archive_data(self, data_name, data):
        "Update data dict and dump it to data file."
        #update data dict
        if data_name in self.archived_variables:
            self.archived_data_dict[data_name] = data
            #dump data dict to data file
            if self.archived_data_dict:
                with open(self._owner.data_file, 'wb') as f:
                    cPickle.dump(self.archived_data_dict, f)

    @staticmethod
    def write2file(filename, line):
        f = open(filename, 'a')
        f.write(line)
        f.close()

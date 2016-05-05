from kynetix import ModelShell


class TableMakerBase(ModelShell):
    def __init__(self, owner):
        """
        Class for create raw data table for all species
        which will be used in kinetic model.
        This class acts as a base class to be inherited
        by other maker class.
        A table maker class should get all species in kinetic model
        to create the first column of table.
        """
        ModelShell.__init__(self, owner)
        self.header = ['surface_name', 'site_name', 'species_name',
                       'DFT_energy', 'infomation']

        #get suface names in setup file
        if hasattr(self._owner, 'surface_name'):
#            self.surfaces = list(getattr(self._owner, 'surface_name'))
            self.surface_name = getattr(self._owner, 'surface_name')
        else:
            self.logger.log(log_type='event',
                            event='attribute_load_failed',
                            attribute='surface_name')

        #species_definition
        if self._owner.species_definitions:
            self.species_definitions = self._owner.species_definitions

        #ref_dict
        if hasattr(self._owner, 'ref_dict'):
            self.ref_dict = self._owner.ref_dict

        self.row_list = [''] * len(self.header)

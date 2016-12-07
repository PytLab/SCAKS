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


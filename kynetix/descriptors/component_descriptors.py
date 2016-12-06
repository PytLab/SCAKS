from kynetix.parsers import *
from kynetix.solvers import *
from kynetix.table_makers import *
from kynetix.correctors import *
from kynetix.plotters import *
from kynetix.descriptors.descriptors import AttrDescriptor

class Component(AttrDescriptor):
    def __init__(self, name, cls_name, default, candidates):
        """
        Descriptor for kinetic model core components.
        """
        super(Component, self).__init__(name, cls_name, default)
        self.candidates = candidates

    def __set__(self, instance, value):
        # Instantialize the corresponding kinetic component.
        if type(value) is not str:
            raise ValueError("{} ({}) is not string".format(self.ori_name, value))
        if value not in self.candidates:
            raise ValueError("{} ({}) is not in {}".format(self.ori_name,
                                                           value,
                                                           self.candidates))
        if value in globals():
            parser_class = globals()[value]
            parser_instance = parser_class(owner=instance)
            instance.__dict__[self.name] = parser_instance
        else:
            raise ValueError("{} class is not found".format(value))


from ..parsers import *
from ..solvers import *
from ..correctors import *
from ..plotters import *
from .descriptors import AttrDescriptor

class Component(AttrDescriptor):
    def __init__(self, name, default, candidates):
        """
        Descriptor for kinetic model core components.
        """
        super(Component, self).__init__(name, default)
        self.candidates = candidates

    def __set__(self, instance, value):
        # Instantialize the corresponding kinetic component.
        if type(value) is not str:
            raise ValueError("{} ({}) is not string".format(self.name, value))
        if value not in self.candidates:
            raise ValueError("{} ({}) is not in {}".format(self.name,
                                                           value,
                                                           self.candidates))
        if value in globals():
            component_class = globals()[value]
            component_instance = component_class(owner=instance)
            private_name = "_{}__{}".format(instance.__class__.__name__, self.name)
            instance.__dict__[private_name] = component_instance
        else:
            raise ValueError("{} class is not found".format(value))


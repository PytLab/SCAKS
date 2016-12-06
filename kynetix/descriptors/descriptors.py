"""
Definitions of attribute descriptors.
"""
import numpy as np

class AttrDescriptor(object):
    def __init__(self, name, default):
        """
        Base descriptor class for the attributes.

        Parameters:
        -----------
        name: The attribute name, str.
        default: The default value when calling __get__() method,
                 could be any type.
        """
        # Get the mangled name.
        self.name = name
        self.default = default

    def __get__(self, instance, owner):
        private_name = "_{}__{}".format(instance.__class__.__name__, self.name)
        if private_name not in instance.__dict__:
            instance.__dict__[private_name] = self.default
        return instance.__dict__[private_name]

    def __set__(self, instance, value):
        private_name = "_{}__{}".format(instance.__class__.__name__, self.name)
        instance.__dict__[private_name] = value


class Type(AttrDescriptor):
    def __init__(self, name, type, default, candidates=None):
        """
        Descriptor for basic type attributes of kinetic model and its tools.

        Parameters:
        -----------
        name: The attribute name, str.
        default: The default value when calling __get__() method,
                 could be any type.
        candidates: All possible values of the attribute.
        """
        super(Type, self).__init__(name, default)
        self.candidates = candidates
        self.type = type

    def __set__(self, instance, value):
        # Check type.
        if type(value) is not self.type:
            msg = "{} ({}) is not a {}".format(self.name, value, self.type)
            raise ValueError(msg)
        # Check if in possible values.
        if self.candidates is not None and value not in self.candidates:
            msg = "{} ({}) is not one of {}".format(self.ori_name, value, self.candidates)
            raise ValueError(msg)
        # Set it.
        private_name = "_{}__{}".format(instance.__class__.__name__, self.name)
        instance.__dict__[private_name] = value


class Integer(Type):
    def __init__(self, name, default, candidates=None):
        super(Integer, self).__init__(name, int, default, candidates)


class Float(Type):
    def __init__(self, name, default, candidates=None):
        super(Float, self).__init__(name, float, default, candidates)

    def __set__(self, instance, value):
        # Overwrite father method.
        if type(value) not in (float, int):
            msg = "{} ({}) is not a {}".format(self.name, value, self.type)
            raise ValueError(msg)
        # Check if in possible values.
        if self.candidates is not None and value not in self.candidates:
            msg = "{} ({}) is not one of {}".format(self.name, value, self.candidates)
            raise ValueError(msg)
        # Set it.
        private_name = "_{}__{}".format(instance.__class__.__name__, self.name)
        instance.__dict__[private_name] = value

class String(Type):
    def __init__(self, name, default, candidates=None):
        super(String, self).__init__(name, str, default, candidates)


class Bool(Type):
    def __init__(self, name, default, candidates=None):
        super(Bool, self).__init__(name, bool, default, candidates)


class Dict(Type):
    def __init__(self, name, default, candidates=None):
        super(Dict, self).__init__(name, dict, default, candidates)


class Sequence(AttrDescriptor):
    def __init__(self, name, default, entry_type=None, candidates=None):
        """
        Descriptor for list type attributes of kinetic model and its tools.
        """
        super(Sequence, self).__init__(name, default)
        self.candidates = candidates
        self.entry_type = entry_type

    def __check(self, value):
        if type(value) not in (list, tuple):
            raise ValueError("{} ({}) is not a sequence".format(self.name, value))
        if self.entry_type is not None:
            for entry in value:
                if type(entry) is not self.entry_type:
                    msg = "{} in {} is not {}".format(entry, self.name, self.entry_type)
                    raise ValueError(msg)
        if self.candidates is not None:
            for entry in value:
                if entry not in self.candidates:
                    msg = "{} in {} is not one of {}".format(entry,
                                                             self.name,
                                                             self.candidates)
                    raise ValueError(msg)

    def __set__(self, instance, value):
        self.__check(value)
        # After check, set it.
        private_name = "_{}__{}".format(instance.__class__.__name__, self.name)
        instance.__dict__[private_name] = value


class FloatList2D(AttrDescriptor):
    def __init__(self, name, default):
        """
        Descriptor for 2D float list attributes of kinetic model.
        """
        super(FloatList2D, self).__init__(name, default)

    def _check(self, value):
        data_array = np.array(value)
        if type(value) is not list:
            raise ValueError("{} ({}) is not a list".format(self.name, value))
        if data_array.dtype.kind != 'f':
            raise ValueError("Types of entry is not float")
        if len(data_array.shape) != 2:
            raise ValueError("{} ({}) is not a 2d float list".format(self.name, value))

    def __set__(self, instance, value):
        self.__check(value)
        private_name = "_{}__{}".format(instance.__class__.__name__, self.name)
        instance.__dict__[private_name] = value


class SpaceVectors(FloatList2D):
    def __init__(self, name, default):
        """
        Descriptor for 3D space vectors list.
        """
        super(SpaceVectors, self).__init__(name, default)

    def __set__(self, instance, value):
        self._check(value)
        data_array = np.array(value)
        if data_array.shape[1] != 3:
            msg = "shape of {} ({}) is not (-1, 3)".format(self.name, value)
            raise ValueError(msg)
        private_name = "_{}__{}".format(instance.__class__.__name__, self.name)
        instance.__dict__[private_name] = value


class Property(object):
    def __init__(self, func):
        self.func = func

    def __get__(self, instance, cls):
        val = self.func(instance)
        attr_name = "_{}__{}".format(instance.__class__.__name__,
                                     self.func.__name__)
        if attr_name not in instance.__dict__:
            setattr(instance, attr_name, val)
        return val

    def __set__(self, instance, value):
        msg ="Changing value of {}.{} is not allowed".format(instance.__class__.__name__,
                                                             self.func.__name__)
        raise AttributeError(msg)


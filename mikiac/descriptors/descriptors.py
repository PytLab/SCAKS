"""
Definitions of attribute descriptors.
"""

import numpy as np
import copy

from ..utilities.check_utilities import check_species_definitions
from ..utilities.check_utilities import check_ref_energies
from ..utilities.check_utilities import check_analysis_interval


class AttrDescriptor(object):
    def __init__(self, name, default, deepcopy=False):
        """
        Base descriptor class for the attributes.

        Parameters:
        -----------
        name: The attribute name, str.
        default: The default value when calling __get__() method,
                 could be any type.
        deepcopy: The flag for returning the deepcopy or not, bool.
        """
        self.name = name
        self.default = default
        self.deepcopy = deepcopy

    def __get__(self, instance, owner):
        private_name = "_{}__{}".format(instance.__class__.__name__, self.name)
        if private_name not in instance.__dict__:
            instance.__dict__[private_name] = self.default
        if self.deepcopy:
            return copy.deepcopy(instance.__dict__[private_name])
        else:
            return instance.__dict__[private_name]

    def _check(self, value):
        # Placeholder for data set checking.
        pass

    def __set__(self, instance, value):
        # Check the value validity.
        self._check(value)
        # Data assignment.
        private_name = "_{}__{}".format(instance.__class__.__name__, self.name)
        if private_name not in instance.__dict__:
            instance.__dict__[private_name] = value
        else:
            msg ="Changing value of {}.{} is not allowed".format(instance.__class__.__name__,
                                                                 self.name)
            raise AttributeError(msg)


class Type(AttrDescriptor):
    def __init__(self, name, type, default, deepcopy=False, candidates=None):
        """
        Descriptor for basic type attributes of kinetic model and its tools.

        Parameters:
        -----------
        name: The attribute name, str.
        default: The default value when calling __get__() method,
                 could be any type.
        candidates: All possible values of the attribute.
        """
        super(Type, self).__init__(name, default, deepcopy)
        self.candidates = candidates
        self.type = type

    def _check(self, value):
        # Check type.
        if type(value) is not self.type:
            msg = "{} ({}) is not a {}".format(self.name, value, self.type)
            raise ValueError(msg)
        # Check if in possible values.
        if self.candidates is not None and value not in self.candidates:
            msg = "{} ({}) is not one of {}".format(self.name, value, self.candidates)
            raise ValueError(msg)


class Integer(Type):
    def __init__(self, name, default, deepcopy=False, candidates=None):
        super(Integer, self).__init__(name, int, default, deepcopy, candidates)


class Float(Type):
    def __init__(self, name, default, deepcopy=False, candidates=None):
        super(Float, self).__init__(name, float, default, deepcopy, candidates)

    def _check(self, value):
        # Overwrite father method.
        if type(value) not in (float, int):
            msg = "{} ({}) is not a {}".format(self.name, value, self.type)
            raise ValueError(msg)
        # Check if in possible values.
        if self.candidates is not None and value not in self.candidates:
            msg = "{} ({}) is not one of {}".format(self.name, value, self.candidates)
            raise ValueError(msg)


class String(Type):
    def __init__(self, name, default, deepcopy=False, candidates=None):
        super(String, self).__init__(name, str, default, deepcopy, candidates)


class Bool(Type):
    def __init__(self, name, default, deepcopy=False, candidates=None):
        super(Bool, self).__init__(name, bool, default, deepcopy, candidates)


class Dict(Type):
    def __init__(self, name, default, deepcopy=False, candidates=None):
        super(Dict, self).__init__(name, dict, default, deepcopy, candidates)


class SpeciesDefinitions(AttrDescriptor):
    def __init__(self, name, default, deepcopy=True):
        super(SpeciesDefinitions, self).__init__(name, default, deepcopy)

    def __set__(self, instance, value):
        check_species_definitions(value)
        super(SpeciesDefinitions, self).__set__(instance, value)


class RefEnergies(AttrDescriptor):
    def __init__(self, name, default, deepcopy=False):
        super(RefEnergies, self).__init__(name, default, deepcopy)

    def _check(self, value):
        check_ref_energies(value)


class AnalysisInterval(AttrDescriptor):
    def __init__(self, name, default, deepcopy=False):
        super(AnalysisInterval, self).__init__(name, default, deepcopy)

    def _check(self, value):
        check_analysis_interval(value)


class Sequence(AttrDescriptor):
    def __init__(self, name, default, deepcopy=False, entry_type=None, candidates=None):
        """
        Descriptor for list type attributes of kinetic model and its tools.
        """
        super(Sequence, self).__init__(name, default, deepcopy)
        self.candidates = candidates
        self.entry_type = entry_type

    def _check(self, value):
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


class FloatList2D(AttrDescriptor):
    def __init__(self, name, default, deepcopy=False):
        """
        Descriptor for 2D float list attributes of kinetic model.
        """
        super(FloatList2D, self).__init__(name, default, deepcopy)

    def _check(self, value):
        data_array = np.array(value)
        if type(value) is not list:
            raise ValueError("{} ({}) is not a list".format(self.name, value))
        if data_array.dtype.kind != 'f':
            raise ValueError("Types of entry is not float")
        if len(data_array.shape) != 2:
            raise ValueError("{} ({}) is not a 2d float list".format(self.name, value))


class SpaceVectors(FloatList2D):
    def __init__(self, name, default, deepcopy=False):
        """
        Descriptor for 3D space vectors list.
        """
        super(SpaceVectors, self).__init__(name, default, deepcopy)

    def _check(self, value):
        super(SpaceVectors, self)._check(value)
        data_array = np.array(value)
        if data_array.shape[1] != 3:
            msg = "shape of {} ({}) is not (-1, 3)".format(self.name, value)
            raise ValueError(msg)


class Property(object):
    def __init__(self, func):
        self.func = func

    def __get__(self, instance, cls):
        return self.func(instance)

    def __set__(self, instance, value):
        msg ="Changing value of {}.{} is not allowed".format(instance.__class__.__name__,
                                                             self.func.__name__)
        raise AttributeError(msg)


# ------------------------------------------------------------------
# Functions and classes for parameters and return value memoization.
# ------------------------------------------------------------------

class HashableDict(dict):
    # Override the __hash__ method of dict to make it hashable.
    def __hash__(self):
        # Make all keys hashable.
        keys = tuple(map(make_hashable, self.keys()))
        # Make all values hashable.
        values = tuple(map(make_hashable, self.values()))
        # Return the hash value of dict.
        hash_value = hash(tuple(sorted(zip(keys, values))))

        return hash_value


class HashableList(list):
    def __hash__(self):
        # Use the hash value of corresponding hashable tuple.
        return hash(tuple(map(make_hashable, self)))


def make_hashable(var):
    """
    Function to make a immutable variable hashable.
    """
    if type(var) in (list, tuple):
        return HashableList(var)
    elif type(var) is dict:
        return HashableDict(var)
    else:
        return var


class Memoized(object):
    """
    Descriptor for returned value memoization.
    """
    def __init__(self, func):
        self.func = func
        self.results = {}

    def __get__(self, instance, cls):
        self.instance = instance
        return self

    def __call__(self, **kwargs):
        # Make all arguments hashable.
        key = make_hashable(kwargs)

        # NOTE: if relative_energies is None, then we should
        #       search the model's relative energies.
        if "relative_energies" in key and key["relative_energies"] is None:
            key["relative_energies"] = make_hashable(self.instance._owner.relative_energies)

        try:
            return self.results[key]
        except KeyError:
            self.results[key] = self.func(self.instance, **kwargs)
            return self.results[key]


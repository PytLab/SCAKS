"""
Module for holding common checking utility functions.
"""

from itertools import combinations
import logging

from ..errors.error import *


def check_sequence(sequence, entry_type=None, param_name="Tested object"):
    """ Check the given object is a list or tuple of a type.

    :param sequence: The sequence object to test
    :type sequence: any

    :param entry_type: The type of entry of list
    :type entry_type: type

    :param param_name: The test parameter name
    :type param_name: str

    :return: The valid sequence object.
    """

    msg = "{} is not a sequence of {}".format(param_name, entry_type)

    if type(sequence) not in (tuple, list):
        raise SetupError(msg)

    for entry in sequence:
        if not isinstance(entry, entry_type):
            raise SetupError(msg)

    return sequence


def check_string(string, string_range=None, param_name="Tested object"):
    """ Check the string type and if it is in the given sequence.

    :param string: The test string
    :type string: str

    :param string_range: All possible string
    :type string_range: list or tuple

    :param param_name: The test parameter name
    :type param_name: str

    :return: The valid string object
    """
    if not isinstance(string, str):
        raise SetupError("{} is not a string.".format(param_name))

    if string not in string_range:
        msg = "{} must be one of {}".format(param_name, string_range)
        raise SetupError(msg)

    return string


def check_species_definitions(species_definitions):
    """ Check the species definitions in setup file.

    :param species_definitions: Species information
    :type species_definitions: dict

    :return: The valid species definitions
    :rtype: dict
    """
    # Check type.
    if not isinstance(species_definitions, dict):
        raise SetupError("species definitions is not a dict.")

    # Check data.
    for species, definitions in species_definitions.items():
        if species.endswith("_g") and "pressure" not in definitions:
            msg = "No pressure info for gas species {}.".format(species)
            raise SetupError(msg)

    return species_definitions

def check_ref_energies(ref_energies):
    """ Check the ref_energies in setup file.

    :param ref_energies: Reference energies
    :type ref_energies: type

    :return: The valid reference energies
    """
    # Check type.
    if not isinstance(ref_energies, dict):
        raise SetupError("ref_energies is not a dict.")

    for key, value in ref_energies.items():
        if not isinstance(key, str):
            raise SetupError("key '{}' is not a string.".format(key))
        if not isinstance(value, float):
            raise SetupError("value '{}' is not a float.".format(value))

    return ref_energies


def check_analysis_interval(analysis_interval):
    """ Check the analysis_interval in setup file.

    :return: The valid analysis interval
    """
    if type(analysis_interval) not in (int, list, tuple):
        raise SetupError("Invalid analysis_interval: int or list is expected.")

    if type(analysis_interval) in (list, tuple):
        for interval in analysis_interval:
            # Check type.
            if type(interval) not in (int, list, tuple):
                raise SetupError("Invalid analysis interval: {}".format(interval))

            # Check custom interval.
            elif type(interval) in (list, tuple):
                if len(interval) != 3:
                    raise SetupError("Length of interval {} not equal to 3".format(interval))
                else:
                    # Check type.
                    if not all([type(i) == int for i in interval]):
                        raise SetupError("Type of element in {} are all not int".format(interval))
                    # Check data validity.
                    else:
                        start, end, step = interval
                        if end < start:
                            raise SetupError("Start larger than end: {}".format(interval))
                        elif step > end - start:
                            raise SetupError("Step larger than span: {}".format(interval))

    # If all test passed, return.
    return analysis_interval


def check_process_dict(process_dict):
    """ Check if the process dict is correct.

    :return: The valid process dict
    :rtype: dict
    """
    # Check keys.
    essential_keys = ("reaction", "coordinates_group", "elements_before",
                      "elements_after", "basis_sites")
    for key in essential_keys:
        if key not in process_dict:
            msg = "key '{}' is not in process_dict {}".format(key, process_dict.keys())
            raise SetupError(msg)

    # Check reaction.
    if not isinstance(process_dict["reaction"], str):
        msg = "reaction must be a string."
        raise SetupError(msg)

    # Check type of group.
    check_sequence(process_dict["coordinates_group"],
                     entry_type=list,
                     param_name="coordinates_group")

    # Check each coordinates in group.
    for coords in process_dict["coordinates_group"]:
        check_sequence(coords,
                         entry_type=list,
                         param_name="coordinates")
        # Check the elements in coordinates.
        for coord in coords:
            check_sequence(coord, float, coord)
            if len(coord) != 3:
                msg = "{} must have 3 entries.".format(coord)
                raise SetupError(msg)

        # Check if coordinates are all different.
        check_process_coordinates(coords)

    # Check elements.
    for key in ["elements_before", "elements_after"]:
        check_sequence(process_dict[key], str, key)

        # Check length.
        for coords in process_dict["coordinates_group"]:
            if len(process_dict[key]) != len(coords):
                msg = "Lengths of {} and {} are different."
                msg = msg.format(process_dict[key], coords)
                raise SetupError(msg)

    # Check basis sites.
    check_sequence(process_dict["basis_sites"], int, "basis_site")

    # All tests passed, return.
    return process_dict


def check_process_coordinates(coordinates):
    """ Function to check all coordinates in a process are different.

    :param coordinates: A list of coordinates
    :type coordinates: 2D list or array of float

    :return: A valid coordinates.
    """
    coord_pairs = combinations(coordinates, 2)

    # Nested function to compare two coordinates.
    def equal(coord1, coord2):
        for e1, e2 in zip(coord1, coord2):
            if abs(e1 - e2) > 1e-6:
                return False
        return True

    # Check coordinates pairs.
    for c1, c2 in coord_pairs:
        if equal(c1, c2):
            msg = "Found equivalent coordinates: {} == {}".format(c1, c2)
            raise SetupError(msg)

    return coordinates


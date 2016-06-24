"""
Module for holding common checking utility functions.
"""

import logging

from kynetix.errors.error import *


def check_list_tuple(sequence, entry_type=None, param_name="Tested object"):
    """
    Check the given object is a list or tuple of a type.

    Parameters:
    -----------
    sequence: The object to test.

    entry_type: The type of entry of list.

    param_name: The test parameter name.

    Returns:
    --------
    sequence: The valid sequence object.
    """

    msg = "{} is not a sequence of {}".format(param_name, entry_type)

    if type(sequence) not in (tuple, list):
        raise SetupError(msg)

    for entry in sequence:
        if not isinstance(entry, entry_type):
            raise SetupError(msg)

    return sequence


def check_species_definitions(species_definitions):
    """
    Check the species definitions in setup file.

    Parameters:
    -----------
    species_definitions: A dict of species information.

    Return:
    -------
    species_definitions: The valid species definitions.
    """
    # Check type.
    if not isinstance(species_definitions, dict):
        raise SetupError("species definitions is not a dict.")

    # Check data.
    for species, definitions in species_definitions.iteritems():
        if species.endswith("_g") and "pressure" not in definitions:
            msg = "No pressure info for gas species {}.".format(species)
            raise SetupError(msg)

    return species_definitions

def check_ref_energies(ref_energies):
    """
    Check the ref_energies in setup file.

    Parameters:
    -----------
    ref_energies: A dict of reference energies.

    Return:
    -------
    ref_energies: The valid reference energies.
    """
    # Check type.
    if not isinstance(ref_energies, dict):
        raise SetupError("ref_energies is not a dict.")

    for key, value in ref_energies.iteritems():
        if not isinstance(key, str):
            raise SetupError("key '{}' is not a string.".format(key))
        if not isinstance(value, float):
            raise SetupError("value '{}' is not a float.".format(value))

    return ref_energies


def check_string(string, string_range=None, param_name="Tested object"):
    """
    Check the string type and if it is in the given sequence.

    Parameters:
    -----------
    string: The test string, str.

    string_range: All possible string, list or tuple.

    param_name: The test parameter name, str.

    Return:
    -------
    string: The valid string object.
    """
    if not isinstance(string, str):
        raise SetupError("{} is not a string.".format(param_name))

    if string not in string_range:
        msg = "{} must be one of {}".format(param_name, string_range)
        raise SetupError(msg)

    return string

def check_process_dict(process_dict):
    """
    Check if the process dict is correct.

    Return:
    -------
    The valid process dict.
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
    check_list_tuple(process_dict["coordinates_group"],
                     entry_type=list,
                     param_name="coordinates_group")

    # Check each coordinates in group.
    for coords in process_dict["coordinates_group"]:
        check_list_tuple(coords,
                         entry_type=list,
                         param_name="coordinates")
        # Check the elements in coordinates.
        for coord in coords:
            check_list_tuple(coord, float, coord)
            if len(coord) != 3:
                msg = "{} must have 3 entries.".format(coord)
                raise SetupError(msg)

    # Check elements.
    for key in ["elements_before", "elements_after"]:
        check_list_tuple(process_dict[key], str, key)

        # Check length.
        for coords in process_dict["coordinates_group"]:
            if len(process_dict[key]) != len(coords):
                msg = "Lengths of {} and {} are different."
                msg = msg.format(process_dict[key], coords)
                raise SetupError(msg)

    # Check basis sites.
    check_list_tuple(process_dict["basis_sites"], int, "basis_site")

    # All tests passed, return.
    return process_dict

table_maker_range = ("CsvMaker", )
parsers_range = ("RelativeEnergyParser", "CsvParser", "KMCParser")
solvers_range = ("KMCSolver", "SteadyStateSolver", "QuasiEquilibriumSolver")
corrector_range = ("ThermodynamicCorrector", )
plotter_range = ("EnergyProfilePlotter", )
rootfinding_range = ("ConstrainedNewton", "MDNewton")

type_rules = {
    "rxn_expressions": (check_list_tuple, str),
    "ref_energies": (check_ref_energies, ),
    "species_definitions": (check_species_definitions, ),
    "surface_name": (str, ),
    "temperature": (int, ),
    "tools": (check_list_tuple, str),
    "parser": (check_string, parsers_range),
    "solver": (check_string, solvers_range),
    "corrector": (check_string, corrector_range),
    "plotter": (check_string, plotter_range),
    "ref_species": (check_list_tuple, str),
    "numerical_representation": (str, ),
    "archived_variables": (check_list_tuple, str),
    "rootfinding": (check_string, rootfinding_range),
    "max_rootfinding_iterations": (int, ),
    "residual_threshold": (float, ),
    "perturbation_size": (float, ),
    "initial_guess_scale_factor": (float, ),
    "stable_criterion": (float, ),
    "perturbation_direction": (str, ),
    "tolerance": (float, ),
    "grid_type": (str, ),
    "gas_thermo_mode": (str, ),
    "decimal_precision": (int, ),
    "data_file": (str, ),
    "table_maker": (check_string, table_maker_range),
    "ref_energies": (dict, ),

    # KMC parameters.
    "cell_vectors": (check_list_tuple, list),
    "basis_sites": (check_list_tuple, list),
    "unitcell_area": (float, ),
    "active_ratio": (float, ),
    "repetitions": (check_list_tuple, int),
    "periodic": (check_list_tuple, bool),
    "nstep": (int, ),
    "seed": (int, ),
    "random_generator": (str, ),
    "analysis": (check_list_tuple, str),
    "analysis_interval": (int, ),
    "trajectory_dump_interval": (int, ),
    "color_dict": (dict, ),
    "circle_attrs": (dict, ),
    "possible_element_types": (str, ),
    "possible_site_types": (str, ),
}


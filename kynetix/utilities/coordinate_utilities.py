"""
Module for doing coordinates group operation.
"""

import copy

import numpy as np


class CoordsGroup(object):

    # Static class varibles.
    tolerance = 1e-4

    def __init__(self, coordinates=None, elements=None):
        """
        Constructor.

        Parameters:
        -----------
        coordinates: A list of coordinates, default value is [].

        elements: A list of elements, list of str,
                  default value is ["V", ...].
        """
        if coordinates is None:
            self.__coords = []
        else:
            self.__coords = coordinates

        if elements is None:
            self.__elements = ["V"]*len(self.__coords)
        else:
            self.__elements = elements

        # Check.
        if len(self.__coords) != len(self.__elements):
            msg = "Length of coordinates({}) and elements({}) are different."
            msg = msg.format(len(self.__coords), len(self.__elements))
            raise ValueError(msg)

    def __check_coordinate(self, coordinate):
        """
        Private helper function to check validty of coordinate.
        """
        msg = "coordinate must be a list of three float numbers."
        if type(coordinate) is not list or len(coordinate) != 3:
            raise ValueError(msg)

        for entry in coordinate:
            if type(entry) is not float:
                raise ValueError(msg)

        return coordinate

    def append(self, coordinate, element=None):
        """
        Function to append a new coordinate.

        Parameters:
        -----------
        coordinate: A list of 3 float.

        element: str, default value is "V".
        """
        # Append coordinates.
        self.__check_coordinate(coordinate)
        self.__coords.append(coordinate)

        # Append elements.
        if element is None:
            self.__elements.append("V")
        else:
            if type(element) is not str:
                raise ValueError("element must be a string.")
            self.__elements.append(element)

    def extend(self, coordinates, elements=None):
        """
        Function to add new coordinates.

        Parameters:
        -----------
        coordinates: A list of coordinates.

        elements: A list of elements, list of str, default value is ["V", ...].
        """
        if elements is None:
            elements = ["V"]*len(coordinates)

        for coordinate, element in zip(coordinates, elements):
            self.append(coordinate, element)

    def move(self, move_vector=None):
        """
        Function to move all coordinates.

        Parameters:
        -----------
        move_vector: Move vector, list of 3 floats, default value is [0, 0, 0]

        Returns:
        --------
        A new moved CoordsGroup object.
        """
        if move_vector is None:
            move_vector = [0.0, 0.0, 0.0]

        # Check the parameter.
        if type(move_vector) is not list and len(move_vector) != 3:
            raise ValueError("move_vector must be a list of three float numbers.")

        move_vector = np.array(move_vector)
        ori_coords = np.array(self.__coords)
        new_coords = (ori_coords + move_vector).tolist()

        return CoordsGroup(new_coords, copy.copy(self.__elements))

    @staticmethod
    def __compare_coords(coord1, coord2):
        """
        Private helper function to compare two coordinates.
        """
        for e1, e2 in zip(coord1, coord2):
            if abs(e1 - e2) > CoordsGroup.tolerance:
                return False
        return True

    def __add__(self, another):
        """
        Operator `+` overload function.
        """
        coords1, coords2 = self.__coords, another.coordinates()
        elems1, elems2 = self.__elements, another.elements()

        merged_coords = copy.copy(coords1)
        merged_elems = copy.copy(elems1)

        # Loop to add different entries.
        for coord, elem in zip(coords2, elems2):
            # Loop to check if the coord is in the first coords.
            in_it = False
            for idx, (ref_coord, ref_elem) in enumerate(zip(coords1, elems1)):
                coords_equal = CoordsGroup.__compare_coords(coord, ref_coord)

                if coords_equal:
                    # Check element, replace or not.
                    if elem != ref_elem:
                        if ref_elem == "V":
                            merged_elems[idx] = elem

                    in_it = True
                    break

            if not in_it:
                merged_coords.append(coord)
                merged_elems.append(elem)

        merged_coords_group = CoordsGroup(merged_coords, merged_elems)

        return merged_coords_group

    def __len__(self):
        return len(self.__coords)

    def coordinates(self):
        """
        Query function to get coordinates.
        """
        return self.__coords

    def elements(self):
        """
        Query function for elements.
        """
        return self.__elements


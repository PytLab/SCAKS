"""
Module for doing coordinates group operation.
"""

import copy


class CoordsGroup(object):

    # Static class varibles.
    tolerance = 1e-4

    def __init__(self, coordinates=None):
        if coordinates is None:
            self.__coords = []
        else:
            self.__coords = coordinates

    def append(self, coordinate):
        """
        Function to append a new coordinate.
        """
        # Check the parameter.
        if type(coordinate) is not list and len(coordinate) != 3:
            raise ValueError("coordinate must be a list of three float numbers.")

        self.__coords.append(coordinate)

    def extend(self, coordinates):
        """
        Function to add new coordinates.
        """
        for coordinate in coordinates:
            self.append(coordinate)

    def move(self, move_vector):
        """
        Function to move all coordinates.
        """
        # Check the parameter.
        if type(move_vector) is not list and len(move_vector) != 3:
            raise ValueError("move_vector must be a list of three float numbers.")

        move_vector = np.array(move_vector)
        ori_coords = np.array(self.__coords)
        new_coords = (ori_coords + move_vector).tolist()

        self.__coords = new_coords

        return new_coord

    @staticmethod
    def __compare_coords(coord1, coord2):
        """
        Private helper function to compare two coordinates.
        """
        for e1, e2 in zip(coord1, coord2):
            if abs(e1 - e2) > CoordsGroup.tolerance:
                return False
        return True

    @staticmethod
    def __merge_coordinates(coords1, coords2):
        """
        Private helper function to merge tow coordinates groups.
        """
        merged_coords = copy.copy(coords1)

        for coord in coords2:
            in_it = False
            for ref_coord in coords1:
                coords_equal = CoordsGroup.__compare_coords(coord, ref_coord)
                if coords_equal:
                    in_it = True
                    break
            if not in_it:
                merged_coords.append(coord)

        return merged_coords

    def __add__(self, another):
        """
        Operator `+` overload function.
        """
        merged_coords = self.__merge_coordinates(self.__coords,
                                                 another.coordinates())
        coords_group = CoordsGroup(merged_coords)

        return coords_group

    def __len__(self):
        return len(self.__coords)

    def coordinates(self):
        """
        Query function to get coordinates.
        """
        return self.__coords


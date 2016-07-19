import logging
import re
import unittest

from kynetix.errors.error import *
from kynetix.utilities.coordinate_utilities import *


class CoordinatesUtilitiesTest(unittest.TestCase):

    def setUp(self):
        # Test case setting.
        self.maxDiff = None

    def test_construction(self):
        " Test CoordsGroup object can be constructed correctly. "
        fcc_coords = [[0.0, 0.0, 0.0],
                      [1.0, 0.0, 0.0],
                      [0.0, 1.0, 0.0],
                      [0.0, 0.5, 0.0],
                      [0.5, 0.0, 0.0],
                      [0.5, 0.5, 0.0],
                      [0.3333333333333333, 0.3333333333333333, 0.0]]

        fcc = CoordsGroup(fcc_coords)

        self.assertListEqual(fcc_coords, fcc.coordinates())
        self.assertEqual(7, len(fcc))

        # Check with elements.
        elements = ['V', 'V', 'V', 'V', 'V', 'V', 'C']
        fcc = CoordsGroup(fcc_coords, elements)

        self.assertListEqual(fcc_coords, fcc.coordinates())
        self.assertListEqual(elements, fcc.elements())
        self.assertEqual(7, len(fcc))

    def test_append(self):
        " Test append() function can work correctly. "
        fcc_coords = [[0.0, 0.0, 0.0],
                      [1.0, 0.0, 0.0],
                      [0.0, 1.0, 0.0],
                      [0.0, 0.5, 0.0],
                      [0.5, 0.0, 0.0],
                      [0.5, 0.5, 0.0],
                      [0.3333333333333333, 0.3333333333333333, 0.0]]

        fcc = CoordsGroup(fcc_coords)
        fcc.append([0.0, 0.7, 0.7])

        ref_coords = [[0.0, 0.0, 0.0],
                      [1.0, 0.0, 0.0],
                      [0.0, 1.0, 0.0],
                      [0.0, 0.5, 0.0],
                      [0.5, 0.0, 0.0],
                      [0.5, 0.5, 0.0],
                      [0.3333333333333333, 0.3333333333333333, 0.0],
                      [0.0, 0.7, 0.7]]
        ret_coords = fcc.coordinates()
        self.assertListEqual(ref_coords, ret_coords)

    def test_append_with_element(self):
        # Check with element.
        fcc_coords = [[0.0, 0.0, 0.0],
                      [1.0, 0.0, 0.0],
                      [0.0, 1.0, 0.0],
                      [0.0, 0.5, 0.0],
                      [0.5, 0.0, 0.0],
                      [0.5, 0.5, 0.0],
                      [0.3333333333333333, 0.3333333333333333, 0.0]]

        fcc = CoordsGroup(fcc_coords)
        fcc.append([0.0, 0.7, 0.7], "C")

        ret_coords = fcc.coordinates()
        ref_coords = [[0.0, 0.0, 0.0],
                      [1.0, 0.0, 0.0],
                      [0.0, 1.0, 0.0],
                      [0.0, 0.5, 0.0],
                      [0.5, 0.0, 0.0],
                      [0.5, 0.5, 0.0],
                      [0.3333333333333333, 0.3333333333333333, 0.0],
                      [0.0, 0.7, 0.7]]
        self.assertListEqual(ref_coords, ret_coords)

        ref_elements = ['V', 'V', 'V', 'V', 'V', 'V', 'V', 'C']
        ret_elements = fcc.elements()
        self.assertListEqual(ref_elements, ret_elements)

    def test_extend(self):
        " Test function extend() can work correctly. "
        fcc_coords = [[0.0, 0.0, 0.0],
                      [1.0, 0.0, 0.0],
                      [0.0, 1.0, 0.0],
                      [0.0, 0.5, 0.0],
                      [0.5, 0.0, 0.0],
                      [0.5, 0.5, 0.0],
                      [0.3333333333333333, 0.3333333333333333, 0.0]]

        fcc = CoordsGroup(fcc_coords)
        fcc.extend([[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]])

        ref_coords = [[0.0, 0.0, 0.0],
                      [1.0, 0.0, 0.0],
                      [0.0, 1.0, 0.0],
                      [0.0, 0.5, 0.0],
                      [0.5, 0.0, 0.0],
                      [0.5, 0.5, 0.0],
                      [0.3333333333333333, 0.3333333333333333, 0.0],
                      [1.0, 1.0, 1.0],
                      [1.0, 1.0, 1.0]]

        ret_coords = fcc.coordinates()

        self.assertListEqual(ref_coords, ret_coords)

    def test_extend_with_elements(self):
        " Test function extend() can work correctly with elements. "
        fcc_coords = [[0.0, 0.0, 0.0],
                      [1.0, 0.0, 0.0],
                      [0.0, 1.0, 0.0],
                      [0.0, 0.5, 0.0],
                      [0.5, 0.0, 0.0],
                      [0.5, 0.5, 0.0],
                      [0.3333333333333333, 0.3333333333333333, 0.0]]

        fcc = CoordsGroup(fcc_coords)
        fcc.extend([[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]], ["C", "C"])

        ref_coords = [[0.0, 0.0, 0.0],
                      [1.0, 0.0, 0.0],
                      [0.0, 1.0, 0.0],
                      [0.0, 0.5, 0.0],
                      [0.5, 0.0, 0.0],
                      [0.5, 0.5, 0.0],
                      [0.3333333333333333, 0.3333333333333333, 0.0],
                      [1.0, 1.0, 1.0],
                      [1.0, 1.0, 1.0]]
        ret_coords = fcc.coordinates()
        self.assertListEqual(ref_coords, ret_coords)

        ref_elements = ['V', 'V', 'V', 'V', 'V', 'V', 'V', 'C', 'C'] 
        ret_elements = fcc.elements()
        self.assertListEqual(ref_elements, ret_elements)

    def test_move(self):
        " Make sure coordinates can move correctly."
        fcc_coords = [[0.0, 0.0, 0.0],
                      [1.0, 0.0, 0.0],
                      [0.0, 1.0, 0.0],
                      [0.0, 0.5, 0.0],
                      [0.5, 0.0, 0.0],
                      [0.5, 0.5, 0.0],
                      [0.3333333333333333, 0.3333333333333333, 0.0]]

        fcc = CoordsGroup(fcc_coords)
        ret_coords = fcc.move([0.0, 0.0, 1.0]).coordinates()

        ref_coords = [[0.0, 0.0, 1.0],
                      [1.0, 0.0, 1.0],
                      [0.0, 1.0, 1.0],
                      [0.0, 0.5, 1.0],
                      [0.5, 0.0, 1.0],
                      [0.5, 0.5, 1.0],
                      [0.3333333333333333, 0.3333333333333333, 1.0]]

        self.assertListEqual(ref_coords, ret_coords)

    def test_add(self):
        " Test operator add. "
        fcc_coords = [[0.0, 0.0, 0.0],
                      [1.0, 0.0, 0.0],
                      [0.0, 1.0, 0.0],
                      [0.0, 0.5, 0.0],
                      [0.5, 0.0, 0.0],
                      [0.5, 0.5, 0.0],
                      [0.3333333333333333, 0.3333333333333333, 0.0]]
        fcc_elems = ['V', 'V', 'V', 'V', 'V', 'V', 'C']

        fcc = CoordsGroup(fcc_coords, fcc_elems)

        top_coords = [[0.0, 0.0, 0.0],
                      [0.0, 0.5, 0.0],
                      [0.3333333333333333, 0.3333333333333333, 0.0],
                      [0.5, 0.0, 0.0],
                      [0.6666666666666666, -0.3333333333333333, 0.0],
                      [0.5, -0.5, 0.0],
                      [0.3333333333333333, -0.6666666666666666, 0.0],
                      [0.0, -0.5, 0.0],
                      [-0.3333333333333333, -0.3333333333333333, 0.0],
                      [-0.5, 0.0, 0.0],
                      [-0.6666666666666666, 0.3333333333333333, 0.0],
                      [-0.5, 0.5, 0.0],
                      [-0.3333333333333333, 0.6666666666666666, 0.0]]
        top_elems = ['C', 'V', 'V', 'V', 'V', 'V', 'V', 'V', 'V', 'V', 'V', 'V', 'V']
        top = CoordsGroup(top_coords, top_elems)

        merge = top + fcc
        ref_coords = [[0.0, 0.0, 0.0],
                      [0.0, 0.5, 0.0],
                      [0.3333333333333333, 0.3333333333333333, 0.0],
                      [0.5, 0.0, 0.0],
                      [0.6666666666666666, -0.3333333333333333, 0.0],
                      [0.5, -0.5, 0.0],
                      [0.3333333333333333, -0.6666666666666666, 0.0],
                      [0.0, -0.5, 0.0],
                      [-0.3333333333333333, -0.3333333333333333, 0.0],
                      [-0.5, 0.0, 0.0],
                      [-0.6666666666666666, 0.3333333333333333, 0.0],
                      [-0.5, 0.5, 0.0],
                      [-0.3333333333333333, 0.6666666666666666, 0.0],
                      [1.0, 0.0, 0.0],
                      [0.0, 1.0, 0.0],
                      [0.5, 0.5, 0.0]]
        ret_coords = merge.coordinates()
        self.assertListEqual(ref_coords, ret_coords)

        ref_elems = ['C', 'V', 'C', 'V', 'V', 'V', 'V', 'V', 'V', 'V', 'V', 'V', 'V', 'V', 'V', 'V']
        ret_elems = merge.elements()
        self.assertListEqual(ref_elems, ret_elems)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(CheckUtilitiesTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


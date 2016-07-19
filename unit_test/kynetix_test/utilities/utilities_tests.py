import unittest

from check_utilities_test import CheckUtilitiesTest
from coordinates_utilities_test import CoordinatesUtilitiesTest

def suite():
    suite = unittest.TestSuite(
        [unittest.TestLoader().loadTestsFromTestCase(CheckUtilitiesTest),
         unittest.TestLoader().loadTestsFromTestCase(CoordinatesUtilitiesTest),
         ])
    return suite


if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())


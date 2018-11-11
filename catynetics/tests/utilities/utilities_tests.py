import unittest

from .coordinates_utilities_test import CoordinatesUtilitiesTest

util_test_cases = [CoordinatesUtilitiesTest]

def suite():
    suite = unittest.TestSuite(
        unittest.TestLoader().loadTestsFromTestCase(tc) for tc in util_test_cases
    )
    return suite


if __name__ == '__main__':
    result = unittest.TextTestRunner(verbosity=2).run(suite())

    if result.errors or result.failures:
        raise Exception('Get errors and failures')


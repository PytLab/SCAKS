import unittest

from .thermodynamic_corrector_test import ThermodynamicCorrectorTest

corrector_test_cases = [ThermodynamicCorrectorTest]

def suite():
    suite = unittest.TestSuite(
        unittest.TestLoader().loadTestsFromTestCase(tc) for tc in corrector_test_cases
    )
    return suite


if __name__ == '__main__':
    result = unittest.TextTestRunner(verbosity=2).run(suite())

    if result.errors or result.failures:
        raise Exception('Get errors and failures')


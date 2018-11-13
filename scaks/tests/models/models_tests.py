import unittest

from .mkm_model_test import MicroKineticModelTest

model_test_cases = [
    MicroKineticModelTest,
]

def suite():
    suite = unittest.TestSuite([
        unittest.TestLoader().loadTestsFromTestCase(tc) for tc in model_test_cases
    ])
    
    return suite

if __name__ == '__main__':
    result = unittest.TextTestRunner(verbosity=2).run(suite())

    if result.errors or result.failures:
        raise Exception('Get errors and failures')


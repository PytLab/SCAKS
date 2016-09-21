import unittest

from mkm_model_test import MicroKineticModelTest
from kmc_model_test import KMCModelTest


def suite():
    suite = unittest.TestSuite(
        [unittest.TestLoader().loadTestsFromTestCase(MicroKineticModelTest),        
         unittest.TestLoader().loadTestsFromTestCase(KMCModelTest)]        
    )
    
    return suite

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())


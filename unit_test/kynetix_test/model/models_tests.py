import unittest

from mkm_model_test import MicroKineticModelTest


def suite():
    suite = unittest.TestSuite(
        [unittest.TestLoader().loadTestsFromTestCase(MicroKineticModelTest)]        
    )
    
    return suite

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())


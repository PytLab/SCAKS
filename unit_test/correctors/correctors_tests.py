import unittest

from thermodynamic_corrector_test import ThermodynamicCorrectorTest


def suite():
    suite = unittest.TestSuite(
         [unittest.TestLoader().loadTestsFromTestCase(ThermodynamicCorrectorTest),]
    )
    return suite


if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())


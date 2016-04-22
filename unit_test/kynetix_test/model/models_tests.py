import unittest

from relative_energy_parser_test import RelativeEnergyParserTest


def suite():
    suite = unittest.TestSuite(
        [unittest.TestLoader().loadTestsFromTestCase(RelativeEnergyParserTest),
         ])
    return suite


if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())

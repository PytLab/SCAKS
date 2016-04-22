import unittest

from .parsers import relative_energy_parser_test


def suite():
    suite = unittest.TestSuite(
        [relative_energy_parser_test.suite()])
    return suite


if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())

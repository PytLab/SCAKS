import unittest

from parser_base_test import ParserBaseTest


def suite():
    suite = unittest.TestSuite(
        [unittest.TestLoader().loadTestsFromTestCase(ParserBaseTest)]
    )
    return suite


if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())


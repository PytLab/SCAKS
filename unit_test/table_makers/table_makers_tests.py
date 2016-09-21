import unittest

from csv_maker_test import CsvMakerTest


def suite():
    suite = unittest.TestSuite(
        [unittest.TestLoader().loadTestsFromTestCase(CsvMakerTest)]
    )
    return suite


if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())


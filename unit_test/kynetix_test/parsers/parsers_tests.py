import unittest

from parser_base_test import ParserBaseTest
from csv_parser_test import CsvParserTest
from relative_energy_parser_test import RelativeEnergyParserTest


def suite():
    suite = unittest.TestSuite(
        [unittest.TestLoader().loadTestsFromTestCase(ParserBaseTest),
         unittest.TestLoader().loadTestsFromTestCase(CsvParserTest),
         unittest.TestLoader().loadTestsFromTestCase(RelativeEnergyParserTest),]
    )
    return suite


if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())


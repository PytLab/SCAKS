import unittest

from .chem_formula_test import ChemFormulaTest
from .chem_state_test import ChemStateTest
from .parser_base_test import ParserBaseTest
from .relative_energy_parser_test import RelativeEnergyParserTest
from .absolute_energy_parser_test import AbsoluteEnergyParserTest
from .kmc_parser_test import KMCParserTest

parser_test_cases = [
    ChemStateTest,
    ChemFormulaTest,
    ParserBaseTest,
    RelativeEnergyParserTest,
    AbsoluteEnergyParserTest,
    KMCParserTest
]


def suite():
    suite = unittest.TestSuite([
        unittest.TestLoader().loadTestsFromTestCase(tc) for tc in parser_test_cases
    ])
    return suite


if __name__ == '__main__':
    result = unittest.TextTestRunner(verbosity=2).run(suite())

    if result.errors or result.failures:
        raise Exception('Get errors and failures')


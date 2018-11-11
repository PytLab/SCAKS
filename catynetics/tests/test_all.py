"""
Module for all test suites for MiKiAC.
"""
import unittest

from ..compatutil import subprocess
from .correctors.correctors_tests import *
from .models.models_tests import *
from .plotters.plotters_tests import *
from .parsers.parsers_tests import *
from .solvers.solvers_tests import *
from .utilities.utilities_tests import *

from . import cleanup

test_cases = (model_test_cases +
              parser_test_cases +
              solver_test_cases +
              plotter_test_cases +
              corrector_test_cases +
              util_test_cases)

def suite():
    """
    Function to get all suites for catynetics.
    """
    suite = unittest.TestSuite([
        unittest.TestLoader().loadTestsFromTestCase(tc) for tc in test_cases
    ])

    return suite


if __name__ == '__main__':
    result = unittest.TextTestRunner(verbosity=2).run(suite())

    cleanup()

    if result.errors or result.failures:
        raise Exception('Get errors and failures')

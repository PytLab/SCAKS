import unittest

from model import models_tests
from parsers import parsers_tests
from solvers import solvers_tests
from utilities import utilities_tests


def suite():
    suite = unittest.TestSuite(
        [models_tests.suite(),
         parsers_tests.suite(),
         solvers_tests.suite(),
         utilities_tests.suite()]
    )
    return suite


if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())


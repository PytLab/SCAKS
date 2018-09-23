"""
Module for all test suites for MiKiAC.
"""
import unittest

from mikiac.compatutil import subprocess
from correctors import correctors_tests
from models import models_tests
from plotters import plotters_tests
from parsers import parsers_tests
from solvers import solvers_tests
from utilities import utilities_tests


def suite():
    """
    Function to get all suites for mikiac.
    """
    all_suites = unittest.TestSuite([models_tests.suite(),
                                     parsers_tests.suite(),
                                     solvers_tests.suite(),
                                     utilities_tests.suite(),
                                     plotters_tests.suite(),
                                     correctors_tests.suite()])
    return all_suites


if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())

    # Remove auto-generated files.
    subprocess.getstatusoutput("for i in `find ./ -name 'auto_*'`; do rm -rf $i; done")
    subprocess.getstatusoutput("for i in `find ./ -name 'out.log'`; do rm -rf $i; done")
    subprocess.getstatusoutput("for i in `find ./ -name 'log'`; do rm -rf $i; done")
    subprocess.getstatusoutput("for i in `find ./ -name '*.pkl'`; do rm -rf $i; done")
    subprocess.getstatusoutput("rm ./*.csv")


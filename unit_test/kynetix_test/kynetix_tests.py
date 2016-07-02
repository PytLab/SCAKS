import commands
import unittest

from correctors import correctors_tests
from model import models_tests
from plotters import plotters_tests
from parsers import parsers_tests
from solvers import solvers_tests
from table_makers import table_makers_tests
from utilities import utilities_tests


def suite():
    suite = unittest.TestSuite([models_tests.suite(),
                                parsers_tests.suite(),
                                solvers_tests.suite(),
                                utilities_tests.suite(),
                                plotters_tests.suite(),
                                correctors_tests.suite(),
                                table_makers_tests.suite()])
    return suite


if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())

    # Remove auto-generated files.
    commands.getstatusoutput("for i in `find ./ -name 'auto_*'`; do rm -rf $i; done")
    commands.getstatusoutput("for i in `find ./ -name 'out.log'`; do rm -rf $i; done")
    commands.getstatusoutput("for i in `find ./ -name '*.pkl'`; do rm -rf $i; done")


import unittest

from solver_base_test import SolverBaseTest
from steady_state_solver_test import SteadyStateSolverTest


def suite():
    suite = unittest.TestSuite(
        [unittest.TestLoader().loadTestsFromTestCase(SolverBaseTest),
         unittest.TestLoader().loadTestsFromTestCase(SteadyStateSolverTest),]
    )
    return suite


if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())


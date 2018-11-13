import unittest

from .mean_field_solver_test import MeanFieldSolverTest
from .solver_base_test import SolverBaseTest
from .steady_state_solver_test import SteadyStateSolverTest

solver_test_cases = [
    MeanFieldSolverTest,
    SolverBaseTest,
    SteadyStateSolverTest,
]

def suite():
    suite = unittest.TestSuite([
        unittest.TestLoader().loadTestsFromTestCase(tc) for tc in solver_test_cases
    ])
    return suite


if __name__ == '__main__':
    result = unittest.TextTestRunner(verbosity=2).run(suite())

    if result.errors or result.failures:
        raise Exception('Get errors and failures')


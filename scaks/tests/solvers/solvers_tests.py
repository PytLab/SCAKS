import unittest

from .kmc_solver_test import KMCSolverTest
from .kmc_coverages_plugin_test import KMCCoveragesPluginTest
from .kmc_frequency_plugin_test import KMCFrequencyPluginTest
from .kmc_tof_plugin_test import KMCTOFPluginTest
from .mean_field_solver_test import MeanFieldSolverTest
from .solver_base_test import SolverBaseTest
from .steady_state_solver_test import SteadyStateSolverTest
from .kmc_redistribution_test import KMCRedistributionTest

solver_test_cases = [
    KMCSolverTest,
    KMCCoveragesPluginTest,
    KMCFrequencyPluginTest,
    KMCTOFPluginTest,
    MeanFieldSolverTest,
    SolverBaseTest,
    SteadyStateSolverTest,
    KMCRedistributionTest
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


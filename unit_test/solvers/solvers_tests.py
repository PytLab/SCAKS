import unittest

from .kmc_solver_test import KMCSolverTest
from .kmc_coverages_plugin_test import KMCCoveragesPluginTest
from .kmc_frequency_plugin_test import KMCFrequencyPluginTest
from .kmc_tof_plugin_test import KMCTOFPluginTest
from .mean_field_solver_test import MeanFieldSolverTest
from .solver_base_test import SolverBaseTest
from .steady_state_solver_test import SteadyStateSolverTest
from .kmc_redistribution_test import KMCRedistributionTest


def suite():
    suite = unittest.TestSuite(
        [unittest.TestLoader().loadTestsFromTestCase(SolverBaseTest),
         unittest.TestLoader().loadTestsFromTestCase(MeanFieldSolverTest),
         unittest.TestLoader().loadTestsFromTestCase(SteadyStateSolverTest),
         unittest.TestLoader().loadTestsFromTestCase(KMCCoveragesPluginTest),
         unittest.TestLoader().loadTestsFromTestCase(KMCTOFPluginTest),
         unittest.TestLoader().loadTestsFromTestCase(KMCFrequencyPluginTest),
         unittest.TestLoader().loadTestsFromTestCase(KMCSolverTest),
         unittest.TestLoader().loadTestsFromTestCase(KMCRedistributionTest),]
    )
    return suite


if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())


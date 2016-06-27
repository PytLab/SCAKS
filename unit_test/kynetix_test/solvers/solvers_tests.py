import unittest

from kmc_solver_test import KMCSolverTest
from kmc_coverages_plugin_test import KMCCoveragesPluginTest
from kmc_frequency_plugin_test import KMCFrequencyPluginTest
from mean_field_solver_test import MeanFieldSolverTest
from solver_base_test import SolverBaseTest
from steady_state_solver_test import SteadyStateSolverTest


def suite():
    suite = unittest.TestSuite(
        [unittest.TestLoader().loadTestsFromTestCase(SolverBaseTest),
         unittest.TestLoader().loadTestsFromTestCase(MeanFieldSolverTest),
         unittest.TestLoader().loadTestsFromTestCase(SteadyStateSolverTest),
         unittest.TestLoader().loadTestsFromTestCase(KMCCoveragesPluginTest),
         unittest.TestLoader().loadTestsFromTestCase(KMCFrequencyPluginTest),
         unittest.TestLoader().loadTestsFromTestCase(KMCSolverTest),]
    )
    return suite


if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())


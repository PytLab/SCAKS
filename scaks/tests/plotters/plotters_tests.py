import unittest

from .energy_profile_plotter_test import EnergyProfilePlotterTest

plotter_test_cases = [EnergyProfilePlotterTest]


def suite():
    suite = unittest.TestSuite([
        unittest.TestLoader().loadTestsFromTestCase(tc) for tc in plotter_test_cases
    ])
    return suite


if __name__ == '__main__':
    result = unittest.TextTestRunner(verbosity=2).run(suite())

    if result.errors or result.failures:
        raise Exception('Get errors and failures')


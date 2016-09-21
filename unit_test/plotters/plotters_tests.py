import unittest

from energy_profile_plotter_test import EnergyProfilePlotterTest


def suite():
    suite = unittest.TestSuite(
         [unittest.TestLoader().loadTestsFromTestCase(EnergyProfilePlotterTest),]
    )
    return suite


if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())


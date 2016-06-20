import unittest

from kynetix.solvers.solver_base import SolverBase


class SolverBaseTest(unittest.TestCase):

    def setUp(self):
        # Test case setting.
        self.maxDiff = None

    def test_get_kTST(self):
        " Test static method get_kTST(). "
        Ga = 1.56
        T = 450

        ref_r = 3.168158762449045e-05
        ret_r = SolverBase.get_kTST(Ga, T)

        self.assertEqual(ref_r, ret_r)

    def test_get_kCT(self):
        " Test statice function get_kCT(). "
        Auc = 9.64e-20
        Ea = 1.56
        p = 1.0
        m = 4.6511864564303997e-26  # absolute mass of CO(kg)
        act_ratio = 0.5
        T = 450
        
        ref_r = 3.82203461455769e-15
        ret_r = SolverBase.get_kCT(Ea, Auc, act_ratio, p, m, T)

        self.assertEqual(ref_r, ret_r)



if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(SolverBaseTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


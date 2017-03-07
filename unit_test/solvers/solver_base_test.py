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
        ret_r = SolverBase.get_kCT(Ea=Ea, Auc=Auc, act_ratio=act_ratio, p=p, m=m, T=T)

        self.assertEqual(ref_r, ret_r)

    def test_get_TST_barrier_from_CT(self):
        " Make sure we can get correct TST barrier from CT rate. "
        kCT = 90235665.7025331
        ref_Ga = 0.5022461602982147
        ret_Ga = SolverBase.get_TST_barrier_from_CT(kCT, 500)

        self.assertEqual(ref_Ga, ret_Ga)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(SolverBaseTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


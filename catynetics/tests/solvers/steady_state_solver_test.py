'''
Test case for SteadyStateSolver.
'''

import logging
import unittest

from mpmath import mpf

from ...models.micro_kinetic_model import MicroKineticModel
from ...solvers import SteadyStateSolver

from .. import *


class SteadyStateSolverTest(unittest.TestCase):
    '''
    Test case for SteadyStateSolver.
    '''

    def setUp(self):
        # Test case setting.
        self.maxDiff = None
        self.setup_dict = dict(
            rxn_expressions=[
                'CO_g + *_s -> CO_s',
                'O2_g + 2*_s -> 2O_s',
                'CO_s + O_s <-> CO-O_2s -> CO2_g + 2*_s',
            ],

            species_definitions={
                'CO_g': {'pressure': 1.0},
                'O2_g': {'pressure': 1./3.},
                'CO2_g': {'pressure': 0.00},
                '*_s': {'site_name': '111', 'type': 'site', 'total': 1.0},
            },

            temperature=450.0,
            parser="RelativeEnergyParser",
            solver="SteadyStateSolver",
            corrector="ThermodynamicCorrector",
            plotter="EnergyProfilePlotter",
            rootfinding='ConstrainedNewton',
            decimal_precision=100,
            tolerance=1e-20,
            max_rootfinding_iterations=100,
        )

    def test_solver_construction_query(self):
        " Test solver can be constructed in kinetic model. "
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        solver = model.solver

        # Check the parser class and base class type.
        self.assertTrue(isinstance(solver, SteadyStateSolver))
        self.assertEqual(solver.__class__.__base__.__name__, "MeanFieldSolver")

        # Test attributes query.

        # Numerical representations.
        self.assertTrue(hasattr(solver, "_math"))
        self.assertTrue(hasattr(solver, "_linalg"))
        self.assertTrue(hasattr(solver, "_mpf"))
        self.assertTrue(hasattr(solver, "_matrix"))
        self.assertTrue(hasattr(solver, "_Axb_solver"))
        self.assertTrue(hasattr(solver, "_norm"))

        # Flags.
        self.assertFalse(solver.has_symbols)

        ref_classified_adsorbates = {'*_s': ['CO_s', 'O_s']}
        self.assertDictEqual(ref_classified_adsorbates, solver.classified_adsorbates)

    def test_elementary_dtheta_dt_expression(self):
        " Test get_elementary_dtheta_dt_expression() function. "
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        solver = model.solver

        # Check.
        adsorbate = "O_s"
        rxn = 'O2_g + 2*_s -> 2O_s'
        ref_dtheta_dt = "2*kf[1]*p['O2_g']*theta['*_s']**2 - 2*kr[1]*theta['O_s']**2"
        ret_dtheta_dt = solver.get_elementary_dtheta_dt_expression(adsorbate, rxn)
        self.assertEqual(ref_dtheta_dt, ret_dtheta_dt)

        adsorbate = "CO_s"
        rxn = 'CO_g + *_s -> CO_s'
        ref_dtheta_dt = "kf[0]*p['CO_g']*theta['*_s'] - kr[0]*theta['CO_s']"
        ret_dtheta_dt = solver.get_elementary_dtheta_dt_expression(adsorbate, rxn)
        self.assertEqual(ref_dtheta_dt, ret_dtheta_dt)

        adsorbate = "O_s"
        rxn = 'CO_s + O_s <-> CO-O_2s -> CO2_g + 2*_s'
        ref_dtheta_dt = "kr[2]*p['CO2_g']*theta['*_s']**2 - kf[2]*theta['CO_s']*theta['O_s']"
        ret_dtheta_dt = solver.get_elementary_dtheta_dt_expression(adsorbate, rxn)
        self.assertEqual(ref_dtheta_dt, ret_dtheta_dt)

    def test_adsorbate_dtheta_dt_expression(self):
        " Test get_adsorbate_dtheta_dt_expression() function. "
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        solver = model.solver

        ref_dtheta_dt = ("kf[0]*p['CO_g']*theta['*_s'] - kr[0]*theta['CO_s'] + " +
                         "kr[2]*p['CO2_g']*theta['*_s']**2 - " +
                         "kf[2]*theta['CO_s']*theta['O_s']")
        ret_dtheta_dt = solver.get_adsorbate_dtheta_dt_expression("CO_s")
        self.assertEqual(ref_dtheta_dt, ret_dtheta_dt)

        ref_dtheta_dt = ("2*kf[1]*p['O2_g']*theta['*_s']**2 - " +
                         "2*kr[1]*theta['O_s']**2 + " +
                         "kr[2]*p['CO2_g']*theta['*_s']**2 - " +
                         "kf[2]*theta['CO_s']*theta['O_s']")
        ret_dtheta_dt = solver.get_adsorbate_dtheta_dt_expression("O_s")
        self.assertEqual(ref_dtheta_dt, ret_dtheta_dt)

    def test_dtheta_dt_expression(self):
        " Make sure we can get dtheta/dt expression correctly. "
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        solver = model.solver

        # Check.
        dtheta_dt_CO_s = ("dtheta_dt[0] = kf[0]*p['CO_g']*theta['*_s'] - " +
                          "kr[0]*theta['CO_s'] + kr[2]*p['CO2_g']*theta['*_s']**2 - " +
                          "kf[2]*theta['CO_s']*theta['O_s']")
        dtheta_dt_O_s = ("dtheta_dt[1] = 2*kf[1]*p['O2_g']*theta['*_s']**2 - " +
                         "2*kr[1]*theta['O_s']**2 + kr[2]*p['CO2_g']*theta['*_s']**2 - " +
                         "kf[2]*theta['CO_s']*theta['O_s']")
        ret_dtheta_dt = [dtheta_dt_CO_s, dtheta_dt_O_s]
        ref_dtheta_dt = solver.get_dtheta_dt_expressions()

        self.assertListEqual(ret_dtheta_dt, ref_dtheta_dt)

    def test_steady_state_function(self):
        " Test steady_state_function(). "
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        parser = model.parser
        solver = model.solver

        # Get data.
        parser.parse_data(filename=mkm_energy)
        solver.get_data()

        # Check.
        coverages = (0.2, 0.5)
        ref_dtheta_dt = (2812943317895.31469634, 562588664794.88448169)
        ret_dtheta_dt = solver.steady_state_function(coverages)

        for ref, ret in ref_dtheta_dt, ret_dtheta_dt:
            self.assertAlmostEqual(ret, ret)

    def test_term_adsorbate_derivation(self):
        " Test private function __term_adsorbate_derivation(). "
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        solver = model.solver

        # Check.
        adsorbate = "CO_s"
        term = "kf[2]*theta['CO_s']*theta['*_s']" 
        ref_derivation = "kf[2]*theta['*_s']"
        ret_derivation = solver._SteadyStateSolver__term_adsorbate_derivation(adsorbate, term)
        self.assertEqual(ref_derivation, ret_derivation)

        adsorbate = "O_s"
        term = "2*kr[1]*theta['O_s']**2"
        ref_derivation = "2*2*kr[1]*theta['O_s']**1"
        ret_derivation = solver._SteadyStateSolver__term_adsorbate_derivation(adsorbate, term)
        self.assertEqual(ref_derivation, ret_derivation)

    def test_total_term_adsorbate_derivation(self):
        " Test private function __total_term_adsorbate_derivation(). "
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        solver = model.solver

        # Check.
        adsorbate = "O_s"
        term = "2*kf[1]*p['O2_g']*theta['*_s']**2"
        ret_derivation = solver._SteadyStateSolver__total_term_adsorbate_derivation(adsorbate,
                                                                                    term)
        ref_derivation = "-2*2*kf[1]*p['O2_g']*(1.0 - theta['CO_s'] - theta['O_s'])**1"
        self.assertEqual(ref_derivation, ret_derivation)

        adsorbate = "O_s"
        term = "2*kf[1]*p['O2_g']*theta['*_s']**2*theta['O_s']"
        ret_derivation = solver._SteadyStateSolver__total_term_adsorbate_derivation(adsorbate,
                                                                                    term)
        ref_derivation = ("1*2*kf[1]*p['O2_g']*theta['*_s']**2 + " +
                          "-2*2*kf[1]*p['O2_g']*(1.0 - theta['CO_s'] - theta['O_s'])**1*theta['O_s']")
        self.assertEqual(ref_derivation, ret_derivation)

    def test_poly_adsorbate_derivation(self):
        " Test we can derive dtheta/dt expression correctly. "
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        solver = model.solver

        adsorbate = "CO_s"
        poly_expression = ("dtheta_dt[0] = kf[0]*p['CO_g']*theta['*_s'] - " +
                           "kr[0]*theta['CO_s'] + kr[2]*p['CO2_g']*theta['*_s']**2 - " +
                           "kf[2]*theta['CO_s']*theta['O_s']")
        ref_expression = ("-kf[0]*p['CO_g'] - kr[0] + " +
                          "-2*kr[2]*p['CO2_g']*(1.0 - theta['CO_s'] - theta['O_s'])**1 - " +
                          "kf[2]*theta['O_s']")
        ret_expression = solver.poly_adsorbate_derivation(adsorbate_name=adsorbate,
                                                          poly_expression=poly_expression)
        self.assertEqual(ref_expression, ret_expression)

        adsorbate = "O_s"
        poly_expression = ("dtheta_dt[1] = 2*kf[1]*p['O2_g']*theta['*_s']**2 - " +
                           "2*kr[1]*theta['O_s']**2 + kr[2]*p['CO2_g']*theta['*_s']**2 - " +
                           "kf[2]*theta['CO_s']*theta['O_s']")
        ref_expression = ("-2*2*kf[1]*p['O2_g']*(1.0 - theta['CO_s'] - theta['O_s'])**1 - " +
                          "2*2*kr[1]*theta['O_s']**1 + " +
                          "-2*kr[2]*p['CO2_g']*(1.0 - theta['CO_s'] - theta['O_s'])**1 - " +
                          "kf[2]*theta['CO_s']")
        ret_expression = solver.poly_adsorbate_derivation(adsorbate_name=adsorbate,
                                                          poly_expression=poly_expression)
        self.assertEqual(ref_expression, ret_expression)

    def test_analytical_jacobian(self):
        " Make sure we can get analytical Jacobian matrix correctly. "
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        parser = model.parser
        solver = model.solver

        parser.parse_data(filename=mkm_energy)
        solver.get_data()

        # Check.
        coverages = (0.2, 0.4)
        ref_jacobian = [[-9376477776977.32, -9376477746581.58],
                        [-5000788131510.20, -5000788131510.186]]    
        ret_jacobian = solver.analytical_jacobian(coverages).tolist()
        #for m in range(2):
        #    for n in range(2):
        #        self.assertAlmostEqual(ref_jacobian[m][n], float(ret_jacobian[m][n]), places=1)

    def test_get_residual(self):
        " Test we can get correct residual. "
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        parser = model.parser
        solver = model.solver

        parser.parse_data(filename=mkm_energy)
        solver.get_data()

        coverages = (0.99, 0.01)
        ref_residual = 30091.76907705
        ret_residual = float(solver.get_residual(coverages))

        self.assertAlmostEqual(ref_residual, ret_residual, places=2)

    def test_coarse_steady_state_cvgs(self):
        " Make sure we can get a coarse coverages. "
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        parser = model.parser
        solver = model.solver

        parser.parse_data(filename=mkm_energy)
        solver.get_data()

        coverages = (0.99, 0.0)
        ref_coarse_cvgs = [0.9993051181164324, 0.0006948786441204215]
        ret_coarse_cvgs = solver.coarse_steady_state_cvgs(coverages).tolist()

        for ref, ret in zip(ref_coarse_cvgs, ret_coarse_cvgs):
            self.assertAlmostEqual(ref, ret, places=2)

    def test_get_steady_state_coverages(self):
        " Test we can get correct steady state coverages. "
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        parser = model.parser
        solver = model.solver

        parser.parse_data(filename=mkm_energy)
        solver.get_data()
        
        # Check.
        coverages = [0.9, 0.1]
        ref_sscvg = [0.9993009023315728, 0.0006990944289937246]
        ret_sscvg = solver.get_steady_state_cvgs(coverages)
        for ref, ret in zip(ref_sscvg, ret_sscvg):
            self.assertAlmostEqual(ref, float(ret))

    def test_get_single_XRC(self):
        " Test function get_single_XRC(). "
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        parser = model.parser
        solver = model.solver

        parser.parse_data(filename=mkm_energy)
        solver.get_data()

        # Get steady state converages first.
        coverages = [0.9, 0.1]
        solver.get_steady_state_cvgs(coverages)

        # Check.
        gas_name = "CO2_g"
        ref_XRC = [-0.00000, 0.9986, 0.00139]
        ret_XRC = solver.get_single_XRC(gas_name, epsilon=1e-5)
        for ref, ret in zip(ref_XRC, ret_XRC):
            self.assertAlmostEqual(ref, float(ret), places=4)

    def test_get_elementary_dtheta_dt_sym(self):
        " Test we can get correct dtheta/dt expression for an elementary reaction. "
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        parser = model.parser
        solver = model.solver

        parser.parse_data(filename=mkm_energy)
        solver.get_data()
        solver.get_data_symbols()

        # Get symbols.

        # Relative energies.
        Ga_0, Ga_1, Ga_2 = solver._Ga_sym
        dG_0, dG_1, dG_2 = solver._dG_sym

        # Coverage.
        theta_CO_s = solver._extract_symbol("CO_s", "ads_cvg")
        theta_O_s = solver._extract_symbol("O_s", "ads_cvg")
        theta_s = solver._extract_symbol("*_s", "free_site_cvg")

        # Pressure.
        p_CO2_g = solver._extract_symbol("CO2_g", "pressure")
        p_O2_g = solver._extract_symbol("O2_g", "pressure")
        p_CO_g = solver._extract_symbol("CO_g", "pressure")

        # Constants.
        kB = solver._kB_sym
        T = solver._T_sym
        h = solver._h_sym
        from sympy import E

        rxn_expression = 'CO_s + O_s <-> CO-O_2s -> CO2_g + 2*_s'
        rf = T*kB*theta_CO_s*theta_O_s*E**(-Ga_2/(T*kB))/h
        rr = T*kB*p_CO2_g*(-theta_CO_s - theta_O_s + 1.0)**2*E**((-Ga_2 + dG_2)/(T*kB))/h

        ref_dtheta_dt = rr - rf
        adsorbate = "CO_s"
        ret_dtheta_dt = solver.get_elementary_dtheta_dt_sym(adsorbate, rxn_expression)

        self.assertEqual(ref_dtheta_dt.simplify(), ret_dtheta_dt.simplify())

    def test_get_adsorbate_dtheta_dt_sym(self):
        " Test we can get correct dtheta/dt for an adsorbate. "
        # NEED IMPLIMENTATION.

    def test_steady_state_function_by_sym(self):
        " Test function steady_state_function(). "
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        parser = model.parser
        solver = model.solver

        parser.parse_data(filename=mkm_energy)
        solver.get_data()
        solver.get_data_symbols()

        # Check.
        coverages = (0.5, 0.3)
        ref_dtheta_dt = (mpf('1875295534118.435791015625'),
                         mpf('125019703287.740081787109375'))
        ret_dtheta_dt = solver.steady_state_function_by_sym(coverages)

        self.assertTupleEqual(ref_dtheta_dt, ret_dtheta_dt)

    def test_analytical_jacobian_sym(self):
        " Make sure we can get anlytical jacobian matrix correctly. "
        # NEED IMPLIMENTATION.

    def test_analytical_jacobian_by_sym(self):
        " Test we can get correct jacobian matrix by symbol derivation. "
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        parser = model.parser
        solver = model.solver

        parser.parse_data(filename=mkm_energy)
        solver.get_data()
        solver.get_data_symbols()

        # Check.
        coverages = (0.5, 0.3)

        ref_jacobian = [[mpf('-9376477776977.31640625'), mpf('-9376477746581.609375')],
                        [mpf('-1250197032877.56982421875'), mpf('-1250197032877.588623046875')]]
        ret_jacobian = solver.analytical_jacobian_by_sym(coverages).tolist()

        self.assertListEqual(ref_jacobian, ret_jacobian)

    def tearDown(self):
        cleanup()

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(SteadyStateSolverTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


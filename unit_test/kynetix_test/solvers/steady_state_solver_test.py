import logging
import re
import unittest

import numpy as np
from mpmath import mpf

from kynetix.model import KineticModel
from kynetix.solvers import *


class SteadyStateSolverTest(unittest.TestCase):

    def setUp(self):
        # Test case setting.
        self.maxDiff = None

    def test_solver_construction_query(self):
        " Test solver can be constructed in kinetic model. "
        # Construction.
        model = KineticModel(setup_file="input_files/steady_state_solver.mkm",
                             verbosity=logging.WARNING)
        solver = model.solver()

        # Check the parser class and base class type.
        self.assertTrue(isinstance(solver, SteadyStateSolver))
        self.assertEqual(solver.__class__.__base__.__name__, "SolverBase")

        # Test attributes query.

        # Default protected attributes.
        self.assertTrue(hasattr(solver, "_perturbation_size"))
        self.assertTrue(hasattr(solver, "_perturbation_direction"))
        self.assertTrue(hasattr(solver, "_numerical_representation"))
        self.assertTrue(hasattr(solver, "_archived_variables"))

        # Numerical representations.
        self.assertTrue(hasattr(solver, "_math"))
        self.assertTrue(hasattr(solver, "_linalg"))
        self.assertTrue(hasattr(solver, "_mpf"))
        self.assertTrue(hasattr(solver, "_matrix"))
        self.assertTrue(hasattr(solver, "_Axb_solver"))
        self.assertTrue(hasattr(solver, "_norm"))

        # Flags.
        self.assertFalse(solver.has_absolute_energy())
        self.assertFalse(solver.has_relative_energy())
        self.assertFalse(solver.has_energy_correction())
        self.assertFalse(solver.has_symbols())

        ref_classified_adsorbates = {'s': ['CO_s', 'O_s']}
        self.assertDictEqual(ref_classified_adsorbates, solver.classified_adsorbates())

    def test_get_data(self):
        " Test solver can get data correctly. "
        # Construction.
        model = KineticModel(setup_file="input_files/steady_state_solver.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        solver = model.solver()

        # Get data before parsing data, an exception would be expected.
        self.assertRaises(IOError, solver.get_data)

        # Parse data.
        parser.parse_data(relative=True, filename="input_files/rel_energy.py")
        solver.get_data()

        self.assertTrue(solver.has_relative_energy())

        # Check pressure.
        ref_pressures = {'CO2_g': mpf('0.0'), 'CO_g': mpf('1.0'), 'O2_g': mpf(1./3.)}
        self.assertDictEqual(ref_pressures, solver.pressures())

        # Check concentrations.
        ref_concentrations = {}
        self.assertDictEqual(ref_concentrations, solver.concentrations())

        # Relative energies.
        ref_relative_energies = {'Gaf': [0.0, 0.0, 1.25],
                                 'Gar': [0.758, 2.64, 0.9259999999999999],
                                 'dG': [-0.758, -2.64, 0.324]}
        self.assertDictEqual(ref_relative_energies, solver.relative_energies())

        # Formation energies.
        self.assertRaises(AttributeError, solver.formation_energies)

        # Parse absolute energies.
        parser.parse_data(relative=False, filename="input_files/rel_energy.py")
        solver.get_data()

        self.assertTrue(solver.has_absolute_energy())

        # Check formation energies.
        ref_formation_energies = {'*_s': mpf('0.0'),
                                  'CO-O_2s': mpf('0.9259999999999999342747969421907328069210052490234375'),
                                  'CO2_g': mpf('0.0'),
                                  'CO_g': mpf('0.0'),
                                  'CO_s': mpf('-0.75800000000000000710542735760100185871124267578125'),
                                  'O2_g': mpf('3.50800000000000000710542735760100185871124267578125'),
                                  'O_s': mpf('0.4339999999999999413802242997917346656322479248046875')}
        self.assertDictEqual(ref_formation_energies, solver.formation_energies())

    def test_elementary_dtheta_dt_expression(self):
        " Test get_elementary_dtheta_dt_expression() function. "
        # Construction.
        model = KineticModel(setup_file="input_files/steady_state_solver.mkm",
                             verbosity=logging.WARNING)
        solver = model.solver()

        # Check.
        adsorbate = "O_s"
        rxn = [['O2_g', '2*_s'], ['2O_s']]
        ref_dtheta_dt = "2*kf[1]*p['O2_g']*theta['*_s']**2 - 2*kr[1]*theta['O_s']**2"
        ret_dtheta_dt = solver.get_elementary_dtheta_dt_expression(adsorbate, rxn)
        self.assertEqual(ref_dtheta_dt, ret_dtheta_dt)

        adsorbate = "CO_s"
        rxn = [['CO_g', '*_s'], ['CO_s']]
        ref_dtheta_dt = "kf[0]*p['CO_g']*theta['*_s'] - kr[0]*theta['CO_s']"
        ret_dtheta_dt = solver.get_elementary_dtheta_dt_expression(adsorbate, rxn)
        self.assertEqual(ref_dtheta_dt, ret_dtheta_dt)

        adsorbate = "O_s"
        rxn = [['CO_s', 'O_s'], ['CO-O_2s'], ['CO2_g', '2*_s']]
        ref_dtheta_dt = "kr[2]*p['CO2_g']*theta['*_s']**2 - kf[2]*theta['CO_s']*theta['O_s']"
        ret_dtheta_dt = solver.get_elementary_dtheta_dt_expression(adsorbate, rxn)
        self.assertEqual(ref_dtheta_dt, ret_dtheta_dt)

    def test_adsorbate_dtheta_dt_expression(self):
        " Test get_adsorbate_dtheta_dt_expression() function. "
        # Construction.
        model = KineticModel(setup_file="input_files/steady_state_solver.mkm",
                             verbosity=logging.WARNING)
        solver = model.solver()

        ref_dtheta_dt = ("kf[0]*p['CO_g']*theta['*_s'] - kr[0]*theta['CO_s'] + " +
                         "kr[2]*p['CO2_g']*theta['*_s']**2 - " +
                         "kf[2]*theta['CO_s']*theta['O_s']")
        ret_dtheta_dt = solver.get_adsorbate_dtheta_dt_expression("CO_s")

        ref_dtheta_dt = ("2*kf[1]*p['O2_g']*theta['*_s']**2 - " +
                         "2*kr[1]*theta['O_s']**2 + " +
                         "kr[2]*p['CO2_g']*theta['*_s']**2 - " +
                         "kf[2]*theta['CO_s']*theta['O_s']")
        ret_dtheta_dt = solver.get_adsorbate_dtheta_dt_expression("O_s")

    def test_dtheta_dt_expression(self):
        " Make sure we can get dtheta/dt expression correctly. "
        # Construction.
        model = KineticModel(setup_file="input_files/steady_state_solver.mkm",
                             verbosity=logging.WARNING)
        solver = model.solver()

        # Check.
        dtheta_dt_CO_s = ("dtheta_dt[0] = kf[0]*p['CO_g']*theta['*_s'] - " +
                          "kr[0]*theta['CO_s'] + kr[2]*p['CO2_g']*theta['*_s']**2 - " +
                          "kf[2]*theta['CO_s']*theta['O_s']")
        dtheta_dt_O_s = ("dtheta_dt[1] = 2*kf[1]*p['O2_g']*theta['*_s']**2 - " +
                         "2*kr[1]*theta['O_s']**2 + kr[2]*p['CO2_g']*theta['*_s']**2 - " +
                         "kf[2]*theta['CO_s']*theta['O_s']")
        ret_dtheta_dt = (dtheta_dt_CO_s, dtheta_dt_O_s)
        ref_dtheta_dt = solver.get_dtheta_dt_expressions()

        self.assertTupleEqual(ret_dtheta_dt, ref_dtheta_dt)

    def test_steady_state_function(self):
        " Test steady_state_function(). "
        # Construction.
        model = KineticModel(setup_file="input_files/steady_state_solver.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        solver = model.solver()

        # Get data.
        parser.parse_data(filename="input_files/rel_energy.py")
        solver.get_data()

        # Check.
        coverages = (0.2, 0.5)
        ref_dtheta_dt = (mpf('2812943317895.314716929552474856840195790042878389374664419807065157914105124045016669624216006091895869'),
                         mpf('562588664794.8844858075949325690317587695978287990072657061366848036962738449361619993106651762296008292'))
        ret_dtheta_dt = solver.steady_state_function(coverages)

        self.assertTupleEqual(ref_dtheta_dt, ret_dtheta_dt)

    def test_term_adsorbate_derivation(self):
        " Test private function __term_adsorbate_derivation(). "
        # Construction.
        model = KineticModel(setup_file="input_files/steady_state_solver.mkm",
                             verbosity=logging.WARNING)
        solver = model.solver()

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

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(SteadyStateSolverTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


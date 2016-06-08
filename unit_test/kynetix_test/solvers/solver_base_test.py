import logging
import re
import unittest

import numpy as np
from mpmath import mpf

from kynetix.model import KineticModel
from kynetix.solvers import *


class SolverBaseTest(unittest.TestCase):

    def setUp(self):
        # Test case setting.
        self.maxDiff = None

    def test_solver_construction_query(self):
        " Test solver can be constructed in kinetic model. "
        # Construction.
        model = KineticModel(setup_file="input_files/solver_base.mkm",
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
        model = KineticModel(setup_file="input_files/solver_base.mkm",
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
                                  'CO-O_2s': mpf('0.9259999999995'),
                                  'CO2_g': mpf('0.0'),
                                  'CO_g': mpf('0.0'),
                                  'CO_s': mpf('-0.7580000000016'),
                                  'O2_g': mpf('3.508000000002'),
                                  'O_s': mpf('0.4340000000011')}
        self.assertDictEqual(ref_formation_energies, solver.formation_energies())

    def test_get_rate_constants(self):
        " Make sure we can get rate constants correctly. "
        # Construction.
        model = KineticModel(setup_file="input_files/solver_base.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        solver = model.solver()

        self.assertRaises(AttributeError, solver.get_rate_constants)

        parser.parse_data(filename="input_files/rel_energy.py")
        solver.get_data()

        # Check rate constants.
        ref_forward_rate_constants = (mpf('9376477746560.0'),
                                      mpf('9376477746560.0'),
                                      mpf('0.09389759708665'))
        ref_reverse_rate_constants = (mpf('30395.72540021'),
                                      mpf('2.542951527174e-17'),
                                      mpf('399.2961611934'))
        ret_forward_rate_constants, ret_reverse_rate_constants = solver.get_rate_constants()
        self.assertTupleEqual(ref_forward_rate_constants, ret_forward_rate_constants)
        self.assertTupleEqual(ref_reverse_rate_constants, ret_reverse_rate_constants)

    def test_boltzmann_coverages(self):
        " Test we can get the Boltzmann converages. "
        # Construction.
        model = KineticModel(setup_file="input_files/solver_base.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        solver = model.solver()

        # Check Boltzmann before parsing absolute energies.
        parser.parse_data(filename="input_files/rel_energy.py", relative=True)
        solver.get_data()

        # An exception would be expected.
        self.assertRaises(IOError, solver.boltzmann_coverages)

        # Parse absolute energies.
        parser.parse_data(filename="input_files/rel_energy.py", relative=False)
        solver.get_data()

        # Check boltzmann coverages.
        ref_coverages = (mpf('0.9999999967549'), mpf('4.468751710442e-14'))
        ret_coverages = solver.boltzmann_coverages()

        self.assertTupleEqual(ref_coverages, ret_coverages)

        # Without empty sites.
        ref_coverages = (mpf('1.0'), mpf('4.468751724917e-14'))
        ret_coverages = solver.boltzmann_coverages(include_empty_site=False)

        self.assertTupleEqual(ref_coverages, ret_coverages)

    def test_elementary_rate_expression(self):
        "Make sure we can get the rate expression for an elementary reaction correctly. "
        # Construction.
        model = KineticModel(setup_file="input_files/solver_base.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        solver = model.solver()

        rxn_expression = 'CO_g + *_s -> CO_s'
        ref_f_expr = "kf[0]*p['CO_g']*theta['*_s']"
        ref_r_expr = "kr[0]*theta['CO_s']"
        ret_f_expr, ret_r_expr = solver.get_elementary_rate_expression(rxn_expression)

        self.assertEqual(ref_f_expr, ret_f_expr)
        self.assertEqual(ref_r_expr, ret_r_expr)

        # Check elementary reaction with TS.
        rxn_expression = 'CO_s + O_s <-> CO-O_2s -> CO2_g + 2*_s'
        ref_f_expr = "kf[2]*theta['CO_s']*theta['O_s']"
        ref_r_expr = "kr[2]*p['CO2_g']*theta['*_s']**2"
        ret_f_expr, ret_r_expr = solver.get_elementary_rate_expression(rxn_expression)

        self.assertEqual(ref_f_expr, ret_f_expr)
        self.assertEqual(ref_r_expr, ret_r_expr)

    def test_rate_expressions(self):
        " Test we can get all rate expressions correctly. "
        # Construction.
        model = KineticModel(setup_file="input_files/solver_base.mkm",
                             verbosity=logging.WARNING)
        solver = model.solver()

        ref_rate_expressions = (["rfs[0] = kf[0]*p['CO_g']*theta['*_s']",
                                 "rfs[1] = kf[1]*p['O2_g']*theta['*_s']**2",
                                 "rfs[2] = kf[2]*theta['CO_s']*theta['O_s']"],
                                ["rrs[0] = kr[0]*theta['CO_s']",
                                 "rrs[1] = kr[1]*theta['O_s']**2",
                                 "rrs[2] = kr[2]*p['CO2_g']*theta['*_s']**2"])
        ret_rate_expressions = solver.get_rate_expressions()

        self.assertTupleEqual(ref_rate_expressions, ret_rate_expressions)

    def test_get_rates(self):
        " Make sure we can get rates correctly. "
        # Construction.
        model = KineticModel(setup_file="input_files/solver_base.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        solver = model.solver()

        # Get data.
        parser.parse_data(filename="input_files/rel_energy.py")
        solver.get_data()

        # Check.
        coverages = (0.5, 0.5)
        ret_rates = solver.get_rates(coverages)
        ref_rates = ((mpf('0.0'), mpf('0.0'), mpf('0.02347439927166')),
                     (mpf('15197.8627001'), mpf('6.357378817934e-18'), mpf('0.0')))
        self.assertTupleEqual(ref_rates, ret_rates)

        # Check net rates.
        ref_net_rates = (mpf('-15197.8627001'),
                         mpf('-6.357378817934e-18'),
                         mpf('0.02347439927166'))
        ret_net_rates = solver.get_net_rates(coverages)
        self.assertTupleEqual(ref_net_rates, ret_net_rates)

    def test_get_reversibilities(self):
        " Make sure we can get the correct reversibilities. "
        # Construction.
        model = KineticModel(setup_file="input_files/solver_base.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        solver = model.solver()

        parser.parse_data(filename="input_files/rel_energy.py")
        solver.get_data()
        rfs, rrs = solver.get_rates((0.2, 0.5))

        ref_reversibilities = [2.1611331548106183e-09, 2.260045114947677e-29, 0.0]
        ret_reversibilities = solver.get_reversibilities(rfs, rrs)

        self.assertListEqual(ref_reversibilities, ret_reversibilities)

    def test_get_tof(self):
        " Test we can get TOF correctly. "
        # Construction.
        model = KineticModel(setup_file="input_files/solver_base.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        solver = model.solver()

        parser.parse_data(filename="input_files/rel_energy.py")
        solver.get_data()

        coverages = (0.2, 0.4)
        ref_tof = [mpf('0.007511807766946'),
                   mpf('-3750591092544.0'),
                   mpf('-500078813148.0')]
        ret_tof = solver.get_tof(coverages)

        self.assertListEqual(ref_tof, ret_tof)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(SolverBaseTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


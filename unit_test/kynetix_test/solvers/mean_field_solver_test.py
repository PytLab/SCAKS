import logging
import re
import unittest

import numpy as np
from mpmath import mpf

from kynetix.model import KineticModel
from kynetix.parsers.rxn_parser import *
from kynetix.solvers import *


class MeanFieldSolverTest(unittest.TestCase):

    def setUp(self):
        # Test case setting.
        self.maxDiff = None

    def test_solver_construction_query(self):
        # {{{
        " Test solver can be constructed in kinetic model. "
        # Construction.
        model = KineticModel(setup_file="input_files/solver_base.mkm",
                             verbosity=logging.WARNING)
        solver = model.solver()

        # Check the parser class and base class type.
        self.assertTrue(isinstance(solver, SteadyStateSolver))
        self.assertEqual(solver.__class__.__base__.__name__, "MeanFieldSolver")

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
        # }}}

    def test_get_data(self):
        # {{{
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
        # }}}

    def test_get_state_energy(self):
        " Test we can get correct state energy. "
        # Construction.
        model = KineticModel(setup_file="input_files/solver_base.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        parser.parse_data(filename="input_files/rel_energy.py")
        solver = model.solver()
        solver.get_data()

        # Check.
        state = ChemState('CO_s + O_s')
        ref_G = solver._G['CO_s'] + solver._G['O_s']
        ret_G = solver._get_state_energy(state)

        self.assertEqual(ref_G, ret_G)

    def test_get_single_relative_energies(self):
        " Make sure we can get correct relative energy for an elementary reaction. "
        # Construction.
        model = KineticModel(setup_file="input_files/solver_base.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        parser.parse_data(filename="input_files/rel_energy.py")
        solver = model.solver()
        solver.get_data()

        # Check.
        rxn_expression = 'CO_g + *_s -> CO_s'
        ref_e = (mpf('0.0'), mpf('0.7580000000016'), mpf('-0.7580000000016'))
        ret_e = solver.get_single_relative_energies(rxn_expression)
        self.assertTupleEqual(ref_e, ret_e)

        rxn_expression = 'O2_g + 2*_s -> 2O_s'
        ref_e = (mpf('0.0'), mpf('2.640000000014'), mpf('-2.640000000014'))
        ret_e = solver.get_single_relative_energies(rxn_expression)
        self.assertTupleEqual(ref_e, ret_e)

        rxn_expression = 'CO_s + O_s <-> CO-O_2s -> CO2_g + 2*_s'
        ref_e = (mpf('1.25'), mpf('0.9259999999995'), mpf('0.3240000000005'))
        ret_e = solver.get_single_relative_energies(rxn_expression)
        self.assertTupleEqual(ref_e, ret_e)

    def test_get_relative_from_absolute(self):
        " Test we can get relative energies from absolute energies correctly. "
        # Construction.
        model = KineticModel(setup_file="input_files/solver_base.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        parser.parse_data(filename="input_files/rel_energy.py")
        solver = model.solver()
        solver.get_data()

        # Check.
        ref_e = {'Gaf': [mpf('0.0'), mpf('0.0'), mpf('1.25')],
                 'Gar': [mpf('0.7580000000016'),
                         mpf('2.640000000014'),
                         mpf('0.9259999999995')],
                 'dG': [mpf('-0.7580000000016'),
                        mpf('-2.640000000014'),
                        mpf('0.3240000000005')]}
        ret_e = solver.get_relative_from_absolute()

        self.assertDictEqual(ref_e, ret_e)

    def test_get_rate_constants(self):
        # {{{
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
                                      mpf('0.09389759708756'))
        ref_reverse_rate_constants = (mpf('30395.7254014'),
                                      mpf('2.542951526972e-17'),
                                      mpf('399.296161212'))
        ret_forward_rate_constants, ret_reverse_rate_constants = solver.get_rate_constants()
        self.assertTupleEqual(ref_forward_rate_constants, ret_forward_rate_constants)
        self.assertTupleEqual(ref_reverse_rate_constants, ret_reverse_rate_constants)
        # }}}

    def test_boltzmann_coverages(self):
        # {{{
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
        # }}}

    def test_elementary_rate_expression(self):
        # {{{
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
        # }}}

    def test_rate_expressions(self):
        # {{{
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
        # }}}

    def test_get_rates(self):
        # {{{
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
        ref_rates = ((mpf('0.0'), mpf('0.0'), mpf('0.02347439927189')),
                     (mpf('15197.8627007'), mpf('6.357378817428e-18'), mpf('0.0')))
        self.assertTupleEqual(ref_rates, ret_rates)

        # Check net rates.
        ref_net_rates = (mpf('-15197.8627007'),
                         mpf('-6.357378817429e-18'),
                         mpf('0.02347439927189'))
        ret_net_rates = solver.get_net_rates(coverages)
        self.assertTupleEqual(ref_net_rates, ret_net_rates)
        # }}}

    def test_get_reversibilities(self):
        # {{{
        " Make sure we can get the correct reversibilities. "
        # Construction.
        model = KineticModel(setup_file="input_files/solver_base.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        solver = model.solver()

        parser.parse_data(filename="input_files/rel_energy.py")
        solver.get_data()
        rfs, rrs = solver.get_rates((0.2, 0.5))

        ref_reversibilities = [2.1611331548919335e-09, 2.260045114764006e-29, 0.0]
        ret_reversibilities = solver.get_reversibilities(rfs, rrs)

        self.assertListEqual(ref_reversibilities, ret_reversibilities)
        # }}}

    def test_get_tof(self):
        # {{{
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

        # Test gas tof.
        ref_tof = mpf('-3750591092544.0')
        ret_tof = solver.get_tof(coverages, gas_name="CO_g")
        self.assertEqual(ref_tof, ret_tof)
        # }}}

    # ----------------------------------------------------------------
    # Symbol tests.

    def test_get_data_symbols(self):
        # {{{
        " Make sure we can get all correct symbols. "
        # Construction.
        model = KineticModel(setup_file="input_files/solver_base.mkm",
                             verbosity=logging.WARNING)
        solver = model.solver()

        solver.get_data_symbols()

        # Check P symbols.
        ref_p_symbols = ('p_CO2_g', 'p_CO_g', 'p_O2_g')
        ret_p_symbols = solver._p_sym

        for p, p_str in zip(ret_p_symbols, ref_p_symbols):
            self.assertEqual(p.name, p_str)

        # Check concentration symbols.
        ref_c_symbols = ()
        ret_c_symbols = solver._c_sym

        for p, p_str in zip(ret_p_symbols, ref_p_symbols):
            self.assertEqual(p.name, p_str)

        # Check adsorbate coverage symbols.
        ref_ads_theta_symbols = ()
        ret_ads_theta_symbols = solver._ads_theta_sym

        for ads_theta, ads_theta_str in zip(ret_ads_theta_symbols, ref_ads_theta_symbols):
            self.assertEqual(ads_theta.name, ads_theta_str)

        # Check free site coverage symbols.
        ref_fsite_theta_symbols = ()
        ret_fsite_theta_symbols = solver._fsite_theta_sym

        for fsite_theta, fsite_theta_str in zip(ret_fsite_theta_symbols, ref_fsite_theta_symbols):
            self.assertEqual(fsite_theta.name, fsite_theta_str)

        # Check free energy coverage symbols.
        ref_G_symbols = ()
        ret_G_symbols = solver._G_sym

        for G, G_str in zip(ret_G_symbols, ref_G_symbols):
            self.assertEqual(G.name, G_str)

        # Check K coverage symbols.
        ref_K_symbols = ()
        ret_K_symbols = solver._K_sym

        for K, K_str in zip(ret_K_symbols, ref_K_symbols):
            self.assertEqual(K.name, K_str)
        # }}}

    def test_extract_symbol(self):
        # {{{
        " Test protected function _extract_symbol(). "
        # Construction.
        model = KineticModel(setup_file="input_files/solver_base.mkm",
                             verbosity=logging.WARNING)
        solver = model.solver()

        solver.get_data_symbols()

        # Pressure symbol.
        ref_pressure = 'p_CO2_g'
        ret_pressure = solver._extract_symbol('CO2_g', 'pressure')
        self.assertEqual(ret_pressure.name, ref_pressure)

        # Adsorbate symbol.
        ref_cvg = 'theta_CO_s'
        ret_cvg = solver._extract_symbol('CO_s', 'ads_cvg')
        self.assertEqual(ret_cvg.name, ref_cvg)

        # Empty site symbol.
        ret_cvg = solver._extract_symbol('s', 'free_site_cvg')
        CO = solver._extract_symbol('CO_s', 'ads_cvg')
        O = solver._extract_symbol('O_s', 'ads_cvg')
        ref_cvg = 1.0 - CO - O
        self.assertEqual(ret_cvg, ref_cvg)

        # Free energy symbol.
        ref_G = 'G_CO_g'
        ret_G = solver._extract_symbol('CO_g', 'free_energy')
        self.assertEqual(ret_G.name, ref_G)
        # }}}

    def test_get_single_barrier_symbols(self):
        # {{{
        " Make sure we can get correct barrier expression for an elementary reaction. "
        # Construction.
        model = KineticModel(setup_file="input_files/solver_base.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        solver = model.solver()

        parser.parse_data(filename="input_files/rel_energy.py")
        solver.get_data()
        solver.get_data_symbols()

        rxn_expression = 'CO_s + O_s <-> CO-O_2s -> CO2_g + 2*_s'
        ret_barrier_symbols = solver.get_single_barrier_symbols(rxn_expression)
        COO_2s = solver._extract_symbol("CO-O_2s","free_energy")
        CO_s = solver._extract_symbol("CO_s","free_energy")
        O_s = solver._extract_symbol("O_s","free_energy")
        CO2_g = solver._extract_symbol("CO2_g","free_energy")
        s = solver._extract_symbol("s","free_energy")
        ref_barrier_symbols = (COO_2s - CO_s - O_s, -2*s + COO_2s - CO2_g)

        self.assertTupleEqual(ref_barrier_symbols, ret_barrier_symbols)
        # }}}

    def test_get_barrier_symbols(self):
        # {{{
        " Make sure we can get all barrier expressions correctly. "
        # Construction.
        model = KineticModel(setup_file="input_files/solver_base.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        solver = model.solver()

        parser.parse_data(filename="input_files/rel_energy.py")
        solver.get_data()
        solver.get_data_symbols()

        # Check.
        ret_Gaf_symbols, ret_Gar_symbols = solver.get_barrier_symbols()

        # Get references.
        COO_2s = solver._extract_symbol("CO-O_2s", "free_energy")
        CO_s = solver._extract_symbol("CO_s", "free_energy")
        O_s = solver._extract_symbol("O_s", "free_energy")
        CO2_g = solver._extract_symbol("CO2_g", "free_energy")
        O2_g = solver._extract_symbol("O2_g", "free_energy")
        CO_g = solver._extract_symbol("CO_g", "free_energy")
        s = solver._extract_symbol("s", "free_energy")

        ref_Gaf_symbols = [0, 0, COO_2s - CO_s - O_s]
        ref_Gar_symbols = [CO_g + s - CO_s,
                           O2_g + 2*s - 2*O_s,
                           COO_2s - CO2_g - 2*s]

        self.assertListEqual(ret_Gaf_symbols, ref_Gaf_symbols)
        self.assertListEqual(ret_Gar_symbols, ref_Gar_symbols)
        # }}}

    def test_get_rate_constant_symbols(self):
        # {{{
        " Test we can get get correct rate constants symbols. "
        # Construction.
        model = KineticModel(setup_file="input_files/solver_base.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        solver = model.solver()

        parser.parse_data(filename="input_files/rel_energy.py")
        solver.get_data()
        solver.get_data_symbols()

        # Get symbols.
        COO_2s = solver._extract_symbol("CO-O_2s", "free_energy")
        CO_s = solver._extract_symbol("CO_s", "free_energy")
        O_s = solver._extract_symbol("O_s", "free_energy")
        CO2_g = solver._extract_symbol("CO2_g", "free_energy")
        O2_g = solver._extract_symbol("O2_g", "free_energy")
        CO_g = solver._extract_symbol("CO_g", "free_energy")
        s = solver._extract_symbol("s", "free_energy")
        kB = solver._kB_sym
        T = solver._T_sym
        h = solver._h_sym
        from sympy import E

        ref_kf_syms = [T*kB/h, T*kB/h, T*kB*E**((-COO_2s + CO_s + O_s)/(T*kB))/h]
        ref_kr_syms = [T*kB*E**((-s - CO_g + CO_s)/(T*kB))/h,
                       T*kB*E**((-2*s - O2_g + 2*O_s)/(T*kB))/h,
                       T*kB*E**((2*s - COO_2s + CO2_g)/(T*kB))/h]

        ret_kf_syms, ret_kr_syms = solver.get_rate_constant_syms()

        self.assertListEqual(ref_kf_syms, ret_kf_syms)
        self.assertListEqual(ref_kr_syms, ret_kr_syms)
        # }}}

    def test_get_equilibrium_constant_symbols(self):
        # {{{
        " Test we can get get correct equilibrium constants symbols. "
        # Construction.
        model = KineticModel(setup_file="input_files/solver_base.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        solver = model.solver()

        parser.parse_data(filename="input_files/rel_energy.py")
        solver.get_data()
        solver.get_data_symbols()

        # Get symbols.
        COO_2s = solver._extract_symbol("CO-O_2s", "free_energy")
        CO_s = solver._extract_symbol("CO_s", "free_energy")
        O_s = solver._extract_symbol("O_s", "free_energy")
        CO2_g = solver._extract_symbol("CO2_g", "free_energy")
        O2_g = solver._extract_symbol("O2_g", "free_energy")
        CO_g = solver._extract_symbol("CO_g", "free_energy")
        s = solver._extract_symbol("s", "free_energy")
        kB = solver._kB_sym
        T = solver._T_sym
        h = solver._h_sym
        from sympy import E

        kf_syms = [T*kB/h, T*kB/h, T*kB*E**((-COO_2s + CO_s + O_s)/(T*kB))/h]
        kr_syms = [T*kB*E**((-s - CO_g + CO_s)/(T*kB))/h,
                   T*kB*E**((-2*s - O2_g + 2*O_s)/(T*kB))/h,
                   T*kB*E**((2*s - COO_2s + CO2_g)/(T*kB))/h]

        ref_K = tuple([kf/kr for kf, kr in zip(kf_syms, kr_syms)])
        ret_K = solver.get_equilibrium_constant_syms()

        self.assertTupleEqual(ref_K, ret_K)
        # }}}

    def test_get_single_rate_sym(self):
        # {{{
        " Make sure we can get correct rate expression for an elementary reaction. "
        # Construction.
        model = KineticModel(setup_file="input_files/solver_base.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        solver = model.solver()

        parser.parse_data(filename="input_files/rel_energy.py")
        solver.get_data()
        solver.get_data_symbols()

        # Get symbols.

        # Free energy.
        G_COO_2s = solver._extract_symbol("CO-O_2s", "free_energy")
        G_CO_s = solver._extract_symbol("CO_s", "free_energy")
        G_O_s = solver._extract_symbol("O_s", "free_energy")
        G_CO2_g = solver._extract_symbol("CO2_g", "free_energy")
        G_O2_g = solver._extract_symbol("O2_g", "free_energy")
        G_CO_g = solver._extract_symbol("CO_g", "free_energy")
        G_s = solver._extract_symbol("s", "free_energy")

        # Coverage.
        c_CO_s = solver._extract_symbol("CO_s", "ads_cvg")
        c_O_s = solver._extract_symbol("O_s", "ads_cvg")
        c_s = solver._extract_symbol("s", "free_site_cvg")

        # Pressure.
        p_CO2_g = solver._extract_symbol("CO2_g", "pressure")
        p_O2_g = solver._extract_symbol("O2_g", "pressure")
        p_CO_g = solver._extract_symbol("CO_g", "pressure")

        # Constants.
        kB = solver._kB_sym
        T = solver._T_sym
        h = solver._h_sym
        from sympy import E

        kf = T*kB*E**((-G_COO_2s + G_CO_s + G_O_s)/(T*kB))/h
        kr = T*kB*E**((2*G_s - G_COO_2s + G_CO2_g)/(T*kB))/h

        rxn_expression = 'CO_s + O_s <-> CO-O_2s -> CO2_g + 2*_s'
        ref_rf = kf*c_CO_s*c_O_s
        ref_rr = kr*p_CO2_g*c_s**2

        ret_rf, ret_rr = solver.get_single_rate_sym(rxn_expression)

        self.assertEqual(ref_rf, ret_rf)
        self.assertEqual(ref_rr, ret_rr)
        # }}}

    def test_get_rate_syms(self):
        " Test we can get rate expressions correctly. "
        # Need Implimentation.

    def test_get_G_sub_dict(self):
        # {{{
        " Test private function _get_G_sub_dict(). "
        # Construction.
        model = KineticModel(setup_file="input_files/solver_base.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        solver = model.solver()

        parser.parse_data(filename="input_files/rel_energy.py")
        solver.get_data()
        solver.get_data_symbols()

        # Check.
        G_COO_2s = solver._extract_symbol("CO-O_2s", "free_energy")
        G_CO_s = solver._extract_symbol("CO_s", "free_energy")
        G_O_s = solver._extract_symbol("O_s", "free_energy")
        G_CO2_g = solver._extract_symbol("CO2_g", "free_energy")
        G_O2_g = solver._extract_symbol("O2_g", "free_energy")
        G_CO_g = solver._extract_symbol("CO_g", "free_energy")
        G_s = solver._extract_symbol("s", "free_energy")

        ref_dict = {G_O2_g: mpf('3.508000000002'),
                    G_CO2_g: mpf('0.0'),
                    G_CO_g: mpf('0.0'),
                    G_COO_2s: mpf('0.9259999999995'),
                    G_s: mpf('0.0'),
                    G_CO_s: mpf('-0.7580000000016'),
                    G_O_s: mpf('0.4340000000011')}

        ret_dict = solver._get_G_subs_dict()

        self.assertDictEqual(ref_dict, ret_dict)

    def test_get_theta_subs_dict(self):
        " Test protected function _get_theta_subs_dict(). "
        # Construction.
        model = KineticModel(setup_file="input_files/solver_base.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        solver = model.solver()

        parser.parse_data(filename="input_files/rel_energy.py")
        solver.get_data()
        solver.get_data_symbols()

        coverages = (0.5, 0.3)

        c_CO_s = solver._extract_symbol("CO_s", "ads_cvg")
        c_O_s = solver._extract_symbol("O_s", "ads_cvg")
        c_s = solver._extract_symbol("s", "free_site_cvg")

        ref_dict = {c_CO_s: 0.5, c_O_s: 0.3}
        ret_dict = solver._get_theta_subs_dict(coverages)
        
        self.assertDictEqual(ref_dict, ret_dict)
        # }}}

    def test_get_p_subs_dict(self):
        # {{{
        " Test protected function _get_p_subs_dict(). "
        # Construction.
        model = KineticModel(setup_file="input_files/solver_base.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        solver = model.solver()

        parser.parse_data(filename="input_files/rel_energy.py")
        solver.get_data()
        solver.get_data_symbols()

        p_CO2_g = solver._extract_symbol("CO2_g", "pressure")
        p_O2_g = solver._extract_symbol("O2_g", "pressure")
        p_CO_g = solver._extract_symbol("CO_g", "pressure")

        ref_dict = {p_O2_g: mpf('0.3333333333321'),
                    p_CO_g: mpf('1.0'),
                    p_CO2_g: mpf('0.0')}
        ret_dict = solver._get_p_subs_dict()

        self.assertDictEqual(ref_dict, ret_dict)
        # }}}

    def test_get_c_sub_dict(self):
        " Test protected function _get_c_sub_dict(). "
        # Need Implimentation.

    def test_get_subs_dict(self):
        # {{{
        " Make sure we can get correct substitution dict for all symbols. "
        # Construction.
        model = KineticModel(setup_file="input_files/solver_base.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        solver = model.solver()

        parser.parse_data(filename="input_files/rel_energy.py")
        solver.get_data()
        solver.get_data_symbols()

        # Get symbols.

        # Free energy.
        G_COO_2s = solver._extract_symbol("CO-O_2s", "free_energy")
        G_CO_s = solver._extract_symbol("CO_s", "free_energy")
        G_O_s = solver._extract_symbol("O_s", "free_energy")
        G_CO2_g = solver._extract_symbol("CO2_g", "free_energy")
        G_O2_g = solver._extract_symbol("O2_g", "free_energy")
        G_CO_g = solver._extract_symbol("CO_g", "free_energy")
        G_s = solver._extract_symbol("s", "free_energy")

        # Coverage.
        c_CO_s = solver._extract_symbol("CO_s", "ads_cvg")
        c_O_s = solver._extract_symbol("O_s", "ads_cvg")
        c_s = solver._extract_symbol("s", "free_site_cvg")

        # Pressure.
        p_CO2_g = solver._extract_symbol("CO2_g", "pressure")
        p_O2_g = solver._extract_symbol("O2_g", "pressure")
        p_CO_g = solver._extract_symbol("CO_g", "pressure")

        # Constants.
        kB = solver._kB_sym
        T = solver._T_sym
        h = solver._h_sym
        
        coverages = (0.5, 0.4)

        ref_dict = {G_O2_g: mpf('3.508000000002'),
                    c_CO_s: 0.5,
                    T: mpf('450.0'),
                    G_CO2_g: mpf('0.0'),
                    G_CO_g: mpf('0.0'),
                    G_COO_2s: mpf('0.9259999999995'),
                    p_O2_g: mpf('0.3333333333321'),
                    c_O_s: 0.4,
                    kB: mpf('8.617332400007e-5'),
                    h: mpf('4.135667662e-15'),
                    G_s: mpf('0.0'),
                    G_CO_s: mpf('-0.7580000000016'),
                    p_CO_g: mpf('1.0'),
                    G_O_s: mpf('0.4340000000011'),
                    p_CO2_g: mpf('0.0')}
        ret_dict = solver.get_subs_dict(coverages)

        self.assertDictEqual(ref_dict, ret_dict)
        # }}}

    def test_get_rate_constants_by_sym(self):
        " Make sure we can get rate constant correctly by symbols derivation. "
        # Construction.
        model = KineticModel(setup_file="input_files/solver_base.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        solver = model.solver()

        parser.parse_data(filename="input_files/rel_energy.py")
        solver.get_data()
        solver.get_data_symbols()

        # Check.
        ref_kfs = [mpf('9376477746560.0'), mpf('9376477746560.0'), mpf('0.09389759709029')]
        ref_krs = [mpf('30395.72540069'), mpf('2.542951527153e-17'), mpf('399.2961612232')]

        ret_kfs, ret_krs = solver.get_rate_constants_by_sym()

        self.assertListEqual(ref_kfs, ret_kfs)
        self.assertListEqual(ref_krs, ret_krs)

    def test_get_rates_by_syms(self):
        " Make sure we can get correct rates values by symbol derivation. "
        # Construction.
        model = KineticModel(setup_file="input_files/solver_base.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        solver = model.solver()

        parser.parse_data(filename="input_files/rel_energy.py")
        solver.get_data()
        solver.get_data_symbols()

        # Check.
        coverages = (0.5, 0.3)

        ref_rfs = (mpf('1875295549312.0'), mpf('125019703287.0'), mpf('0.01408463956352'))
        ref_rrs = (mpf('15197.86270034'), mpf('2.288656374422e-18'), mpf('0.0'))

        ret_rfs, ret_rrs = solver.get_rates_by_sym(cvgs_tuple=coverages)

        self.assertTupleEqual(ref_rfs, ret_rfs)
        self.assertTupleEqual(ref_rrs, ret_rrs)

    def test_get_net_rate_syms(self):
        " Make sure we can get correct net rate symbols for all elementary reactions. "
        # Need Implimentation.

    def test_get_net_rates_by_sym(self):
        " Test net rates calculating by symbol derivation. "
        # Construction.
        model = KineticModel(setup_file="input_files/solver_base.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        solver = model.solver()

        parser.parse_data(filename="input_files/rel_energy.py")
        solver.get_data()
        solver.get_data_symbols()

        # Check.
        coverages = (0.5, 0.3)

        ref_net_rates = (mpf('1875295534112.0'),
                         mpf('125019703287.0'),
                         mpf('0.01408463956352'))
        ret_net_rates = solver.get_net_rates_by_sym(coverages)

        self.assertTupleEqual(ref_net_rates, ret_net_rates)

    def test_get_tof_syms(self):
        " Test we can get TOF symbols correctly. "
        # NEED IMPLIMENTATION.

    def test_get_tof_by_sym(self):
        " Make sure we can get correct TOF value by symbols derivation. "
        # Construction.
        model = KineticModel(setup_file="input_files/solver_base.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        solver = model.solver()

        parser.parse_data(filename="input_files/rel_energy.py")
        solver.get_data()
        solver.get_data_symbols()

        # Check.
        coverages = (0.5, 0.3)

        ref_tof = (mpf('0.01408463956352'),
                   mpf('-1875295534112.0'),
                   mpf('-125019703287.0'))
        ret_tof = solver.get_tof_by_sym(coverages)

        self.assertTupleEqual(ref_tof, ret_tof)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(MeanFieldSolverTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


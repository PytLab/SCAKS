import logging
import re
import unittest

import numpy as np
from mpmath import mpf

from kynetix.models.micro_kinetic_model import MicroKineticModel
from kynetix.parsers.rxn_parser import *
from kynetix.solvers import *

from unit_test import *


class MeanFieldSolverTest(unittest.TestCase):

    def setUp(self):
        # Test case setting.
        # {{{
        self.maxDiff = None
        self.setup_dict = dict(
            rxn_expressions = [
                'CO_g + *_s -> CO_s',
                'O2_g + 2*_s -> 2O_s',
                'CO_s + O_s <-> CO-O_2s -> CO2_g + 2*_s',
            ],

            species_definitions = {
                'CO_g': {'pressure': 1.0},
                'O2_g': {'pressure': 1./3.},
                'CO2_g': {'pressure': 0.00},
                '*_s': {'site_name': '111', 'type': 'site', 'total': 1.0},
            },

            temperature = 450.0,
            parser = "RelativeEnergyParser",
            solver = "SteadyStateSolver",
            corrector = "ThermodynamicCorrector",
            plotter = "EnergyProfilePlotter",
            rootfinding = 'ConstrainedNewton',
            decimal_precision = 10,
            tolerance = 1e-20,
            max_rootfinding_iterations = 100,
        )

        self.abs_setup_dict = dict(
            rxn_expressions = [
                'CO_g + *_s -> CO_s',
                'O2_g + 2*_s -> 2O_s',
                'CO_s + O_s <-> CO-O_2s -> CO2_g + 2*_s',
            ],

            species_definitions = {
                'CO_g': {'pressure': 1.0},
                'O2_g': {'pressure': 1./3.},
                'CO2_g': {'pressure': 0.00},
                '*_s': {'site_name': '111', 'type': 'site', 'total': 1.0},
            },

            temperature = 450.0,
            parser = "AbsoluteEnergyParser",
            solver = "SteadyStateSolver",
            corrector = "ThermodynamicCorrector",
            plotter = "EnergyProfilePlotter",
            rootfinding = 'ConstrainedNewton',
            decimal_precision = 10,
            tolerance = 1e-20,
            max_rootfinding_iterations = 100,
        )
        # }}}

    def test_solver_construction_query(self):
        # {{{
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
        # }}}

    def test_get_data(self):
        # {{{
        " Test solver can get data correctly. "
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        parser = model.parser
        solver = model.solver

        # Parse data.
        parser.parse_data(filename=mkm_energy)
        solver.get_data()

        # Check pressure.
        ref_pressures = {'CO2_g': mpf('0.0'), 'CO_g': mpf('1.0'), 'O2_g': mpf(1./3.)}
        self.assertDictEqual(ref_pressures, solver.pressures)

        # Check concentrations.
        ref_concentrations = {}
        self.assertDictEqual(ref_concentrations, solver.concentrations)
        # }}}

    def test_get_rate_constants(self):
        # {{{
        " Make sure we can get rate constants correctly. "
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        parser = model.parser
        solver = model.solver

        self.assertRaises(AttributeError, solver.get_rate_constants)

        parser.parse_data(filename=mkm_energy)
        solver.get_data()

        # Check rate constants.
        ref_forward_rate_constants = [9376477746581.562, 9376477746581.562, 0.09389759708784133]
        ref_reverse_rate_constants = [30395.72540148798, 2.5429515269621107e-17, 399.2961612111053]
        ret_forward_rate_constants, ret_reverse_rate_constants = solver.get_rate_constants()
        self.assertListEqual(ref_forward_rate_constants, ret_forward_rate_constants)
        self.assertListEqual(ref_reverse_rate_constants, ret_reverse_rate_constants)
        # }}}

    def test_boltzmann_coverages(self):
        # {{{
        " Test we can get the Boltzmann converages. "
        # Construction.
        model = MicroKineticModel(setup_dict=self.abs_setup_dict,
                                  logger_level=logging.WARNING)
        parser = model.parser
        solver = model.solver

        # Parse absolute energies.
        parser.parse_data(filename=mkm_abs_energy)
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
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        parser = model.parser
        solver = model.solver

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
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        solver = model.solver

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
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        parser = model.parser
        solver = model.solver

        # Get data.
        parser.parse_data(filename=mkm_energy)
        solver.get_data()

        # Check.
        coverages = (0.5, 0.5)
        ret_rates = solver.get_rates(coverages)
        ref_rates = ((0.0, 0.0, 0.02347439927189),
                     (15197.8627007, 6.357378817428e-18, 0.0))
        for m in range(2):
            for n in range(3):
                self.assertAlmostEqual(ref_rates[m][n], ret_rates[m][n])

        # Check net rates.
        ref_net_rates = (-15197.8627007, -6.357378817429e-18, 0.02347439927189)
        ret_net_rates = solver.get_net_rates(coverages)
        for ref, ret in zip(ref_net_rates, ret_net_rates):
            self.assertAlmostEqual(ref, ret)
        # }}}

    def test_get_reversibilities(self):
        # {{{
        " Make sure we can get the correct reversibilities. "
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        parser = model.parser
        solver = model.solver

        parser.parse_data(filename=mkm_energy)
        solver.get_data()
        rfs, rrs = solver.get_rates((0.2, 0.5))

        ref_reversibilities = [2.1611331548919335e-09, 2.260045114764006e-29, 0.0]
        ret_reversibilities = solver.get_reversibilities(rfs, rrs)

        for ref, ret in zip(ref_reversibilities, ret_reversibilities):
            self.assertAlmostEqual(ref, ret)
        # }}}

    def test_get_tof(self):
        # {{{
        " Test we can get TOF correctly. "
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        parser = model.parser
        solver = model.solver

        parser.parse_data(filename=mkm_energy)
        solver.get_data()

        coverages = (0.2, 0.4)
        ref_tof = [0.007511807766946,
                   -3750591092544.0,
                   -500078813148.0]
        ret_tof = solver.get_tof(coverages)

        for ref, ret in zip(ref_tof, ret_tof):
            self.assertAlmostEqual(ref, float(ret))

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
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        solver = model.solver

        solver.get_data_symbols()

        # Check P symbols.
        ref_p_symbols = ('p_CO2_g', 'p_CO_g', 'p_O2_g')
        ret_p_symbols = solver._p_sym

        for p, p_str in zip(ret_p_symbols, ref_p_symbols):
            self.assertEqual(p.name, p_str)

        # Check concentration symbols.
        ref_c_symbols = ()
        ret_c_symbols = solver._c_sym

        for c, c_str in zip(ret_p_symbols, ref_p_symbols):
            self.assertEqual(c.name, c_str)

        # Check adsorbate coverage symbols.
        ref_ads_theta_symbols = ("theta_CO_s", "theta_O_s")
        ret_ads_theta_symbols = solver._ads_theta_sym

        for ads_theta, ads_theta_str in zip(ret_ads_theta_symbols, ref_ads_theta_symbols):
            self.assertEqual(ads_theta.name, ads_theta_str)

        # Check free site coverage symbols.
        self.assertEqual(1.0, solver._fsite_theta_sym[0] + sum(solver._ads_theta_sym))

        # Check free energy coverage symbols.
        ref_Ga_symbols = ("Ga_0", "Ga_1", "Ga_2")
        ret_Ga_symbols = solver._Ga_sym

        for Ga, Ga_str in zip(ret_Ga_symbols, ref_Ga_symbols):
            self.assertEqual(Ga.name, Ga_str)

        ref_dG_symbols = ("dG_0", "dG_1", "dG_2")
        ret_dG_symbols = solver._dG_sym
        for dG, dG_str in zip(ret_dG_symbols, ref_dG_symbols):
            self.assertEqual(dG.name, dG_str)

        # Check K coverage symbols.
        ref_K_symbols = ("K_0", "K_1", "K_2")
        ret_K_symbols = solver._K_sym

        for K, K_str in zip(ret_K_symbols, ref_K_symbols):
            self.assertEqual(K.name, K_str)
        # }}}

    def test_extract_symbol(self):
        # {{{
        " Test protected function _extract_symbol(). "
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        solver = model.solver

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
        ret_cvg = solver._extract_symbol('*_s', 'free_site_cvg')
        CO = solver._extract_symbol('CO_s', 'ads_cvg')
        O = solver._extract_symbol('O_s', 'ads_cvg')
        ref_cvg = 1.0 - CO - O
        self.assertEqual(ret_cvg, ref_cvg)

        # }}}

    def test_get_rate_constant_symbols(self):
        # {{{
        " Test we can get get correct rate constants symbols. "
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict,
                                  logger_level=logging.WARNING)
        parser = model.parser
        solver = model.solver

        parser.parse_data(filename=mkm_energy)
        solver.get_data()
        solver.get_data_symbols()

        # Get symbols.
        kB = solver._kB_sym
        T = solver._T_sym
        h = solver._h_sym
        Ga_0, Ga_1, Ga_2 = solver._Ga_sym
        dG_0, dG_1, dG_2 = solver._dG_sym
        from sympy import E

        ref_kf_syms = [T*kB*E**(-Ga_0/(T*kB))/h,
                       T*kB*E**(-Ga_1/(T*kB))/h,
                       T*kB*E**(-Ga_2/(T*kB))/h]
        ref_kr_syms = [T*kB*E**((-Ga_0 + dG_0)/(T*kB))/h,
                       T*kB*E**((-Ga_1 + dG_1)/(T*kB))/h,
                       T*kB*E**((-Ga_2 + dG_2)/(T*kB))/h]

        ret_kf_syms, ret_kr_syms = solver.get_rate_constant_syms()

        self.assertListEqual(ref_kf_syms, ret_kf_syms)
        self.assertListEqual(ref_kr_syms, ret_kr_syms)
        # }}}

    def test_get_equilibrium_constant_symbols(self):
        # {{{
        " Test we can get get correct equilibrium constants symbols. "
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        parser = model.parser
        solver = model.solver

        parser.parse_data(filename=mkm_energy)
        solver.get_data()
        solver.get_data_symbols()

        # Get symbols.
        Ga_0, Ga_1, Ga_2 = solver._Ga_sym
        dG_0, dG_1, dG_2 = solver._dG_sym
        kB = solver._kB_sym
        T = solver._T_sym
        h = solver._h_sym
        from sympy import E

        ref_K = (E**(-Ga_0/(T*kB))*E**(-(-Ga_0 + dG_0)/(T*kB)),
                 E**(-Ga_1/(T*kB))*E**(-(-Ga_1 + dG_1)/(T*kB)),
                 E**(-Ga_2/(T*kB))*E**(-(-Ga_2 + dG_2)/(T*kB)))
        ret_K = solver.get_equilibrium_constant_syms()

        for ref, ret in zip(ref_K, ret_K):
            self.assertEqual(ref.simplify(), ret.simplify())
        # }}}

    def test_get_single_rate_sym(self):
        # {{{
        " Make sure we can get correct rate expression for an elementary reaction. "
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        parser = model.parser
        solver = model.solver

        parser.parse_data(filename=mkm_energy)
        solver.get_data()
        solver.get_data_symbols()

        # Get symbols.

        # Coverage.
        theta_CO_s = solver._extract_symbol("CO_s", "ads_cvg")
        theta_O_s = solver._extract_symbol("O_s", "ads_cvg")
        theta_s = solver._extract_symbol("*_s", "free_site_cvg")

        # Pressure.
        p_CO2_g = solver._extract_symbol("CO2_g", "pressure")
        p_O2_g = solver._extract_symbol("O2_g", "pressure")
        p_CO_g = solver._extract_symbol("CO_g", "pressure")

        Ga_0, Ga_1, Ga_2 = solver._Ga_sym
        dG_0, dG_1, dG_2 = solver._dG_sym

        # Constants.
        kB = solver._kB_sym
        T = solver._T_sym
        h = solver._h_sym
        from sympy import E

        rxn_expression = 'CO_s + O_s <-> CO-O_2s -> CO2_g + 2*_s'
        ref_rf = T*kB*theta_CO_s*theta_O_s*E**(-Ga_2/(T*kB))/h
        ref_rr = T*kB*p_CO2_g*(-theta_CO_s - theta_O_s + 1.0)**2*E**((-Ga_2 + dG_2)/(T*kB))/h

        ret_rf, ret_rr = solver.get_single_rate_sym(rxn_expression)

        self.assertEqual(ref_rf.simplify(), ret_rf.simplify())
        self.assertEqual(ref_rr.simplify(), ret_rr.simplify())
        # }}}

    def test_get_rate_syms(self):
        " Test we can get rate expressions correctly. "
        # {{{
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        parser = model.parser
        solver = model.solver

        parser.parse_data(filename=mkm_energy)
        solver.get_data()
        solver.get_data_symbols()

        # Get symbols.

        # Coverage.
        theta_CO_s = solver._extract_symbol("CO_s", "ads_cvg")
        theta_O_s = solver._extract_symbol("O_s", "ads_cvg")
        theta_s = solver._extract_symbol("*_s", "free_site_cvg")

        # Pressure.
        p_CO2_g = solver._extract_symbol("CO2_g", "pressure")
        p_O2_g = solver._extract_symbol("O2_g", "pressure")
        p_CO_g = solver._extract_symbol("CO_g", "pressure")

        Ga_0, Ga_1, Ga_2 = solver._Ga_sym
        dG_0, dG_1, dG_2 = solver._dG_sym

        # Constants.
        kB = solver._kB_sym
        T = solver._T_sym
        h = solver._h_sym
        from sympy import E

        ref_forwards = [
            T*kB*p_CO_g*(-theta_CO_s - theta_O_s + 1.0)*E**(-Ga_0/(T*kB))/h,
            T*kB*p_O2_g*(-theta_CO_s - theta_O_s + 1.0)**2*E**(-Ga_1/(T*kB))/h,
            T*kB*theta_CO_s*theta_O_s*E**(-Ga_2/(T*kB))/h
        ]

        ref_reverses = [
            T*kB*theta_CO_s*E**((-Ga_0 + dG_0)/(T*kB))/h,
            T*kB*theta_O_s**2*E**((-Ga_1 + dG_1)/(T*kB))/h,
            T*kB*p_CO2_g*(-theta_CO_s - theta_O_s + 1.0)**2*E**((-Ga_2 + dG_2)/(T*kB))/h
        ]

        ret_forwards, ret_reverses = solver.get_rate_syms()

        for ref, ret in zip(ref_forwards, ret_forwards):
            self.assertEqual(ref.simplify(), ret.simplify())

        for ref, ret in zip(ref_reverses, ret_reverses):
            self.assertEqual(ref.simplify(), ret.simplify())
        # }}}

    def test_get_G_sub_dict(self):
        # {{{
        " Test private function _get_G_sub_dict(). "
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        parser = model.parser
        solver = model.solver

        parser.parse_data(filename=mkm_energy)
        solver.get_data()
        solver.get_data_symbols()

        # Check.
        Ga_0, Ga_1, Ga_2 = solver._Ga_sym
        dG_0, dG_1, dG_2 = solver._dG_sym

        ref_dict = {Ga_0: 0.0,
                    dG_1: -2.64,
                    Ga_1: 0.0,
                    dG_0: -0.758,
                    Ga_2: 1.25,
                    dG_2: 0.324}

        ret_dict = solver._get_G_subs_dict()

        for key in ref_dict:
            self.assertEqual(ref_dict[key], ret_dict[key])

    def test_get_theta_subs_dict(self):
        " Test protected function _get_theta_subs_dict(). "
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        parser = model.parser
        solver = model.solver

        parser.parse_data(filename=mkm_energy)
        solver.get_data()
        solver.get_data_symbols()

        coverages = (0.5, 0.3)

        c_CO_s = solver._extract_symbol("CO_s", "ads_cvg")
        c_O_s = solver._extract_symbol("O_s", "ads_cvg")
        c_s = solver._extract_symbol("*_s", "free_site_cvg")

        ref_dict = {c_CO_s: 0.5, c_O_s: 0.3}
        ret_dict = solver._get_theta_subs_dict(coverages)
        
        self.assertDictEqual(ref_dict, ret_dict)
        # }}}

    def test_get_p_subs_dict(self):
        # {{{
        " Test protected function _get_p_subs_dict(). "
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        parser = model.parser
        solver = model.solver

        parser.parse_data(filename=mkm_energy)
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
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        parser = model.parser
        solver = model.solver

        parser.parse_data(filename=mkm_energy)
        solver.get_data()
        solver.get_data_symbols()

        # Get symbols.

        # Free energy.
        Ga_0, Ga_1, Ga_2 = solver._Ga_sym
        dG_0, dG_1, dG_2 = solver._dG_sym

        G_dict = {Ga_0: 0.0,
                  dG_1: -2.64,
                  Ga_1: 0.0,
                  dG_0: -0.758,
                  Ga_2: 1.25,
                  dG_2: 0.324}

        # Coverage.
        c_CO_s = solver._extract_symbol("CO_s", "ads_cvg")
        c_O_s = solver._extract_symbol("O_s", "ads_cvg")
        c_s = solver._extract_symbol("*_s", "free_site_cvg")

        # Pressure.
        p_CO2_g = solver._extract_symbol("CO2_g", "pressure")
        p_O2_g = solver._extract_symbol("O2_g", "pressure")
        p_CO_g = solver._extract_symbol("CO_g", "pressure")

        # Constants.
        kB = solver._kB_sym
        T = solver._T_sym
        h = solver._h_sym
        
        coverages = (0.5, 0.4)

        ref_dict = {c_CO_s: 0.5,
                    T: mpf('450.0'),
                    p_O2_g: mpf('0.3333333333321'),
                    c_O_s: 0.4,
                    kB: mpf('8.617332400007e-5'),
                    h: mpf('4.135667662e-15'),
                    Ga_0: 0.0,
                    dG_1: -2.64,
                    Ga_1: 0.0,
                    dG_0: -0.758,
                    Ga_2: 1.25,
                    dG_2: 0.324}
        ret_dict = solver.get_subs_dict(coverages)

        for key in ref_dict:
            self.assertAlmostEqual(float(ref_dict[key]), float(ret_dict[key]))
        # }}}

    def test_get_rate_constants_by_sym(self):
        " Make sure we can get rate constant correctly by symbols derivation. "
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        parser = model.parser
        solver = model.solver

        parser.parse_data(filename=mkm_energy)
        solver.get_data()
        solver.get_data_symbols()

        # Check.
        ref_kfs = [mpf('9376477746560.0'), mpf('9376477746560.0'), mpf('0.09389759709029')]
        ref_krs = [mpf('30395.72540212'), mpf('2.542951527113e-17'), mpf('399.2961612195')]

        ret_kfs, ret_krs = solver.get_rate_constants_by_sym()

        self.assertListEqual(ref_kfs, ret_kfs)
        self.assertListEqual(ref_krs, ret_krs)

    def test_get_rates_by_syms(self):
        " Make sure we can get correct rates values by symbol derivation. "
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        parser = model.parser
        solver = model.solver

        parser.parse_data(filename=mkm_energy)
        solver.get_data()
        solver.get_data_symbols()

        # Check.
        coverages = (0.5, 0.3)

        ref_rfs = (mpf('1875295549312.0'), mpf('125019703287.0'), mpf('0.01408463956352'))
        ref_rrs = (mpf('15197.86270106'), mpf('2.288656374397e-18'), mpf('0.0'))

        ret_rfs, ret_rrs = solver.get_rates_by_sym(cvgs_tuple=coverages)

        self.assertTupleEqual(ref_rfs, ret_rfs)
        self.assertTupleEqual(ref_rrs, ret_rrs)

    def test_get_net_rate_syms(self):
        " Make sure we can get correct net rate symbols for all elementary reactions. "
        # {{{
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        parser = model.parser
        solver = model.solver

        parser.parse_data(filename=mkm_energy)
        solver.get_data()
        solver.get_data_symbols()

        # Free energy.
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

        ref_net_rates = [
            (T*kB*p_CO_g*(-theta_CO_s - theta_O_s + 1.0)*E**(-Ga_0/(T*kB))/h - 
             T*kB*theta_CO_s*E**((-Ga_0 + dG_0)/(T*kB))/h),
            (T*kB*p_O2_g*(-theta_CO_s - theta_O_s + 1.0)**2*E**(-Ga_1/(T*kB))/h -
             T*kB*theta_O_s**2*E**((-Ga_1 + dG_1)/(T*kB))/h),
            (-T*kB*p_CO2_g*(-theta_CO_s - theta_O_s + 1.0)**2*E**((-Ga_2 + dG_2)/(T*kB))/h +
             T*kB*theta_CO_s*theta_O_s*E**(-Ga_2/(T*kB))/h)
        ]
        ret_net_rates = solver.get_net_rate_syms()

        for ref, ret in zip(ref_net_rates, ret_net_rates):
            self.assertEqual(ref.simplify(), ret.simplify())
        # }}}

    def test_get_net_rates_by_sym(self):
        " Test net rates calculating by symbol derivation. "
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        parser = model.parser
        solver = model.solver

        parser.parse_data(filename=mkm_energy)
        solver.get_data()
        solver.get_data_symbols()

        # Check.
        coverages = (0.5, 0.3)

        ref_net_rates = (mpf('1875295534112.0'),
                         mpf('125019703287.0'),
                         mpf('0.01408463956352'))
        ret_net_rates = solver.get_net_rates_by_sym(coverages)

        for ref, ret in zip(ref_net_rates, ret_net_rates):
            self.assertAlmostEqual(ref, ret)

    def test_get_tof_syms(self):
        " Test we can get TOF symbols correctly. "
        # {{{
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        parser = model.parser
        solver = model.solver

        parser.parse_data(filename=mkm_energy)
        solver.get_data()
        solver.get_data_symbols()

        # Free energy.
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

        ref_tofs = (
            (-1.0*T*kB*p_CO2_g*(-theta_CO_s - theta_O_s + 1.0)**2*E**((-Ga_2 + dG_2)/(T*kB))/h +
             1.0*T*kB*theta_CO_s*theta_O_s*E**(-Ga_2/(T*kB))/h),
            (-1.0*T*kB*p_CO_g*(-theta_CO_s - theta_O_s + 1.0)*E**(-Ga_0/(T*kB))/h +
             1.0*T*kB*theta_CO_s*E**((-Ga_0 + dG_0)/(T*kB))/h),
            (-1.0*T*kB*p_O2_g*(-theta_CO_s - theta_O_s + 1.0)**2*E**(-Ga_1/(T*kB))/h +
             1.0*T*kB*theta_O_s**2*E**((-Ga_1 + dG_1)/(T*kB))/h)
        )

        ret_tofs = solver.get_tof_syms()

        for ref, ret in zip(ref_tofs, ret_tofs):
            self.assertEqual(ref.simplify(), ret.simplify())
        # }}}

    def test_get_tof_by_sym(self):
        " Make sure we can get correct TOF value by symbols derivation. "
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        parser = model.parser
        solver = model.solver

        parser.parse_data(filename=mkm_energy)
        solver.get_data()
        solver.get_data_symbols()

        # Check.
        coverages = (0.5, 0.3)

        ref_tof = (mpf('0.01408463956352'),
                   mpf('-1875295534112.0'),
                   mpf('-125019703287.0'))
        ret_tof = solver.get_tof_by_sym(coverages)

        for ref, ret in zip(ref_tof, ret_tof):
            self.assertAlmostEqual(ref, ret)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(MeanFieldSolverTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


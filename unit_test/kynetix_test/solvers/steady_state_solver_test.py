import logging
import re
import unittest

import numpy as np
from numpy import matrix
import mpmath as mp
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

#    def test_get_data(self):
#        " Test solver can get data correctly. "
#        # Construction.
#        model = KineticModel(setup_file="input_files/steady_state_solver.mkm",
#                             verbosity=logging.WARNING)
#        parser = model.parser()
#        solver = model.solver()
#
#        # Get data before parsing data, an exception would be expected.
#        self.assertRaises(IOError, solver.get_data)
#
#        # Parse data.
#        parser.parse_data(relative=True, filename="input_files/rel_energy.py")
#        solver.get_data()
#
#        self.assertTrue(solver.has_relative_energy())
#
#        # Check pressure.
#        ref_pressures = {'CO2_g': mpf('0.0'), 'CO_g': mpf('1.0'), 'O2_g': mpf(1./3.)}
#        self.assertDictEqual(ref_pressures, solver.pressures())
#
#        # Check concentrations.
#        ref_concentrations = {}
#        self.assertDictEqual(ref_concentrations, solver.concentrations())
#
#        # Relative energies.
#        ref_relative_energies = {'Gaf': [0.0, 0.0, 1.25],
#                                 'Gar': [0.758, 2.64, 0.9259999999999999],
#                                 'dG': [-0.758, -2.64, 0.324]}
#        self.assertDictEqual(ref_relative_energies, solver.relative_energies())
#
#        # Formation energies.
#        self.assertRaises(AttributeError, solver.formation_energies)
#
#        # Parse absolute energies.
#        parser.parse_data(relative=False, filename="input_files/rel_energy.py")
#        solver.get_data()
#
#        self.assertTrue(solver.has_absolute_energy())
#
#        # Check formation energies.
#        ref_formation_energies = {'*_s': mpf('0.0'),
#                                  'CO-O_2s': mpf('0.9259999999999999342747969421907328069210052490234375'),
#                                  'CO2_g': mpf('0.0'),
#                                  'CO_g': mpf('0.0'),
#                                  'CO_s': mpf('-0.75800000000000000710542735760100185871124267578125'),
#                                  'O2_g': mpf('3.50800000000000000710542735760100185871124267578125'),
#                                  'O_s': mpf('0.4339999999999999413802242997917346656322479248046875')}
#        self.assertDictEqual(ref_formation_energies, solver.formation_energies())

    def test_elementary_dtheta_dt_expression(self):
        " Test get_elementary_dtheta_dt_expression() function. "
        # Construction.
        model = KineticModel(setup_file="input_files/steady_state_solver.mkm",
                             verbosity=logging.WARNING)
        solver = model.solver()

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
        ref_dtheta_dt = (mpf('2812943317895.314696341954248072244057093672649102054657986641525702269301823174878392219699918541664374'),
                         mpf('562588664794.8844816900751974899822236555742319376607896488580965637603379490166467802314542723252088763'))
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

    def test_total_term_adsorbate_derivation(self):
        " Test private function __total_term_adsorbate_derivation(). "
        # Construction.
        model = KineticModel(setup_file="input_files/steady_state_solver.mkm",
                             verbosity=logging.WARNING)
        solver = model.solver()

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
        model = KineticModel(setup_file="input_files/steady_state_solver.mkm",
                             verbosity=logging.WARNING)
        solver = model.solver()

        adsorbate = "CO_s"
        poly_expression = ("dtheta_dt[0] = kf[0]*p['CO_g']*theta['*_s'] - " +
                           "kr[0]*theta['CO_s'] + kr[2]*p['CO2_g']*theta['*_s']**2 - " +
                           "kf[2]*theta['CO_s']*theta['O_s']")
        ref_expression = ("-kf[0]*p['CO_g'] - kr[0] + " +
                          "-2*kr[2]*p['CO2_g']*(1.0 - theta['CO_s'] - theta['O_s'])**1 - " +
                          "kf[2]*theta['O_s']")
        ret_expression = solver.poly_adsorbate_derivation(adsorbate, poly_expression)
        self.assertEqual(ref_expression, ret_expression)

        adsorbate = "O_s"
        poly_expression = ("dtheta_dt[1] = 2*kf[1]*p['O2_g']*theta['*_s']**2 - " +
                           "2*kr[1]*theta['O_s']**2 + kr[2]*p['CO2_g']*theta['*_s']**2 - " +
                           "kf[2]*theta['CO_s']*theta['O_s']")
        ref_expression = ("-2*2*kf[1]*p['O2_g']*(1.0 - theta['CO_s'] - theta['O_s'])**1 - " +
                          "2*2*kr[1]*theta['O_s']**1 + " +
                          "-2*kr[2]*p['CO2_g']*(1.0 - theta['CO_s'] - theta['O_s'])**1 - " +
                          "kf[2]*theta['CO_s']")
        ret_expression = solver.poly_adsorbate_derivation(adsorbate, poly_expression)
        self.assertEqual(ref_expression, ret_expression)

    def test_analytical_jacobian(self):
        " Make sure we can get analytical Jacobian matrix correctly. "
        # Construction.
        model = KineticModel(setup_file="input_files/steady_state_solver.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        solver = model.solver()

        parser.parse_data(filename="input_files/rel_energy.py")
        solver.get_data()

        # Check.
        coverages = (0.2, 0.4)
        ref_jacobian = \
        [[mpf('-9376477776977.325460526815116254348794135605545524270101419146551990140566123176340401768058541165373754'),
          mpf('-9376477746581.581279519417568267191072063225136043385050709573275995070283061588170200884029270582686877')],
         [mpf('-5000788131510.204225705501803201033400924944119994425853606776028424154830007847261086340040492359548807'),
          mpf('-5000788131510.18544618608423497452955329311275738716732389010571374197652467107703380508721422584175113')]]    
        ret_jacobian = solver.analytical_jacobian(coverages).tolist()
        self.assertListEqual(ref_jacobian, ret_jacobian)

    def test_get_residual(self):
        " Test we can get correct residual. "
        # Construction.
        model = KineticModel(setup_file="input_files/steady_state_solver.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        solver = model.solver()

        parser.parse_data(filename="input_files/rel_energy.py")
        solver.get_data()

        coverages = (0.99, 0.01)
        ref_residual = mpf('30091.76907705931107958367882296829850017413372279726741798085781696430451251481888596068989202277621021')
        ret_residual = solver.get_residual(coverages)

        self.assertEqual(ref_residual, ret_residual)

    def test_coarse_steady_state_cvgs(self):
        " Make sure we can get a coarse coverages. "
        # Construction.
        model = KineticModel(setup_file="input_files/steady_state_solver.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        solver = model.solver()

        parser.parse_data(filename="input_files/rel_energy.py")
        solver.get_data()

        coverages = (0.99, 0.0)
        ref_coarse_cvgs = [0.9993009034360529, 0.0006990933245135515]
        ret_coarse_cvgs = solver.coarse_steady_state_cvgs(coverages).tolist()

        self.assertListEqual(ref_coarse_cvgs, ret_coarse_cvgs)

    def test_get_steady_state_coverages(self):
        " Test we can get correct steady state coverages. "
        # Construction.
        model = KineticModel(setup_file="input_files/steady_state_solver.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        solver = model.solver()

        parser.parse_data(filename="input_files/rel_energy.py")
        solver.get_data()
        
        # Check.
        coverages = solver.boltzmann_coverages()
        ref_sscvg = (mpf('0.9993009023313341077071690952744161899993498376199216404808490009152529045923214583626561308597735533188'),
                     mpf('0.0006990944292324176713411953522859201607970088893776515794721879333704919691278789195839768935388121262119'))
        ret_sscvg = solver.get_steady_state_cvgs(coverages)
        self.assertTupleEqual(ref_sscvg, ret_sscvg)

        # Check error.
        ref_error = mpf('2.915131960917183181735290581323356891541606225932002083392991290212256934502240148437948718330981045238e-38')
        ret_error = solver.error()
        self.assertEqual(ret_error, ref_error)

    def test_get_intermediates_Gs(self):
        " Test private function __get_intermediates_Gs(). "
        # Construction.
        model = KineticModel(setup_file="input_files/steady_state_solver.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        solver = model.solver()

        parser.parse_data(filename="input_files/rel_energy.py")
        solver.get_data()

        # Check.
        ref_Gs = [mpf('-0.75800000000000000710542735760100185871124267578125'),
                  mpf('0.4339999999999999413802242997917346656322479248046875'),
                  mpf('0.9259999999999999342747969421907328069210052490234375')]
        ret_Gs = solver._SteadyStateSolver__get_intermediates_Gs()

        self.assertListEqual(ref_Gs, ret_Gs)

    def test_get_Gs_tof(self):
        " Test private function __get_Gs_tof(). "
        # Construction.
        model = KineticModel(setup_file="input_files/steady_state_solver.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        solver = model.solver()

        parser.parse_data(filename="input_files/rel_energy.py")
        solver.get_data()

        # Get steady state converages first.
        coverages = solver.boltzmann_coverages()
        solver.get_steady_state_cvgs(coverages)

        # Check.
        Gs = solver._SteadyStateSolver__get_intermediates_Gs()
        ref_tof = [mpf('0.0000655973959734845716839853127464104731569988368580121920613237008751507516356019917597493208387158845367'),
                   mpf('-0.00006559739597348503704886124510603927478715353925738869127523791123005360979166894188698600090288385521978'),
                   mpf('-0.00003279869798674251852443062255301963739357676962869434563761895561502680489741843434563928385867591142223')]
        ret_tof = solver._SteadyStateSolver__get_Gs_tof(Gs)

        self.assertListEqual(ref_tof, ret_tof)

    def test_get_Gs_gas_tof(self):
        " Test private function __get_Gs_tof(). "
        # Construction.
        model = KineticModel(setup_file="input_files/steady_state_solver.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        solver = model.solver()

        parser.parse_data(filename="input_files/rel_energy.py")
        solver.get_data()

        # Get steady state converages first.
        coverages = solver.boltzmann_coverages()
        solver.get_steady_state_cvgs(coverages)

        # Check gas tof.
        Gs = solver._SteadyStateSolver__get_intermediates_Gs()
        ref_tof = mpf('0.0000655973959734845716839853127464104731569988368580121920613237008751507516356019917597493208387158845367')
        ret_tof = solver._SteadyStateSolver__get_Gs_tof(Gs, gas_name="CO2_g")
        self.assertEqual(ref_tof, ret_tof)

    def test_get_single_XTRC(self):
        " Test function get_single_XTRC(). "
        # Construction.
        model = KineticModel(setup_file="input_files/steady_state_solver.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        solver = model.solver()

        parser.parse_data(filename="input_files/rel_energy.py")
        solver.get_data()

        # Get steady state converages first.
        coverages = solver.boltzmann_coverages()
        solver.get_steady_state_cvgs(coverages)

        # Check.
        gas_name = "CO2_g"
        ref_XTRC = [mpf('-1.140775679060479546680777929992790525034725825880253939251504297700926329015933617537522978478527774183'),
                    mpf('-1.140775679060479546680777929992790525034725825880253939251504297700926329015933617537522978478527774183'),
                    mpf('0.881465200934955213215196217181568318297209213874500010342977516784412704475283998243820141954241091872')]
        ret_XTRC = solver.get_single_XTRC(gas_name)
        self.assertListEqual(ref_XTRC, ret_XTRC)

    def test_get_XTRC(self):
        " Test function get_XTRC(). "
        # NEED IMPLIMENTATION.

    def test_get_single_XRC(self):
        " Test function get_single_XRC(). "
        # Construction.
        model = KineticModel(setup_file="input_files/steady_state_solver.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        solver = model.solver()

        parser.parse_data(filename="input_files/rel_energy.py")
        solver.get_data()

        # Get steady state converages first.
        coverages = solver.boltzmann_coverages()
        solver.get_steady_state_cvgs(coverages)

        # Check.
        gas_name = "CO2_g"
        ref_XRC = [mpf('-0.000000004330301262036152923211709673711721112213085538797673380141841669630093740288896353967710982069299592622'),
                   mpf('0.9986021755886114784434712184789699872827416819468954527192794077349320805208548704540326869330347344589'),
                   mpf('0.001398549092163431169588440111213386822318548225009429894534084532021244154191984863118384815049813695996')]
        ret_XRC = solver.get_single_XRC(gas_name, epsilon=1e-5)
        self.assertListEqual(ref_XRC, ret_XRC)

#    def test_get_single_XRC_multi_thread(self):
#        " Test function get_single_XRC() in multi-threads. "
#        # Construction.
#        model = KineticModel(setup_file="input_files/steady_state_solver.mkm",
#                             verbosity=logging.WARNING)
#        parser = model.parser()
#        solver = model.solver()
#
#        parser.parse_data(filename="input_files/rel_energy.py")
#        solver.get_data()
#
#        # Get steady state converages first.
#        coverages = solver.boltzmann_coverages()
#        solver.get_steady_state_cvgs(coverages)
#
#        # Check.
#        gas_name = "CO2_g"
#        ref_XRC = [mpf('-0.000000004330301262036152923211709673711721112213085538797673380141841669630093740288896353967710982069299592622'),
#                   mpf('0.9986021755886114784434712184789699872827416819468954527192794077349320805208548704540326869330347344589'),
#                   mpf('0.001398549092163431169588440111213386822318548225009429894534084532021244154191984863118384815049813695996')]
#        ret_XRC = solver.get_single_XRC(gas_name, epsilon=1e-5, parallel=True)
#        for ref, ret in zip(ref_XRC, ret_XRC):
#            self.assertAlmostEqual(ref, ret)

    def test_get_elementary_dtheta_dt_sym(self):
        " Test we can get correct dtheta/dt expression for an elementary reaction. "
        # Construction.
        model = KineticModel(setup_file="input_files/steady_state_solver.mkm",
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
        rf = kf*c_CO_s*c_O_s
        rr = kr*p_CO2_g*c_s**2

        ref_dtheta_dt = rr - rf
        adsorbate = "CO_s"
        ret_dtheta_dt = solver.get_elementary_dtheta_dt_sym(adsorbate, rxn_expression)

        self.assertEqual(ref_dtheta_dt, ret_dtheta_dt)

    def test_get_adsorbate_dtheta_dt_sym(self):
        " Test we can get correct dtheta/dt for an adsorbate. "
        # NEED IMPLIMENTATION.

    def test_steady_state_function_by_sym(self):
        " Test function steady_state_function(). "
        # Construction.
        model = KineticModel(setup_file="input_files/steady_state_solver.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        solver = model.solver()

        parser.parse_data(filename="input_files/rel_energy.py")
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
        model = KineticModel(setup_file="input_files/steady_state_solver.mkm",
                             verbosity=logging.WARNING)
        parser = model.parser()
        solver = model.solver()

        parser.parse_data(filename="input_files/rel_energy.py")
        solver.get_data()
        solver.get_data_symbols()

        # Check.
        coverages = (0.5, 0.3)

        ref_jacobian = [[mpf('-9376477776977.31640625'), mpf('-9376477746581.609375')],
                        [mpf('-1250197032877.56982421875'), mpf('-1250197032877.588623046875')]]
        ret_jacobian = solver.analytical_jacobian_by_sym(coverages).tolist()

        self.assertListEqual(ref_jacobian, ret_jacobian)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(SteadyStateSolverTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


import logging
import re
import unittest

import numpy as np
from mpmath import mpf

from kynetix.models.micro_kinetic_model import MicroKineticModel
from kynetix.correctors import *

from unit_test import *


class ThermodynamicCorrectorTest(unittest.TestCase):

    def setUp(self):
        # Test case setting.
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
            decimal_precision = 100,
            tolerance = 1e-20,
            max_rootfinding_iterations = 100,
        )

    def test_construction_and_query(self):
        " Test plotter construction and query. "
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        corrector = model.corrector

        self.assertTrue(isinstance(corrector, ThermodynamicCorrector))

    def test_shomate_correction(self):
        " Test we can get correct energy correction by Shomate equation. "
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.ERROR)
        corrector = model.corrector
        # Check.
        gas = "CO_g"
        ret_delta = corrector.shomate_correction(gas)
        ref_delta = -0.84190990132027166
        self.assertEqual(ret_delta, ref_delta)

        gas = "O2_g"
        ret_delta = corrector.shomate_correction(gas)
        ref_delta = -0.87627116516400394
        self.assertEqual(ret_delta, ref_delta)

        gas = "O4_g"
        # A warning would be expected.
        self.assertEqual(corrector.shomate_correction(gas), 0.0)

        species = "O-O_s"
        self.assertEqual(corrector.shomate_correction(species), 0.0)

    def test_entropy_correction(self):
        " Make sure we can get correct entropy correction. "
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.ERROR)
        corrector = model.corrector

        # Check.
        gas = "CO_g"
        ret_delta = corrector.entropy_correction(gas)
        ref_delta = -1.1538116935108251
        self.assertEqual(ref_delta, ret_delta)

        gas = "O2_g"
        ret_delta = corrector.entropy_correction(gas)
        ref_delta = -1.2250150716175705
        self.assertEqual(ref_delta, ret_delta)

        gas = "O4_g"
        # A warning would be expected.
        self.assertEqual(corrector.entropy_correction(gas), 0.0)

        species = "O-O_s"
        self.assertEqual(corrector.entropy_correction(species), 0.0)

#    def test_solvers_correction_energy(self):
#        " Test solver's correction energy function. "
#        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
#        parser = model.parser
#        parser.parse_data(filename=mkm_energy)
#        solver = model.solver
#        solver.get_data()
#
#        ref_e = {'*_s': mpf('0.0'),
#                 'CO-O_2s': mpf('0.9259999999999999342747969421907328069210052490234375'),
#                 'CO2_g': mpf('0.0'),
#                 'CO_g': mpf('0.0'),
#                 'CO_s': mpf('-0.75800000000000000710542735760100185871124267578125'),
#                 'O2_g': mpf('3.50800000000000000710542735760100185871124267578125'),
#                 'O_s': mpf('0.4339999999999999413802242997917346656322479248046875')}
#        ret_e = solver._G
#
#        self.assertDictEqual(ref_e, ret_e)
#
#        # Correction.
#        solver.correct_absolute_energies()
#        ref_e = {'*_s': mpf('0.0'),
#                 'CO-O_2s': mpf('0.9259999999999999342747969421907328069210052490234375'),
#                 'CO2_g': mpf('-0.89590348800350272373549387339153327047824859619140625'),
#                 'CO_g': mpf('-0.84190990132027165859796014046878553926944732666015625'),
#                 'CO_s': mpf('-0.75800000000000000710542735760100185871124267578125'),
#                 'O2_g': mpf('2.63172883483599606702085793585865758359432220458984375'),
#                 'O_s': mpf('0.4339999999999999413802242997917346656322479248046875')}
#        ret_e = solver._G
#
#        self.assertDictEqual(ref_e, ret_e)
#        self.assertTrue(solver.absolute_corrected)

    def test_relative_energies_correction(self):
        " Test solver can correct its relative energies with help of corrector. "
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        parser = model.parser
        parser.parse_data(filename=mkm_energy)
        solver = model.solver
        solver.get_data()

        ret_energies = model.relative_energies

        model.corrector.correct_relative_energies(ret_energies)

        ref_energies = {'Gaf': [0.083909901320271651, 0.0, 1.25],
                        'Gar': [0.0, 1.7637288348359963, 1.8219034880035028],
                        'dG': [0.083909901320271651, -1.7637288348359963, -0.57190348800350277]}

        self.assertDictEqual(ref_energies, ret_energies)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(ThermodynamicCorrectorTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


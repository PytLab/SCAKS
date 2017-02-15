import logging
import unittest
import os

from kynetix.models.micro_kinetic_model import MicroKineticModel
from kynetix.parsers import *

from unit_test import *


class MicroKineticModelTest(unittest.TestCase):

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
        )

    def test_mkm_construction_query(self):
        " Test micro kinetic model can be constructed with parser. "
        # Test construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)

        # Load data in setup file.
        self.assertEqual(model.corrector.__class__.__name__, self.setup_dict["corrector"])
        self.assertEqual(model.plotter.__class__.__name__, self.setup_dict["plotter"])
        self.assertListEqual(model.rxn_expressions, self.setup_dict["rxn_expressions"])
        self.assertEqual(model.temperature, self.setup_dict["temperature"])
        self.assertEqual(model.logger_level, logging.WARNING)
        self.assertEqual(model.decimal_precision, 100)

        self.assertTrue(isinstance(model.parser, RelativeEnergyParser))

    def test_generate_relative_energies_file(self):
        " Test we can generate relative energies input file correctly. "
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)

        abs_path = os.getcwd()
        filename = "{}/rel_energy.py".format(abs_path)
        model.generate_relative_energies_file(filename)

        ref_content = ("# Relative Energies for all elementary reactions.\n" +
                       "Ga, dG = [], []\n\n" +
                       "# CO_g + *_s -> CO_s\n" +
                       "Ga.append()\ndG.append()\n\n" +
                       "# O2_g + 2*_s -> 2O_s\n" +
                       "Ga.append()\ndG.append()\n\n" +
                       "# CO_s + O_s <-> CO-O_2s -> CO2_g + 2*_s\n" +
                       "Ga.append()\ndG.append()\n\n")

        with open(filename, "r") as f:
            ret_content = f.read()
        self.assertEqual(ref_content, ret_content)

        os.remove(filename)

    def test_generate_absolute_energies_file(self):
        " Test we can generate absolute energies input file correctly. "
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)

        abs_path = os.getcwd()
        filename = "{}/abs_energy.py".format(abs_path)
        model.generate_absolute_energies_file(filename)

        ref_content = ("# Absolute energies for all species.\n" +
                       "absolute_energies = {\n\n" +
                       "    'CO2_g': 0.0, # eV\n\n" +
                       "    'CO_g': 0.0, # eV\n\n" +
                       "    'O2_g': 0.0, # eV\n\n" +
                       "    'CO_s': 0.0, # eV\n\n" +
                       "    'O_s': 0.0, # eV\n\n" +
                       "    'CO-O_2s': 0.0, # eV\n\n" +
                       "    '*_s': 0.0, # eV\n\n" +
                       "}\n\n")

        with open(filename, "r") as f:
            ret_content = f.read()
        self.assertEqual(ref_content, ret_content)

        os.remove(filename)

    def test_run(self):
        " Test micro kinetic model can run correctly. "
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        init_cvgs = [0.9, 0.1]
        model.parser.parse_data(filename=mkm_energy)
        model.solver.get_data()
        model.run(init_cvgs=init_cvgs)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(MicroKineticModelTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


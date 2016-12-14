import commands
import logging
import os
import unittest

from kynetix.models.micro_kinetic_model import MicroKineticModel
from kynetix.functions import *
from kynetix.table_makers import *

from unit_test import *


class CsvMakerTest(unittest.TestCase):

    def setUp(self):
        self.maxDiff = None
        self.setup_dict = dict(
            rxn_expressions = [
                'HCOOH_g + 2*_s <-> HCOO-H_s + *_s -> HCOO_s + H_s',
                'HCOO_s + *_s <-> H-COO_s + *_s -> COO_s + H_s',
                'COO_s -> CO2_g + *_s',
                '2H_s <-> H-H_s + *_s -> H2_g + 2*_s',
            ],

            species_definitions = {
                'HCOOH_g': {'pressure': 0.06},
                'H2_g': {'pressure': 0.03},
                'CO2_g': {'pressure': 0.03},
                's': {'site_name': '111', 'type': 'site', 'total': 1.0},
            },
            ref_energies = {
                'C': -8.218910000000001,
                'H': -3.379537,                         
                'O': -7.460926000000001,
                's': -177.41,
            },

            temperature = 450.0,
            table_maker = "CsvMaker",
            parser = "CsvParser",
            solver = "SteadyStateSolver",
        )

    def test_construction(self):
        " Test table maker can be constructed correctly. "
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        table_maker = model.table_maker

        # Check.
        self.assertTrue(table_maker, CsvMaker)

    def test_init_table(self):
        " Test we can initialize table correctly. "
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        table_maker = model.table_maker

        # Init a table.
        table_maker.init_table()

    def test_get_formation_energy(self):
        " Test private function __get_formation_energy(). "
        # Construction.
        model = MicroKineticModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        table_maker = model.table_maker

        # Check.
        ref_e = 0.0
        ret_e = table_maker._CsvMaker__get_formation_energy("H2_g", -6.759074)
        self.assertEqual(ref_e, ret_e)

        ref_e = 0.039835999999979776
        ret_e = table_maker._CsvMaker__get_formation_energy("HCOO-H_s", -207.27)
        self.assertEqual(ref_e, ret_e)

        ref_e = 0.0
        ret_e = table_maker._CsvMaker__get_formation_energy("s", -177.41)
        self.assertEqual(ref_e, ret_e)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(CsvMakerTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


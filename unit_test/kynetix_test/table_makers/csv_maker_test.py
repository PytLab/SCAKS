import commands
import logging
import os
import unittest

from kynetix.model import KineticModel
from kynetix.functions import *
from kynetix.table_makers import *


class CsvMakerTest(unittest.TestCase):

    def setUp(self):
        self.maxDiff = None

    def test_construction(self):
        " Test table maker can be constructed correctly. "
        model = KineticModel(setup_file="input_files/csv_maker.mkm",
                             verbosity=logging.WARNING)
        table_maker = model.table_maker()

        # Check.
        self.assertTrue(table_maker, CsvMaker)

    def test_init_table(self):
        " Test we can initialize table correctly. "
        model = KineticModel(setup_file="input_files/csv_maker.mkm",
                             verbosity=logging.WARNING)
        table_maker = model.table_maker()

        # Init a table.
        table_maker.init_table()

    def test_get_formation_energy(self):
        " Test private function __get_formation_energy(). "
        # Construction.
        model = KineticModel(setup_file="input_files/csv_maker.mkm",
                             verbosity=logging.WARNING)
        table_maker = model.table_maker()

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


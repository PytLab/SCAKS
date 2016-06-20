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

    def test_update_table(self):
        " Make sure we can update table correctly. "
        # Construction.
        model = KineticModel(setup_file="input_files/csv_maker.mkm",
                             verbosity=logging.WARNING)
        table_maker = model.table_maker()
        table_maker.update_table(infile="input_files/in_energy.csv",
                                 outfile="input_files/out_energy.csv",
                                 remove=False)

        ref_content = """species_type,species_name,DFT_energy,formation_energy,frequencies,information
gas,CO2_g,-22.975866,0.16489600000000237,[],None
gas,H2_g,-6.759074,0.0,[],None
gas,HCOOH_g,-29.899836,0.0,[],None
intermediate,COO_s,-200.23,0.320762000000002,[],None
intermediate,HCOO_s,-204.00,-0.06970100000000912,[],None
intermediate,H_s,-181.36,-0.5704630000000179,[],None
transition state,H-COO_s,-203.18,0.7502989999999841,[],None
transition state,H-H_s,-184.35,-0.18092599999999948,[],None
transition state,HCOO-H_s,-207.27,0.039835999999979776,[],None
slab,s,-177.41,0.0,[],None
"""
        with open("input_files/out_energy.csv") as f:
            ret_content = f.read()

        self.assertEqual(ref_content, ret_content)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(CsvMakerTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


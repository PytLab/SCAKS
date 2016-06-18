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

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(CsvMakerTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


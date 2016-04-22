'''
Module for testing the relative energy parser class.
'''

import unittest

from pynetics import model


class MicroKineticsModelTest(unittest.TestCase):

    def setUp(self):
        #create an instance of KineticModel
        self.km = model.KineticModel(setup_file='mkm_model_test.mkm',
                                     logging_level=30)
        self.parser = self.km.parser
        self.maxDiff = None


if __name__ == '__main__':
    unittest.main()

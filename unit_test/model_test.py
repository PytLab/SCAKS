import unittest
import sys


class TestKineticModel(unittest.TestCase):

    def setUp(self):
        #create an instance of KineticModel
        sys.path.append('E:\BaiduYun\MyDocuments\ECUST+\Python\script')
        from kinetic import model
        self.km = model.KineticModel(setup_file='methanation.mkm')

    def test_set_logger(self):
        "Make sure set_logger method takes effects"
        self.assertTrue(hasattr(self.km, 'logger'))

        #should raise an exception for an AttributeError
        self.assertRaises(AttributeError)

    def test_load(self):
        "Make sure the load() method has loaded variables in setup file"

        #make sure all vars in setup file are parsed in
        globs, locs = {}, {}
        execfile(self.km.setup_file, globs, locs)
        for var in locs:
            self.assertTrue(hasattr(self.km, var))
            self.assertRaises(AttributeError)
        #make sure all tools are parsed in
        tools_list = ['parser', 'table_maker']
        for tool in tools_list:
            self.assertTrue(hasattr(self.km, tool))
            self.assertRaises(AttributeError)

#if __name__ == '__main__':
#    unittest.main()
suite = unittest.TestLoader().loadTestsFromTestCase(TestKineticModel)
unittest.TextTestRunner(verbosity=2).run(suite)
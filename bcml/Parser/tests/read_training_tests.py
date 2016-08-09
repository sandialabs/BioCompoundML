"""

This contains the unit tests for the read_training module.

"""


import unittest
import os
from Parser import read_training as rt

"""This file is pulled from the data directory"""
_dir = os.path.dirname(__file__)
_input_file = os.path.abspath(os.path.join(_dir, os.pardir, os.pardir, 'data', 'RON_min.txt'))
_count = 11
_name = '1-Butene'
_pubchem_id = '7844'
_RON = '98.8'
_id = 'PubChem'
_pred = 'RON'


class ReadTrainingTests(unittest.TestCase):

    def setUp(self):
        """Create an instance of the read_training Read class"""
        print "Initializing test"
        self.read = rt.Read(_input_file, predictor=_pred, user=False, id_name=_id)

    def tearDown(self):
        """Delete data structure"""
        print "Clearing out file"
        del self.read

    def test_training_compounds(self):
        """Test the values associated with the test file"""
        print "Finding compound"
        self.assertEqual(len(self.read.compounds), _count)
        test_compound = self.read.compounds[0]
        print "Testing name"
        self.assertEqual(test_compound['Name'], _name)
        print "Testing ID"
        self.assertEqual(test_compound['PubChem'], _pubchem_id)
        print "Testing RON"
        self.assertEqual(test_compound['RON'], _RON)
        print "Testing with User"
        del self.read
        self.read = rt.Read(_input_file, predictor=_pred, user=True, id_name=_id)
        test_compound = self.read.compounds[0]
        self.assertTrue('userhash' in test_compound)
        print "Testing with Predictor"
        self.assertFalse(_pred in test_compound['userhash'])
        self.assertEqual(self.read.predictors[0], _RON)

if __name__ == '__main__':
    unittest.main()

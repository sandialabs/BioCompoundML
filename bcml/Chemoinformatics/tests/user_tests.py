"""

This contains the unit tests for the read_training module.

"""


from __future__ import print_function
import unittest
from Chemoinformatics import user
from collections import defaultdict, OrderedDict
from copy import deepcopy


class Object(object):
    pass


class UserTest(unittest.TestCase):

    def setUp(self):
        """Create an instance of the Update class"""
        print("Initializing test")
        compound = OrderedDict({'19502': {'userhash': {'PubChem': '19502', 'Name': '1-Ethyl-3-Methylcyclopentane'}}, '7459': {'userhash': {'PubChem': '7459', 'Name': '1-Isopropyl-4-methylcyclohexane'}}, '8004': {'userhash': {'PubChem': '8004', 'Name': '1-Pentene'}}, '11597': {'userhash': {'PubChem': '11597', 'Name': '1-Hexene'}}, '107252': {'userhash': {'PubChem': '107252', 'Name': '1-Methyl-2-propylcyclohexane'}}, '11610': {'userhash': {'PubChem': '11610', 'Name': '1-Heptene'}}, '35411': {'userhash': {'PubChem': '35411', 'Name': '1-Methyl-1-ethylcyclohexane'}}, '136729': {'userhash': {'PubChem': '136729', 'Name': '1-Methyl-2-ethylcyclopentane'}}, '7844': {'userhash': {'PubChem': '7844', 'Name': '1-Butene'}}, '8125': {'userhash': {'PubChem': '8125', 'Name': '1-Octene'}}, '11549': {'userhash': {'PubChem': '11549', 'Name': '1,1-Dimethylcyclohexane'}}})
        self.test_data = Object()
        self.test_data.compound = deepcopy(compound)
        self.original = Object()
        self.original.compound = deepcopy(compound)

    def tearDown(self):
        """Delete data structure"""
        print("Clearing out file")
        del self.test_data
        del self.original

    def test_fingerprint_functions(self):
        """Test the values associated with the test file"""
        print("Run Fingerprint functions")
        self.test_data = user.Update(self.test_data, remove_static=True, verbose=True)
        self.test_data.update()
        self.assertTrue('7844' in self.test_data.compound)
        self.assertTrue('userhash' in self.test_data.compound['7844'])
        self.assertTrue('PubChem' in self.test_data.compound['7844']['userhash'])
        self.assertEquals('7844', self.test_data.compound['7844']['userhash']['PubChem'])
        variable_length = len(self.test_data.compound['7844']['userhash'].keys())
        self.assertEquals(variable_length, 1)
        self.original = user.Update(self.original, remove_static=False, verbose=True)
        self.original.update()
        total_length = len(self.original.compound['7844']['userhash'].keys())
        self.assertEquals(total_length, 2)

if __name__ == '__main__':
    unittest.main()

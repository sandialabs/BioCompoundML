"""

This contains the unit tests for the read_training module.

"""


from __future__ import print_function
import unittest
from Distance import distance as ds


def dictitems(dict):
    if sys.version_info[0]>=3:
        return dict.items()
    else:
        return dictitems(dict)


class DistanceTests(unittest.TestCase):

    def setUp(self):
        """Create an instance of the read_training Read class"""
        print("Initializing test")
        input_data = {'19502': {'CC1CCC(C)CC1': '0', '[#1]-C=C-[#1]': '0',
                                'C-C-C-C-C-C': '1', 'C-C-C-C-C-C(C)-C': '1',
                                '>= 8 C': '1', 'C-C(C)-C-C-C': '1',
                                'C(~C)(~C)(~C)': '1', '>= 16 H': '1',
                                'C-C-C(C)-C-C': '1', 'CC1CC(C)CC1': '1',
                                'CC1C(C)CCCC1': '0', 'C-C-C=C': '0', 'C-C(C)-C-C': '1',
                                'C(~C)(~C)(~C)(~H)': '1', 'C-C-C-C-C': '1',
                                'C(~C)(~C)(~C)(~C)': '0', 'C-C(C)-C(C)-C': '0',
                                'C=C-C-C-C': '0', 'C=C': '0', 'CC1C(C)CCC1': '0',
                                'C(-C)(-H)(=C)': '0', 'C-C(C)(C)-C-C': '0',
                                '[#1]-C-C=C-[#1]': '0',
                                '>= 1 saturated or aromatic carbon-only ring size 5': '1',
                                '>= 1 saturated or aromatic carbon-only ring size 6': '0',
                                'C=C-C-C': '0', 'C-C-C-C-C-C-C-C': '0', 'C-C-C-C-C-C-C': '1',
                                'C(-H)(=C)': '0', '>= 1 any ring size 5': '1',
                                '>= 1 any ring size 6': '0', 'C(-C)(=C)': '0'},
                      '7459': {'CC1CCC(C)CC1': '1', '[#1]-C=C-[#1]': '0', 'C-C-C-C-C-C': '1',
                               'C-C-C-C-C-C(C)-C': '1', '>= 8 C': '1', 'C-C(C)-C-C-C': '1',
                               'C(~C)(~C)(~C)': '1', '>= 16 H': '1', 'C-C-C(C)-C-C': '1',
                               'CC1CC(C)CC1': '0', 'CC1C(C)CCCC1': '0', 'C-C-C=C': '0',
                               'C-C(C)-C-C': '1', 'C(~C)(~C)(~C)(~H)': '1', 'C-C-C-C-C': '1',
                               'C(~C)(~C)(~C)(~C)': '0', 'C-C(C)-C(C)-C': '1', 'C=C-C-C-C': '0',
                               'C=C': '0', 'CC1C(C)CCC1': '0', 'C(-C)(-H)(=C)': '0', 'C-C(C)(C)-C-C': '0',
                               '[#1]-C-C=C-[#1]': '0',
                               '>= 1 saturated or aromatic carbon-only ring size 5': '0',
                               '>= 1 saturated or aromatic carbon-only ring size 6': '1',
                               'C=C-C-C': '0', 'C-C-C-C-C-C-C-C': '1', 'C-C-C-C-C-C-C': '1',
                               'C(-H)(=C)': '0', '>= 1 any ring size 5': '0',
                               '>= 1 any ring size 6': '1', 'C(-C)(=C)': '0'},
                      '136729': {'CC1CCC(C)CC1': '0', '[#1]-C=C-[#1]': '0', 'C-C-C-C-C-C': '1',
                                 'C-C-C-C-C-C(C)-C': '0', '>= 8 C': '1', 'C-C(C)-C-C-C': '1',
                                 'C(~C)(~C)(~C)': '1', '>= 16 H': '1', 'C-C-C(C)-C-C': '1',
                                 'CC1CC(C)CC1': '0', 'CC1C(C)CCCC1': '0', 'C-C-C=C': '0',
                                 'C-C(C)-C-C': '1', 'C(~C)(~C)(~C)(~H)': '1', 'C-C-C-C-C': '1',
                                 'C(~C)(~C)(~C)(~C)': '0', 'C-C(C)-C(C)-C': '1', 'C=C-C-C-C': '0',
                                 'C=C': '0', 'CC1C(C)CCC1': '1', 'C(-C)(-H)(=C)': '0','C-C(C)(C)-C-C': '0',
                                 '[#1]-C-C=C-[#1]': '0',
                                 '>= 1 saturated or aromatic carbon-only ring size 5': '1',
                                 '>= 1 saturated or aromatic carbon-only ring size 6': '0',
                                 'C=C-C-C': '0', 'C-C-C-C-C-C-C-C': '1', 'C-C-C-C-C-C-C': '1',
                                 'C(-H)(=C)': '0', '>= 1 any ring size 5': '1',
                                 '>= 1 any ring size 6': '0', 'C(-C)(=C)': '0'},
                      '11597': {'CC1CCC(C)CC1': '0', '[#1]-C=C-[#1]': '1', 'C-C-C-C-C-C': '0',
                                'C-C-C-C-C-C(C)-C': '0', '>= 8 C': '0', 'C-C(C)-C-C-C': '0',
                                'C(~C)(~C)(~C)': '0', '>= 16 H': '0', 'C-C-C(C)-C-C': '0',
                                'CC1CC(C)CC1': '0', 'CC1C(C)CCCC1': '0', 'C-C-C=C': '1',
                                'C-C(C)-C-C': '0', 'C(~C)(~C)(~C)(~H)': '0', 'C-C-C-C-C': '1',
                                'C(~C)(~C)(~C)(~C)': '0', 'C-C(C)-C(C)-C': '0', 'C=C-C-C-C': '1',
                                'C=C': '1', 'CC1C(C)CCC1': '0', 'C(-C)(-H)(=C)': '1', 'C-C(C)(C)-C-C': '0',
                                '[#1]-C-C=C-[#1]': '1',
                                '>= 1 saturated or aromatic carbon-only ring size 5': '0',
                                '>= 1 saturated or aromatic carbon-only ring size 6': '0',
                                'C=C-C-C': '1', 'C-C-C-C-C-C-C-C': '0', 'C-C-C-C-C-C-C': '0',
                                'C(-H)(=C)': '1', '>= 1 any ring size 5': '0',
                                '>= 1 any ring size 6': '0', 'C(-C)(=C)': '1'}}
        self.fingerprint_vector = list()
        self.key_list = list()
        for (key, value) in dictitems(input_data):
            self.fingerprint_vector.append(value)
            self.key_list.append(key)

    def tearDown(self):
        """Delete data structure"""
        print("Clearing out file")
        del self.fingerprint_vector

    def test_distance_functions(self):
        """Test the values associated with the test file"""
        print("Run standard distance")
        distance = ds.Distance(self.fingerprint_vector, self.key_list)
        """Test jaccard distance between 0 and 1"""
        self.assertAlmostEqual(0.83636364, distance.distance[0, 1])
        (m, n) = distance.distance.shape
        self.assertEqual(4, m)
        self.assertEqual(4, n)
        distance = ds.Distance(self.fingerprint_vector, self.key_list, method='cosine')
        self.assertAlmostEqual(7.18750000e-01, distance.distance[0, 1])
        distance = ds.Distance(self.fingerprint_vector, self.key_list, n_jobs=1, method='euclidean')
        self.assertAlmostEqual(6.78232998, distance.distance[0, 1])

if __name__ == '__main__':
    unittest.main()

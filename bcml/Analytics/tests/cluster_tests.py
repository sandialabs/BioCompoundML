"""

This contains the unit tests for the cluster module.

"""


from __future__ import print_function
import unittest
from Analytics import cluster as cl
from Parser import build_training as bt
from Train import train_model as tm
from copy import deepcopy
import numpy as np
import sys
from sklearn.preprocessing import Imputer
from Distance import distance as ds


_split_value = 70
_random = 12345
_n_iter = 10


def dictitems(dict):
    if sys.version_info[0]>=3:
        return dict.items()
    else:
        return dict.iteritems()


class Object(object):
    pass


class ClusterTests(unittest.TestCase):

    def setUp(self):
        """Create an instance of the read_training Read class"""
        print("Initializing test")
        input_data = {'7844': {'CC1CCC(C)CC1': '0', '[#1]-C=C-[#1]': '0',
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
                      '19502': {'CC1CCC(C)CC1': '1', '[#1]-C=C-[#1]': '0', 'C-C-C-C-C-C': '1',
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
                      '11610': {'CC1CCC(C)CC1': '0', '[#1]-C=C-[#1]': '0', 'C-C-C-C-C-C': '1',
                                 'C-C-C-C-C-C(C)-C': '0', '>= 8 C': '1', 'C-C(C)-C-C-C': '1',
                                 'C(~C)(~C)(~C)': '1', '>= 16 H': '1', 'C-C-C(C)-C-C': '1',
                                 'CC1CC(C)CC1': '0', 'CC1C(C)CCCC1': '0', 'C-C-C=C': '0',
                                 'C-C(C)-C-C': '1', 'C(~C)(~C)(~C)(~H)': '1', 'C-C-C-C-C': '1',
                                 'C(~C)(~C)(~C)(~C)': '0', 'C-C(C)-C(C)-C': '1', 'C=C-C-C-C': '0',
                                 'C=C': '0', 'CC1C(C)CCC1': '1', 'C(-C)(-H)(=C)': '0', 'C-C(C)(C)-C-C': '0',
                                 '[#1]-C-C=C-[#1]': '0',
                                 '>= 1 saturated or aromatic carbon-only ring size 5': '1',
                                 '>= 1 saturated or aromatic carbon-only ring size 6': '0',
                                 'C=C-C-C': '0', 'C-C-C-C-C-C-C-C': '1', 'C-C-C-C-C-C-C': '1',
                                 'C(-H)(=C)': '0', '>= 1 any ring size 5': '1',
                                 '>= 1 any ring size 6': '0', 'C(-C)(=C)': '0'}}
        self.fingerprint_vector = list()
        self.key_list = list()
        for (key, value) in dictitems(input_data):
            self.fingerprint_vector.append(value)
            self.key_list.append(key)
        self.distance = ds.Distance(self.fingerprint_vector, self.key_list).distance
        compound = {'7844': {'predictor': '98.8', 'experimentalhash': {u'Density': 0.577, u'Vapor Density': 1.93, u'Boiling Point': -6.47, u'Rotatable Bond Count': 1.0, u'XLogP3': 2.4, u'Melting Point': -185.3, u'Flash Point': False, u'Undefined Atom Stereocenter Count': 0.0, u'Auto-Ignition': 725.0, u'Molecular Weight': 56.10632, u'LogP': 2.4, u'Complexity': 14.0, u'Vapor Pressure': 2253.0, u'Heavy Atom Count': 4.0, u'Exact Mass': 56.0626, u'Monoisotopic Mass': 56.0626}}, '19502': {'predictor': '57.6', 'experimentalhash': {u'Rotatable Bond Count': 1.0, u'Heavy Atom Count': 8.0, u'Undefined Atom Stereocenter Count': 3.0, u'Molecular Weight': 112.21264, u'Complexity': 66.4, u'Exact Mass': 112.125201, u'Monoisotopic Mass': 112.125201}}, '11610': {'predictor': '54.5', 'experimentalhash': {u'Density': 0.697, u'Vapor Density': 0.7, u'Boiling Point': 93.6, u'Rotatable Bond Count': 4.0, u'XLogP3': 4.0, u'Melting Point': -119.7, u'Flash Point': 32.0, u'Undefined Atom Stereocenter Count': 0.0, u'Auto-Ignition': 500.0, u'Molecular Weight': 98.18606, u'LogP': 3.99, u'Complexity': 37.3, u'Vapor Pressure': 59.3, u'Heavy Atom Count': 7.0, u'Exact Mass': 98.10955, u'Monoisotopic Mass': 98.10955}}}
        self.test_data = Object()
        self.test_data.compound = deepcopy(compound)
        np.random.seed(_random)
        train = bt.Process(self.test_data, split_value=_split_value)
        X = train.train
        imp = Imputer(missing_values='NaN', strategy='mean', axis=0)
        imp.fit(X)
        t = imp.transform(X)
        train.train = t
        self.model = tm.Train(train)
        self.model.train_model()

    def tearDown(self):
        """Delete data structure"""
        print("Clearing out file")
        del self.test_data

    def testCluster(self):
        np.random.seed(_random)
        print("Testing Clustering")
        cluster = cl.Clustering(self.test_data.compound, seed=_random)
        cluster.cluster_training(self.model)
        self.assertEqual(0.9375, cluster.p_feat_matrix[0][1])
        cluster.cluster_training(self.model, distance=self.distance)
        self.assertEqual(0.4, cluster.p_feat_matrix[0][1])


if __name__ == '__main__':
    unittest.main()

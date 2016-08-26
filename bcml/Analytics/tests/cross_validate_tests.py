"""

This contains the unit tests for the cross_validate module.

"""


from __future__ import print_function
import unittest
from Analytics import cross_validate as cv
from Parser import build_training as bt
from Train import train_model as tm
from copy import deepcopy
import numpy as np
from sklearn.preprocessing import Imputer

_split_value = 70
_random = 12345
_n_iter = 10


class Object(object):
    pass


class CrossValidateTests(unittest.TestCase):

    def setUp(self):
        """Create an instance of the read_training Read class"""
        print("Initializing test")
        compound = {'7844': {'predictor': '98.8', 'experimentalhash': {u'Density': 0.577, u'Vapor Density': 1.93, u'Boiling Point': -6.47, u'Rotatable Bond Count': 1.0, u'XLogP3': 2.4, u'Melting Point': -185.3, u'Flash Point': False, u'Undefined Atom Stereocenter Count': 0.0, u'Auto-Ignition': 725.0, u'Molecular Weight': 56.10632, u'LogP': 2.4, u'Complexity': 14.0, u'Vapor Pressure': 2253.0, u'Heavy Atom Count': 4.0, u'Exact Mass': 56.0626, u'Monoisotopic Mass': 56.0626}}, '19502': {'predictor': '57.6', 'experimentalhash': {u'Rotatable Bond Count': 1.0, u'Heavy Atom Count': 8.0, u'Undefined Atom Stereocenter Count': 3.0, u'Molecular Weight': 112.21264, u'Complexity': 66.4, u'Exact Mass': 112.125201, u'Monoisotopic Mass': 112.125201}}, '11610': {'predictor': '54.5', 'experimentalhash': {u'Density': 0.697, u'Vapor Density': 0.7, u'Boiling Point': 93.6, u'Rotatable Bond Count': 4.0, u'XLogP3': 4.0, u'Melting Point': -119.7, u'Flash Point': 32.0, u'Undefined Atom Stereocenter Count': 0.0, u'Auto-Ignition': 500.0, u'Molecular Weight': 98.18606, u'LogP': 3.99, u'Complexity': 37.3, u'Vapor Pressure': 59.3, u'Heavy Atom Count': 7.0, u'Exact Mass': 98.10955, u'Monoisotopic Mass': 98.10955}}, '7855': {'predictor': '98.8', 'experimentalhash': {u'Density': 0.577, u'Vapor Density': 1.93, u'Boiling Point': -6.47, u'Rotatable Bond Count': 1.0, u'XLogP3': 2.4, u'Melting Point': -185.3, u'Flash Point': False, u'Undefined Atom Stereocenter Count': 0.0, u'Auto-Ignition': 725.0, u'Molecular Weight': 56.10632, u'LogP': 2.4, u'Complexity': 14.0, u'Vapor Pressure': 2253.0, u'Heavy Atom Count': 4.0, u'Exact Mass': 56.0626, u'Monoisotopic Mass': 56.0626}}, '19503': {'predictor': '57.6', 'experimentalhash': {u'Rotatable Bond Count': 1.0, u'Heavy Atom Count': 8.0, u'Undefined Atom Stereocenter Count': 3.0, u'Molecular Weight': 112.21264, u'Complexity': 66.4, u'Exact Mass': 112.125201, u'Monoisotopic Mass': 112.125201}}, '11611': {'predictor': '54.5', 'experimentalhash': {u'Density': 0.697, u'Vapor Density': 0.7, u'Boiling Point': 93.6, u'Rotatable Bond Count': 4.0, u'XLogP3': 4.0, u'Melting Point': -119.7, u'Flash Point': 32.0, u'Undefined Atom Stereocenter Count': 0.0, u'Auto-Ignition': 500.0, u'Molecular Weight': 98.18606, u'LogP': 3.99, u'Complexity': 37.3, u'Vapor Pressure': 59.3, u'Heavy Atom Count': 7.0, u'Exact Mass': 98.10955, u'Monoisotopic Mass': 98.10955}}, '7864': {'predictor': '98.8', 'experimentalhash': {u'Density': 0.577, u'Vapor Density': 1.93, u'Boiling Point': -6.47, u'Rotatable Bond Count': 1.0, u'XLogP3': 2.4, u'Melting Point': -185.3, u'Flash Point': False, u'Undefined Atom Stereocenter Count': 0.0, u'Auto-Ignition': 725.0, u'Molecular Weight': 56.10632, u'LogP': 2.4, u'Complexity': 14.0, u'Vapor Pressure': 2253.0, u'Heavy Atom Count': 4.0, u'Exact Mass': 56.0626, u'Monoisotopic Mass': 56.0626}}}
        self.test_data = Object()
        self.test_data.compound = deepcopy(compound)

    def tearDown(self):
        """Delete data structure"""
        print("Clearing out file")
        del self.test_data

    def testCrossValidate(self):
        np.random.seed(seed=_random)
        train = bt.Process(self.test_data, split_value=_split_value)
        X = train.train
        imp = Imputer(missing_values='NaN', strategy='mean', axis=0)
        imp.fit(X)
        t = imp.transform(X)
        train.train = t
        model = tm.Train(train)
        model.train_model()
        print("Test cross_validate")
        cross = cv.Analysis(model, seed=_random, verbose=True)
        cross.cross_validate(n_iter=_n_iter)
        self.assertEqual(0.95, cross.acc)
        self.assertEqual(0.95, cross.prec)
        self.assertEqual(1.0, cross.recall)
        self.assertEqual(0.95, cross.roc)
        cross.feature_importances()
        self.assertAlmostEqual(0.0586, cross.feats[0][0], 3)

if __name__ == '__main__':
    unittest.main()

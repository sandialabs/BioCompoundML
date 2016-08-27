"""

This contains the unit tests for the build_training module.

"""


from __future__ import print_function
import unittest
from Parser import build_training as bt
from copy import deepcopy
import numpy as np
try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict

"""This file is pulled from the data directory"""
_split_value = 70
_distance = np.array([[0., 0.04872647, 0.00678733], [0.04872647, 0., 0.04656319], [0.00678733, 0.04656319, 0.]])
_verbose = True

def list_check(list1, list2):
	for i, val in enumerate(list1):
		if list1[i] != list2[i]:
			return False
	return True


class Object(object):
    pass


class BuildTrainingTests(unittest.TestCase):

    def setUp(self):
        """Create an instance of the read_training Read class"""
        print("Initializing test")
        compound = {'7844': {'predictor': '98.8', 'experimentalhash': {u'Density': 0.577, u'Vapor Density': 1.93, u'Boiling Point': -6.47, u'Rotatable Bond Count': 1.0, u'XLogP3': 2.4, u'Melting Point': -185.3, u'Flash Point': False, u'Undefined Atom Stereocenter Count': 0.0, u'Auto-Ignition': 725.0, u'Molecular Weight': 56.10632, u'LogP': 2.4, u'Complexity': 14.0, u'Vapor Pressure': 2253.0, u'Heavy Atom Count': 4.0, u'Exact Mass': 56.0626, u'Monoisotopic Mass': 56.0626}}, '19502': {'predictor': '57.6', 'experimentalhash': {u'Rotatable Bond Count': 1.0, u'Heavy Atom Count': 8.0, u'Undefined Atom Stereocenter Count': 3.0, u'Molecular Weight': 112.21264, u'Complexity': 66.4, u'Exact Mass': 112.125201, u'Monoisotopic Mass': 112.125201}}, '11610': {'predictor': '54.5', 'experimentalhash': {u'Density': 0.697, u'Vapor Density': 0.7, u'Boiling Point': 93.6, u'Rotatable Bond Count': 4.0, u'XLogP3': 4.0, u'Melting Point': -119.7, u'Flash Point': 32.0, u'Undefined Atom Stereocenter Count': 0.0, u'Auto-Ignition': 500.0, u'Molecular Weight': 98.18606, u'LogP': 3.99, u'Complexity': 37.3, u'Vapor Pressure': 59.3, u'Heavy Atom Count': 7.0, u'Exact Mass': 98.10955, u'Monoisotopic Mass': 98.10955}}}
        compound = OrderedDict(sorted(compound.items(), key=lambda t: t[0]))
        self.test_data = Object()
        self.test_data.compound = deepcopy(compound)

    def tearDown(self):
        """Delete data structure"""
        print("Clearing out file")
        del self.test_data

    def testBuildTraining(self):
        train = bt.Process(self.test_data, split_value=_split_value)
        distance = _distance
        a = train.train[0][0] * (1./distance[0][1])
        b = train.train[2][0] * (1./distance[2][1])
        denominator = (1./distance[0][1]) + (1./distance[2][1])
        value = (a+b)/denominator
        print("Test imputing")
        train.impute_values(distance, verbose=_verbose, k=2)
        self.assertAlmostEquals(value, train.train[1][0])
        print("Testing Boruta")
        compound = {'19502': {'predictor': '57.6', 'experimentalhash': {u'Rotatable Bond Count': 1.0, u'XLogP3-AA': 3.8, u'Heavy Atom Count': 8.0, u'Undefined Atom Stereocenter Count': 3.0, u'Molecular Weight': 112.21264, u'Complexity': 66.4, u'Exact Mass': 112.125201, u'Monoisotopic Mass': 112.125201}}, '7459': {'predictor': '67.3', 'experimentalhash': {u'Density': 0.8039, u'Boiling Point': 170.9, u'Rotatable Bond Count': 1.0, u'Melting Point': -89.84, u'Undefined Atom Stereocenter Count': 0.0, u'Molecular Weight': 140.2658, u'XLogP3-AA': 4.5, u'Complexity': 86.2, u'Heavy Atom Count': 10.0, u'Exact Mass': 140.156501, u'Monoisotopic Mass': 140.156501}}, '8004': {'predictor': '87.9', 'experimentalhash': {u'Density': 0.6405, u'Vapor Density': 2.42, u'Boiling Point': 29.9, u'Rotatable Bond Count': 2.0, u'Melting Point': -165.2, u'Flash Point': -18.0, u'Undefined Atom Stereocenter Count': 0.0, u'Auto-Ignition': 527.0, u'Molecular Weight': 70.1329, u'XLogP3-AA': 2.4, u'Complexity': 21.2, u'Vapor Pressure': 635.0, u'Heavy Atom Count': 5.0, u'Exact Mass': 70.07825, u'Monoisotopic Mass': 70.07825}}, '11597': {'predictor': '76.4', 'experimentalhash': {u'Density': 0.6731, u'Vapor Density': 3.0, u'Boiling Point': 63.4, u'Rotatable Bond Count': 3.0, u'XLogP3': 3.4, u'Melting Point': -139.7, u'Flash Point': 20.0, u'Undefined Atom Stereocenter Count': 0.0, u'Auto-Ignition': 487.0, u'Molecular Weight': 84.15948, u'LogP': 3.39, u'Complexity': 29.0, u'Vapor Pressure': 183.7, u'Heavy Atom Count': 6.0, u'Exact Mass': 84.0939, u'Monoisotopic Mass': 84.0939}}, '107252': {'predictor': '29.9', 'experimentalhash': {u'XLogP3-AA': 4.8, u'Exact Mass': 140.156501, u'Monoisotopic Mass': 140.156501, u'Undefined Atom Stereocenter Count': 2.0, u'Complexity': 86.0, u'Rotatable Bond Count': 2.0, u'Molecular Weight': 140.2658, u'Heavy Atom Count': 10.0}}, '11610': {'predictor': '54.5', 'experimentalhash': {u'Density': 0.697, u'Vapor Density': 0.7, u'Boiling Point': 93.6, u'Rotatable Bond Count': 4.0, u'XLogP3': 4.0, u'Melting Point': -119.7, u'Flash Point': 32.0, u'Undefined Atom Stereocenter Count': 0.0, u'Auto-Ignition': 500.0, u'Molecular Weight': 98.18606, u'LogP': 3.99, u'Complexity': 37.3, u'Vapor Pressure': 59.3, u'Heavy Atom Count': 7.0, u'Exact Mass': 98.10955, u'Monoisotopic Mass': 98.10955}}, '35411': {'predictor': '68.7', 'experimentalhash': {u'XLogP3-AA': 4.4, u'Exact Mass': 126.140851, u'Monoisotopic Mass': 126.140851, u'Undefined Atom Stereocenter Count': 0.0, u'Complexity': 78.0, u'Rotatable Bond Count': 1.0, u'Molecular Weight': 126.23922, u'Heavy Atom Count': 9.0}}, '136729': {'predictor': '57.6', 'experimentalhash': {u'Rotatable Bond Count': 1.0, u'XLogP3-AA': 3.8, u'Heavy Atom Count': 8.0, u'Undefined Atom Stereocenter Count': 2.0, u'Molecular Weight': 112.21264, u'Complexity': 66.4, u'Exact Mass': 112.125201, u'Monoisotopic Mass': 112.125201}}, '7844': {'predictor': '98.8', 'experimentalhash': {u'Density': 0.577, u'Vapor Density': 1.93, u'Boiling Point': -6.47, u'Rotatable Bond Count': 1.0, u'XLogP3': 2.4, u'Melting Point': -185.3, u'Flash Point': False, u'Undefined Atom Stereocenter Count': 0.0, u'Auto-Ignition': 725.0, u'Molecular Weight': 56.10632, u'LogP': 2.4, u'Complexity': 14.0, u'Vapor Pressure': 2253.0, u'Heavy Atom Count': 4.0, u'Exact Mass': 56.0626, u'Monoisotopic Mass': 56.0626}}, '8125': {'predictor': '28.7', 'experimentalhash': {u'Density': 0.7149, u'Vapor Density': 3.87, u'Boiling Point': 121.2, u'Rotatable Bond Count': 5.0, u'XLogP3': 4.6, u'Melting Point': -101.7, u'Flash Point': 70.0, u'Undefined Atom Stereocenter Count': 0.0, u'Auto-Ignition': 446.0, u'Molecular Weight': 112.21264, u'LogP': 4.57, u'Complexity': 46.0, u'Vapor Pressure': 17.4, u'Heavy Atom Count': 8.0, u'Exact Mass': 112.125201, u'Monoisotopic Mass': 112.125201}}, '11549': {'predictor': '87.3', 'experimentalhash': {u'XLogP3-AA': 3.9, u'Exact Mass': 112.125201, u'Monoisotopic Mass': 112.125201, u'Undefined Atom Stereocenter Count': 0.0, u'Complexity': 64.1, u'Rotatable Bond Count': 0.0, u'Molecular Weight': 112.21264, u'Heavy Atom Count': 8.0}}}
        test_data = Object()
        test_data.compound = deepcopy(compound)
        train = bt.Process(test_data, split_value=_split_value)
        train.feature_selection(seed=12345, verbose=False)
        support = [False, False, True, False, False, False, False, False, False, False, False, False, False, False, False, False, False]
        self.assertTrue(list_check(support, train.feature_support.tolist()))
        self.assertEquals('Complexity', train.feature_names[0])
        '''
        Need an appropriate training set to test median split
        test_data = Object()
        test_data.compound = deepcopy(compound)
        train = bt.Process(test_data, split_value=False)
        train.feature_selection(seed=12345, verbose=True)
        print(train.feature_support)
        support = [False, False, False, False, False, False, False, False, False, False, True, False, False, False, False, False, False]
        self.assertTrue(list_check(support, train.feature_support.tolist()))
        self.assertEquals('Complexity', train.feature_names[0])
        '''

if __name__ == '__main__':
    unittest.main()

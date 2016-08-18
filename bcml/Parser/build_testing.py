"""
This process takes in the training dataset and outputs
a data structure that includes the name of the molecule,
the predictor, and the CAS number
Attributes:
    input_file (str): This is the training file that
    is read by the output
    Instance (class): This is a private class which
    structures each instance
    Model (class): This is a public class with the
    total structure of the set
"""

from __future__ import print_function
import numpy as np
from KNNImpute.knnimpute import (
    knn_impute_optimistic,
)
from sklearn.preprocessing import Imputer

#_possible_features = ("binhash", "padelhash", "experimentalhash")
_possible_features = ('experimentalhash', 'binhash', 'padelhash', 'userhash')


def dictitems(dict):
    if sys.version_info[0]>=3:
        return dict.items()
    else:
        return dictitems(dict)


def verbose_print(verbose, line):
    if verbose:
        print(line)


def _get_feature_names(compounds):
    """This function handles collecting the feature names"""
    feature_names = {}
    for compound in compounds:
        for feature in _possible_features:
            if feature in compound.keys():
                keys = compound[feature].keys()
                for feat in keys:
                    feature_names[feat] = 1
                    compound[feat] = compound[feature][feat]
    return (compounds, feature_names.keys())


class Process(object):
    """This file reads a training file"""
    def load_testing_set(self):
        """This function takes the features and
        compounds and loads them into a numpy array
        """
        for index, value in np.ndenumerate(self.test):
            compound = self.compounds[index[0]]
            feature = self.features[index[1]]
            if feature in compound.keys() and compound[feature] is not "" and compound[feature] != "NULL":
                self.test[index] = float(compound[feature])
            else:
                self.test[index] = np.nan

    def impute_values(self, distance=False, k=5, verbose=True):
        """This function handles the missing values from
        the training set and estimates their value, based on
        the mean and reloads them into the training set"""
        verbose_print(verbose, 'Imputing using KNN strategy')
        X = self.test
        missing_mask = np.isnan(X)
        '''First impute uing knn optimistic'''
        impute = knn_impute_optimistic(X, missing_mask=missing_mask,
                                       distance=distance, k=k)
        X = impute.astype(float)
        '''For features with a small number of features, use mean
        imputation to remove NaN values'''
        imp = Imputer(missing_values='NaN', strategy='mean', axis=0)
        imp.fit(X)
        self.test = imp.transform(X)

    '''def impute_values(self, distance=False):
        """This function handles the missing values from
        the training set and estimates their value, based on
        the mean and reloads them into the training set
        A smarter way of handling this is to impute based on the
        distance of each value from the combined training and
        test sets.
        """
        for index, value in np.ndenumerate(self.test):
            column = index[1]
            if (np.isnan(value) is True) or (np.isinf(value) is True) or (np.can_cast(value, np.float64) is False) or (np.can_cast(value, np.float32) is False):
                self.test[index] = self.impute[column]'''

    def __init__(self, testing_data, features, training_data):
        """This initialization function handles the heavy
        work of loading the features and processing the
        compounds"""
        self.input = testing_data
        self.features = features
        self.impute = training_data.mean(axis=0)
        self.columns = len(self.features)
        self.rows = len(testing_data.compound)
        self.test = np.zeros((self.rows, self.columns,), dtype=np.float64)
        compounds = []
        self.test_names = []
        for id, compound in dictitems(self.input.compound):
            compounds.append(compound)
            self.test_names.append(id)
        (self.compounds, self.feature_names) = _get_feature_names(compounds)
        self.load_testing_set()
        #self.impute_values()

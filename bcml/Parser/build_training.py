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
import warnings
import sys
from sklearn.preprocessing import Imputer
from Boruta import boruta_py
from sklearn.ensemble import RandomForestClassifier
from KNNImpute.knnimpute import (
    knn_impute_optimistic,
)
from scipy import stats
from math import sqrt
try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict


def dictitems(dict):
    if sys.version_info[0] >= 3:
        return dict.items()
    else:
        return dict.iteritems()


def verbose_print(verbose, line):
    if verbose:
        print(line)

_possible_features = ('experimentalhash', 'binhash', 'padelhash', 'userhash')


def _convert_predictor(predictors, split_value):
    """"This function discretizes the predictors. This is currently
    set to handle a binary class, future versions will handle multi-class
    predictors"""
    preds = np.zeros(len(predictors), dtype=int)
    for i, value in np.ndenumerate(predictors):
        if value >= split_value:
            preds[i] = 1
    return preds


def _get_feature_names(compounds):
    """This function handles collecting the feature names"""
    feature_names = {}
    experimental_features = []
    for compound in compounds:
        for feature in sorted(_possible_features):
            if feature == 'experimentalhash':
                keys = sorted(compound[feature].keys())
                for feat in keys:
                    experimental_features.append(feat)
            if feature in compound.keys():
                keys = sorted(compound[feature].keys())
                for feat in keys:
                    feature_names[feat] = 1
                    compound[feat] = compound[feature][feat]
    return (compounds, sorted(feature_names.keys()), experimental_features)


class Process(object):
    """This file reads a training file"""
    def _load_training_set(self):
        """This function takes the features and
        compounds and loads them into a numpy array
        """
        for index, value in np.ndenumerate(self.train):
            compound = self.compounds[index[0]]
            feature = list(self.feature_names)[index[1]]
            if (feature in compound.keys()) and (compound[feature] is not "")\
               and (compound[feature] != "NULL")\
               and (compound[feature] != "False"):
                self.train[index] = float(compound[feature])
            else:
                self.train[index] = np.nan

    def check_errors(self, distance, verbose):
        """This function identifies values that may be in error in
        the training set"""
        feature_mask = []
        distance[distance == np.inf] = 0
        for i, feature in enumerate(self.feature_names):
            if feature in self.experimental_features:
                feature_mask.append(i)
        mask = np.ones_like(self.train)
        mask[:, feature_mask] = 0
        masked_training = np.ma.masked_array(self.train, mask)
        column_means = np.ma.mean(masked_training, axis=0)
        column_std = np.ma.std(masked_training, axis=0)
        column_sqrt = sqrt(masked_training.shape[1])
        critical = stats.t.isf([0.0001], masked_training.shape[1])[0]
        n = masked_training.shape[1]
        tau = (critical * (n - 1)) / (sqrt(n)*sqrt(n - 2. + (critical**2)))
        print(tau)
        weighted_sum = np.sum(distance, axis=0)
        for i in feature_mask:
            x_vector = masked_training[:,i]
            weighted_vector = distance.dot(x_vector)
            weighted_average = weighted_vector / weighted_sum
            gamma = abs(x_vector - weighted_average)
            weighted_sd = (distance[:,i] * ((x_vector - weighted_average)**2)) / (weighted_sum - ((distance[:,i]**2)/weighted_sum)) 
            d = tau * column_std[i]
            error_values = np.where(gamma > (tau * column_std[i]))
            print(self.feature_names[i], error_values)


    def feature_selection(self, verbose, seed=False):
        """This function runs Boruta feature selection to remove
        unimportant features from the data"""
        if verbose:
            verbose = 2
        '''Boruta cannot handle missing values. Either run impute_values
        before feature_selection, or the following function runs mean
        imputation prior to running Boruta'''
        if np.any(np.isnan(self.train)):
            warnings.warn('RandomForestClassifier requires no missing data,\
                           features being imputed by mean')
            X = self.train
            imp = Imputer(missing_values='NaN', strategy='mean', axis=0)
            imp.fit(X)
            self.train = imp.transform(X)
        rf = RandomForestClassifier(n_jobs=-1, oob_score=True, max_depth=5)
        feat_selector = boruta_py.BorutaPy(rf, n_estimators='auto',
                                           verbose=verbose, seed=seed)
        feat_selector.fit(self.train, self.predictors)
        self.feature_support = feat_selector.support_
        filtered_names = [i for indx, i in enumerate(self.feature_names) if self.feature_support[indx]]
        self.feature_names = filtered_names
        self.train = feat_selector.transform(self.train)

    def impute_values(self, distance, verbose, k=5):
        """This function handles the missing values from
        the training set and estimates their value, based on
        the mean and reloads them into the training set"""
        verbose_print(verbose, 'Imputing using KNN strategy')
        X = self.train
        missing_mask = np.isnan(X)
        '''First impute uing knn optimistic'''
        impute = knn_impute_optimistic(X, missing_mask=missing_mask,
                                       distance=distance, k=k)
        X = impute.astype(float)
        '''For features with a small number of features, use mean
        imputation to remove NaN values'''
        imp = Imputer(missing_values='NaN', strategy='mean', axis=0)
        imp.fit(X)
        self.train = imp.transform(X)

    def __init__(self, model_input, split_value=False, verbose=False):
        """This initialization function handles the heavy
        work of loading the features and processing the
        compounds"""
        self.input = model_input
        compounds = []
        predictors = []
        weights = []
        self.input.compound = OrderedDict(sorted(self.input.compound.items(), key=lambda t: t[0]))
        for id, compound in dictitems(self.input.compound):
            compounds.append(self.input.compound[id])
            predictors.append(self.input.compound[id]['predictor'])
            weights.append(self.input.compound[id]['weight'])
        predictor_values = np.array(predictors, '|S4').astype(np.float)
        weight_values = np.array(weights, '|S4').astype(np.float)
        self.weights = weight_values
        if split_value:
            self.predictors = _convert_predictor(predictor_values, split_value)
        else:
            print_line = "Splitting at " + str(np.median(predictor_values))
            verbose_print(verbose, print_line)
            self.predictors = _convert_predictor(predictor_values,
                                                 np.median(predictor_values))
        self.rows = len(self.predictors)
        self.compounds, self.feature_names, self.experimental_features = _get_feature_names(compounds)
        self.columns = len(self.feature_names)
        '''Initialize the training array'''
        self.train = np.zeros((self.rows, self.columns,), dtype=np.float64)
        '''Load the training array'''
        self._load_training_set()

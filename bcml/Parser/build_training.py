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


def dictitems(dict):
    if sys.version_info[0]>=3:
        return dict.items()
    else:
        return dictitems(dict)


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
    def _load_training_set(self):
        """This function takes the features and
        compounds and loads them into a numpy array
        """
        for index, value in np.ndenumerate(self.train):
            compound = self.compounds[index[0]]
            feature = self.feature_names[index[1]]
            if (feature in compound.keys()) and (compound[feature] is not "")\
               and (compound[feature] != "NULL") and (compound[feature] != "False"):
                self.train[index] = float(compound[feature])
            else:
                self.train[index] = np.nan

    def feature_selection(self, verbose, seed=False):
        """This function runs Boruta feature selection to remove
        unimportant features from the data"""
        if verbose:
            verbose = 2
        '''Boruta cannot handle missing values. Either run impute_values
        before feature_selection, or the following function runs mean
        imputation prior to running Boruta'''
        if np.any(np.isnan(self.train)):
            warnings.warn('RandomForestClassifier requires no missing data, features being imputed by mean')
            X = self.train
            imp = Imputer(missing_values='NaN', strategy='mean', axis=0)
            imp.fit(X)
            self.train = imp.transform(X)
        rf = RandomForestClassifier(n_jobs=-1, max_depth=5)
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


    def __init__(self, model_input, split_value=False):
        """This initialization function handles the heavy
        work of loading the features and processing the
        compounds"""
        self.input = model_input
        compounds = []
        predictors = []
        
        for id, compound in dictitems(self.input.compound):
            compounds.append(self.input.compound[id])
            predictors.append(self.input.compound[id]['predictor'])
        predictor_values = np.array(map(float, predictors))
        if split_value:
            self.predictors = _convert_predictor(predictor_values, split_value)
        else:
            self.predictors = _convert_predictor(predictor_values, np.median(predictor_values))
        self.rows = len(self.predictors)
        self.compounds, self.feature_names = _get_feature_names(compounds)
        self.columns = len(self.feature_names)
        '''Initialize the training array'''
        self.train = np.zeros((self.rows, self.columns,), dtype=np.float64)
        '''Load the training array'''
        self._load_training_set()

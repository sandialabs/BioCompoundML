from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import OneHotEncoder, Imputer
from sklearn.cross_validation import train_test_split
import numpy as np
import warnings


class Train(object):
    """docstring for TrainModel"""
    def preprocess_model(self):
        '''This allows preprocessing using logistic regression'''
        X_train, X_train_lr, y_train, y_train_lr = train_test_split(self.train,
                                                                    self.predictors,
                                                                    test_size=0.5)
        encode = OneHotEncoder()
        logistic = LogisticRegression()
        self.clf = RandomForestClassifier(n_estimators=512,
                                          oob_score=True, n_jobs=-1)
        self.clf.fit(X_train, y_train)
        encode.fit(self.clf.apply(X_train))
        self.predmodel = logistic.fit(encode.transform(self.clf.apply(X_train_lr)), y_train_lr)

    def train_model(self):
        '''This is standard model training'''
        '''For RandomForestClassifier to work their must be no nan values, one way of
        handling this is to use the --impute option. This uses mean imputation, which
        is the least information imputer, imputation is done by feature
        '''
        if np.any(np.isnan(self.train)):
            warnings.warn('RandomForestClassifier requires no missing data, features being imputed by mean')
            X = self.train
            imp = Imputer(missing_values='NaN', strategy='mean', axis=0)
            imp.fit(X)
            self.train = imp.transform(X)
        self.clf = RandomForestClassifier(n_estimators=512,
                                          oob_score=True, n_jobs=-1)
        self.predmodel = self.clf.fit(self.train, self.predictors)

    def __init__(self, train):
        self.train = train.train
        self.predictors = train.predictors
        self.features = train.feature_names

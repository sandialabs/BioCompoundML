"""
This class handles distance calculation

Attributes:
    input: list of features
    output: transformed numpy array
"""
from sklearn.feature_extraction import DictVectorizer
from sklearn.metrics.pairwise import pairwise_distances


class Distance(object):
    def build_matrix(self, feature_dict):
        self.distance = pairwise_distances(X=feature_dict, metric=self.metric,
                                           n_jobs=self.n_jobs)

    def __init__(self, input_data, key_data, method='jaccard', n_jobs=-1):
        self.input = input_data
        vector = DictVectorizer(sparse=False)
        vectorized = vector.fit_transform(input_data)
        self.keys = key_data
        self.metric = method
        self.n_jobs = n_jobs
        self.build_matrix(vectorized)

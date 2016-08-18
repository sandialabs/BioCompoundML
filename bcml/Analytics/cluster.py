'''
This module handles the clustering of compounds in the training set
and also the processing of test molecules.
'''
from __future__ import print_function
from scipy import stats
import numpy as np
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.metrics.pairwise import pairwise_distances
from sklearn import cluster
from sklearn.ensemble import RandomTreesEmbedding
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.covariance import EllipticEnvelope
from sklearn.svm import OneClassSVM


class Clustering():
    def __init__(self, compounds, output=False, seed=False):
        self.compounds = compounds
        self.count = 0
        self.count_1 = 0
        self.output = output
        self.tools = clustertools()
        if self.output is not False:
            self.figures = clusterfigures(self.compounds)
        self.testcompound = []
        self.seed = seed

    def cluster_training(self, train, distance=False):
        '''
        This is the basic clustering function
        '''
        self.train_matrix = train.train
        '''
        Step one is to make sure that their is a distance matrix in place.
        It is best to feed an existing distance matrix if one is available.
        '''
        if distance is False:
            self.p_feat_matrix = self.tools.pairwise_distance_matrix(train.train, 'jaccard')
        else:
            self.p_feat_matrix = distance
        '''
        Step two is to cluster your data using a random trees embedding. This a
        random ensemble of trees. This is a transformation on the data, into a
        high dimensional, sparse space
        '''
        self.clf = RandomTreesEmbedding(n_estimators=512, random_state=self.seed, max_depth=5)
        #self.clf.fit(self.train_matrix)
        X_transformed = self.clf.fit_transform(self.train_matrix)
        '''
        Step three performs truncated SVD (similar to PCA). It operates on the sample
        vectors directly, rather than the covariance matrix. It takes the first two
        components. Essentially this reduces the sparse embedding to a low dimensional
        representation.
        '''
        self.svd = TruncatedSVD(n_components=2)
        self.svd.clf = self.svd.fit(X_transformed)
        self.model = self.svd.clf.transform(X_transformed)
        '''
        The next step is to take the transformed model and the original dataset and
        determine the max silhouette_score of clusters
        '''
        (self.cluster_assignment,
         self.cluster_num,
         self.cluster_score) = self.tools.identify_accurate_number_of_clusters(self.model, self.compounds)
        self.individualclusters = []
        '''
        The individual datapoints are assessed with regard to the best clustering scheme
        '''
        for i in range(self.cluster_num):
            self.individualclusters.append([])
            for j in range(len(self.cluster_assignment)):
                if self.cluster_assignment[j] == i:
                    self.individualclusters[i].append(self.model[j, :])
            self.individualclusters[i] = np.array(self.individualclusters[i])
        '''
        Finally, this clustering scheme is used to generate a one class Support
        Vector Machine decision boundary.
        '''
        (self.clf_OCSVM,
         self.OCSVM_model) = self.tools.determine_test_similarity(self.individualclusters)

    def cluster_testing(self, testing):
        '''Create RandomTreesEmbedding of data'''
        clf = RandomTreesEmbedding(n_estimators=512, random_state=self.seed, max_depth=5)
        '''Fit testing data to training model'''
        clf.fit = self.clf.fit(testing)
        X_transformed = self.clf.fit_transform(testing)
        n_components = 2
        '''SVD transform data'''
        svd = TruncatedSVD(n_components=n_components)
        svd.clf = svd.fit(X_transformed)
        svd.model = svd.clf.transform(X_transformed)
        '''Train transformed data using original model'''
        train_transformed = clf.fit.transform(self.train_matrix)
        train_model = svd.clf.transform(train_transformed)
        '''Generate One Class SVM rejection criteria'''
        (clf_OCSVM_t, OCSVMmodel_t) = self.tools.determine_testing_data_similarity(train_model)
        predicted = []
        '''Remove testing compounds outside rejection margin'''
        for i in range(len(svd.model)):
            p = OCSVMmodel_t.predict(svd.model[i, :].reshape(1, -1))
            pred = OCSVMmodel_t.decision_function(svd.model[i, :].reshape(1, -1)).ravel()
            if (p == 1):
                predicted.append(i)
        return predicted


class clusterfigures():
    def __init__(self, compounds):
        pass


class clustertools():
    def __init__(self):
        pass

    def pairwise_distance_matrix(self, feature_matrix, dist_matrix_type):
        pair_distance_matrix = pairwise_distances(feature_matrix,
                                                  metric=dist_matrix_type)
        return pair_distance_matrix

    def PCA(self, feature_matrix):
        model = PCA(n_components=2).fit(feature_matrix)
        transformed = model.transform(feature_matrix)
        return model, transformed

    def identify_accurate_number_of_clusters(self, model, compounds, max_range=3):
        silhouette_avg_scores = []
        for n_cluster in range(2, max_range):
            assigned_cluster = cluster.KMeans(n_clusters=n_cluster,
                                              n_init=20).fit_predict(model)
            silhouette_avg = silhouette_score(model, assigned_cluster)
            silhouette_avg_scores.append(silhouette_avg)
        max_silhouette_score = max(silhouette_avg_scores)
        index_max_score = silhouette_avg_scores.index(max_silhouette_score)
        final_cluster_num = range(2, max_range)[index_max_score]
        final_assigned_cluster = cluster.KMeans(n_init=20,
                                                n_clusters=final_cluster_num).fit_predict(model)
        final_sample_sil_vals = silhouette_samples(model, final_assigned_cluster)
        return final_assigned_cluster, final_cluster_num, final_sample_sil_vals

    def model_2_determine_test_data_similarity(self,model):
        clf_EE={}
        model_EE={}
        for i in range(len(model)):
            clf=EllipticEnvelope(contamination=0.01,support_fraction=1)
            clf_EE[i]=clf
            EEmodel=clf.fit(model[i])
            model_EE[i]=EEmodel
        return clf_EE,model_EE

    def determine_test_similarity(self, model):
        clf_OCSVM = {}
        model_OCSVM = {}
        for i in range(len(model)):
            clf = OneClassSVM(kernel='rbf', nu=0.1, gamma=.023)
            clf_OCSVM[i] = clf
            OCSVMmodel = clf.fit(model[i])
            model_OCSVM[i] = OCSVMmodel
        return clf_OCSVM, model_OCSVM

    def determine_testing_data_similarity(self, model):
        clf = OneClassSVM(kernel='rbf', nu=0.03, gamma=.05)
        clf_OCSVM = clf
        OCSVMmodel = clf.fit(model)
        model_OCSVM = OCSVMmodel
        return clf_OCSVM, model_OCSVM

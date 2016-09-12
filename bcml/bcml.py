'''Main class'''
from __future__ import print_function
import argparse
from Parser import read_training as rt
from Parser import read_testing as rtest
from Chemoinformatics import fingerprints as fp
from Chemoinformatics import chemofeatures as cf
from Chemoinformatics import experimental as exp
from Analytics import cluster as cl
from Chemoinformatics import user as usr
from Parser import build_training as bt
from Parser import build_testing as btest
from Train import train_model as tm
from Analytics import cross_validate as cv
import pickle
import numpy as np
from Distance import distance as ds
try:
    from urllib.parse import urlparse
except ImportError:
    from urlparse import urlparse
from PubChemUtils import pubchempy_utils as pcp
from copy import deepcopy
import sys
import warnings

_chunks = 10
_id = 'PubChem'
_test_results = "tested_metacyc.out"
_try_count = 3
_plot = False


class Training(object):
    pass


class ArgumentExceptions(Exception):
    def __init__(self, arg):
        self.msg = arg

    def __str__(self):
        return repr(self.msg)


def dictitems(dict):
    if sys.version_info[0] >= 3:
        return dict.items()
    else:
        return dict.iteritems()


def verbose_print(verbose, line):
    if verbose:
        print(line)


def parse_arguments():
    '''
    The point of this program is to apply machine learning to a set of
    chemicals with know properties Default API (training)
    Usage:
    python train_model.py --input training.txt --basic --model model.pkl
    python train_model.py --input properties.txt --user --model model.pkl
    python train_model.py --input props.txt --basic --user --model model.pkl
    :params i/input The required input file
    :params train If True then run in training mode
    :params test If True then run in testing model_file
    :params m/model The output name of the model
    :params d/datain The name of the input datafile: either .model or .features
    from prior runs
    :params o/dataout The name of the saved pkl files. These will be appended
    with .model and .features
    :params p/pred The name of the parameter in the input file that is trained,
    which allows other features to be provided in this file
    :params x/proxy The http for the proxy, if one is used
    :params cluster Cluster the training data
    :params cluster_testing Filter test compounds based on training clusters
    :params s/split_value Split the data in classes based on the provided value
    :params r/random User provided random seed for numpy randomization
    :params v/verbose Specify the level of output
    :params e/experimental If True then extract experimental PubChem features
    :params f/fingerprint If True then extract fingerprint PubChem features
    :params c/chemofeatures If True then extract sdfs from PubChem and
    descriptors from Padel-Descriptor
    :params u/user If True then extract additional features from input file
    :returns The parsed arguments
    It is important to know that there is a hierarcy of processing
    (e.g., training):
                                  |
                              input_file
                                  |
                                random
                                  |
                                proxy
                                /
                            train
                              |
                             pred
                            /    \
                         datain   no-datain
                         /    \         \
                       model  features   fingerprint (Y/N)
                       /        \         \
            cross-validate  build_training user (Y/N)
                                  \         \
                                 train       experimental (Y/N)
                                    \         \
                             cross-validate    chemofeatures (Y/N)
                                                 /       \
                                            median   split_value
                                                \        /
                                                 dataout
                                                    |
                                                  train
                                                    |
                                                 dataout
                                                    |
                                                 cluster
    '''
    parser = argparse.ArgumentParser(description="This program trains,\
                                     evaluates and tests a Random \
                                     Forest Classifer using list of chemicals\
                                     and a measured property.")
    grp_input = parser.add_mutually_exclusive_group(required=True)
    grp_input.add_argument('-i', '--input', help='Input file', type=str)
    grp_input.add_argument('-d', '--datain', help="Data input file", type=str)
    parser.add_argument('-t', '--test_input', help="Data input for test",
                        type=str)
    parser.add_argument('--train', help="Train the model", action="store_true")
    parser.add_argument('--test', help="Test the model", action="store_true")
    parser.add_argument('-m', '--model', help="Model output file", type=str)
    parser.add_argument('-o', '--dataout', help="Data output file", type=str)
    parser.add_argument('--pred', help="Prediction variable", type=str)
    parser.add_argument('-x', '--proxy', help="Proxy", type=str)
    parser.add_argument('--cluster', help="Run clustering",
                        action="store_true")
    parser.add_argument('-s', '--split_value', help="Value to split training \
                        data, if not specified median value of data is used",
                        type=float)
    parser.add_argument('--random', help="Seed for randomization", type=int)
    parser.add_argument('-v', '--verbose', help="Level of output",
                        action="store_true")
    parser.add_argument('-e', '--experimental',
                        help="Extract experimental features",
                        action="store_true")
    parser.add_argument('-f', '--fingerprint',
                        help="Extract fingerprint features",
                        action="store_true")
    parser.add_argument('-c', '--chemofeatures',
                        help="Extract PaDEL-Descriptors",
                        action="store_true")
    parser.add_argument('-u', '--user',
                        help="Extract user-provided features from input file",
                        action="store_true")
    parser.add_argument('--distance', help="Create a distance matrix of data",
                        action="store_true")
    parser.add_argument('--impute',
                        help="Impute missing values using K-NN imputation",
                        action="store_true")
    parser.add_argument('--selection',
                        help="Run Boruta Feature Selection to reduce\
                        uncharacterizing features", action="store_true")
    parser.add_argument('--cv',
                        help="Cross-Validate model using 50\% hold-out data",
                        action="store_true")
    parser.add_argument('--weight',
                        help="Insert sample weights",
                        action="store_true")
    return parser.parse_args()


def define_random_seed(args):
    '''Define random seed'''
    if args.random:
        '''
        Users may add a random seed. This helps reproducibility in parts of
        the machine learning that require reandomization
        :param seed: int provided in the program arguments
        '''
        assert args.random > 0, "Random seed must be a positive integer"
        assert isinstance(args.random, int), "Random seed must be a positive integer"
        seed = args.random
        np.random.seed(seed)
        return seed
    else:
        return False


def define_proxy(args):
    '''Define the proxy URL'''
    if args.proxy:
        '''
        Users may add a proxy for interacting with the web.
        :param proxy: web address
        '''
        result = urlparse(args.proxy)
        assert result.scheme, "Proxy must be a web address"
        return args.proxy
    else:
        return False


def existing_training_model(args, seed):
    '''
    Input data are used to prepopulate the model, based on prior runs
    :param datain Input data files are stored as models (if .model) or
    features (if .features) they are pickeled model and feature files
    '''
    verbose_print(args.verbose, "Using existing model")
    trained_model = Training()
    if args.datain.endswith('.cluster'):
        with open(args.datain, 'rb') as fid:
            cluster = pickle.load(fid)
            trained_model.cluster = cluster
            trained_model.model = cluster.model
    elif args.datain.endswith('.model'):
        with open(args.datain, 'rb') as fid:
            model = pickle.load(fid)
            trained_model.model = model
        if args.cluster:
            cluster = cl.Clustering(model.training_data.compound,
                                    seed=args.random)
            cluster.cluster_training(model)
            trained_model.cluster = cluster
    elif args.datain.endswith('.features'):
        with open(args.datain, 'rb') as fid:
            train = pickle.load(fid)
            train = bt.Process(train, seed=seed)
            model = tm.Train(train)
            model.train_model()
    if args.cv:
        verbose_print(args.verbose, "Running cross-validational analysis")
        cv.Analysis(model, seed)
    return trained_model


def check_features(args):
    '''
    This function checks what kinds of features have been defined for this
    run
    '''
    [user, experimental, chemofeatures, fingerprint] = [False] * 4
    if args.user:
        user = True
    if args.experimental:
        experimental = True
    if args.chemofeatures:
        chemofeatures = True
    if args.fingerprint:
        fingerprint = True
    return [user, experimental, chemofeatures, fingerprint]


def add_pubchem_features(compounds, args, user=False, proxy=False,
                         fingerprint=False, experimental=False,
                         chemofeatures=False, id_name='PubChem', chunks=False):
    '''This function loads the pubchem features and downloads the data'''
    predictors = False
    if compounds.predictors:
        predictors = compounds.predictors
    if compounds.weights:
        weights = compounds.weights
    collected_data = pcp.Collect(compounds.compounds, fingerprint=fingerprint,
                                 xml=experimental, sdf=chemofeatures,
                                 proxy=args.proxy, user=user, id_name=id_name,
                                 chunks=chunks, try_count=_try_count,
                                 verbose=args.verbose, predictors=predictors,
                                 weights=weights)
    return collected_data


def collect_distance_matrix(collected_data):
    '''If clustering the data, retain the original training_data, without
    removing the redundant features'''
    original_data = deepcopy(collected_data)
    original_data = fp.Update(original_data, remove_static=False)
    original_data.update()
    fingerprint_vector = list()
    key_list = list()
    for key, value in dictitems(original_data.compound):
        fingerprint_vector.append(value['binhash'])
        key_list.append(key)
    distance = ds.Distance(fingerprint_vector, key_list).distance
    return distance


def extract_features(collected_data, args, user=False, fingerprint=False,
                     experimental=False, chemofeatures=False,
                     remove_static=True):
    '''This function collects and extracts features from the raw user and
    PubChem data'''
    if fingerprint is True:
        '''This part of the code converts the CACTVS fingerprints into features
        it also removes the fully redundant features, which for many chemical
        properties may be a large portion'''
        collected_data = fp.Update(collected_data, remove_static=remove_static,
                                   verbose=args.verbose)
        collected_data.update()
    if user is True:
        '''This part of the code extracts the user defined features'''
        collected_data = usr.Update(collected_data,
                                    remove_static=remove_static,
                                    verbose=args.verbose)
        collected_data.update()
    if experimental is True:
        '''This part of the code extracts PubChem experimental and computed
        features from the PubChem xml files'''
        collected_data = exp.Update(collected_data,
                                    remove_static=remove_static,
                                    verbose=args.verbose)
        collected_data.update()
    if chemofeatures is True:
        '''This code runs PaDEL-Descriptor and extracts relevant 2-D
        and 3-D chemical descriptors'''
        collected_data = cf.Update(collected_data,
                                   remove_static=remove_static,
                                   verbose=args.verbose)
        collected_data.update(padel=True)
    return collected_data


def report_model_validation(model, args):
    verbose_print(args.verbose, "Running cross-validational analysis")
    analysis = cv.Analysis(model, args.random)
    analysis.cross_validate()
    if _plot:
        analysis.plot_random_ROC()
    analysis.feature_importances()


def train_model(args, seed, proxy, pred):
    trained_model = Training()
    if args.train:
        '''
        The program is essentially run in one of two mutually exclusive modes
        (training or test)
        :param train if True, being parsing and training model file
        '''
        verbose_print(args.verbose, "Training model")
        if args.datain:
            warnings.warning("WARNING: The pickle datatype is inherently\
                             insecure. A quick question: do you trust the\
                             source of your model? Pickle files can contain\
                             corrupt code and executable commands.\
                             They can take over your computer and install\
                             malicious code on your computer or server. Use\
                             caution! Your best bet is to train your own\
                             models and run those! Use --datain at your own\
                             risk")
            continue_program = raw_input("Press [Y/y] if you want to continue")
            if continue_program in ['Y', 'y']:
                trained_model = existing_training_model(args, seed)
            else:
                exit()
        else:
            distance = False
            training_data = False
            verbose_print(args.verbose, "Reading training set")
            (user, experimental, chemofeatures, fingerprint) = check_features(args)
            if (args.distance is True) or (args.cluster is True) or (args.impute is True):
                '''These functions all require a distance matrix, which is best
                collected using the fingerprint data'''
                fingerprint = True
            training = rt.Read(args.input, pred, user=user, id_name=_id,
                               weights=args.weight)
            '''This block of code generally works on feature collection and
            parsing, including the removal of fully redundant features. The
            difference between remove_static=True and False is whether or not
            to get rid of fully redundant features. Since the distance matrix
            is the same, regardless, it is run using original data'''
            training_data = add_pubchem_features(training, args, user=user,
                                                 proxy=proxy,
                                                 fingerprint=fingerprint,
                                                 experimental=experimental,
                                                 chemofeatures=chemofeatures,
                                                 id_name=_id, chunks=_chunks)
            if (args.cluster is True) or (args.distance is True) or (args.impute is True):
                verbose_print(args.verbose, "Creating distance matrix")
                '''Collect distance matrix using the original dataset'''
                distance = collect_distance_matrix(training_data)
            '''Extract features from the user and PubChem data'''
            verbose_print(args.verbose, "Extracting features")
            training_data = extract_features(training_data, args, user=user,
                                             fingerprint=fingerprint,
                                             experimental=experimental,
                                             chemofeatures=chemofeatures,
                                             remove_static=True)
            '''Discretize the y-values for the the classification process.
            If no split value is provided then the default for the program
            is to break the value at the median
            '''
            if training_data.compound:
                train = bt.Process(training_data, split_value=args.split_value,
                                   verbose=args.verbose)
                if args.impute is True:
                    train.impute_values(distance=distance,
                                        verbose=args.verbose)
                if args.selection is True:
                    train.feature_selection(verbose=args.verbose,
                                            seed=args.random)
                '''If dataout parameter is set, it prints to pickle a file
                containing the features that were extracted. In later runs
                this can be specified as the data input using the datain
                parameter
                '''
                if args.dataout:
                    features_file = args.dataout + ".features"
                    with open(features_file, 'wb') as fid:
                        pickle.dump(train, fid)
                '''This is where the model is actually trained in the tm module'''
                model = tm.Train(train)
                model.train_model()
                trained_model.model = model
                '''If dataout parameter is set, it prints to pickle a file
                containing the RF model. In later runs this can be specified
                as the data input using the datain parameter
                '''
                if args.dataout:
                    model_file = args.dataout + ".model"
                    with open(model_file, 'wb') as fid:
                        pickle.dump(model, fid)
                if args.cv:
                    report_model_validation(model, args)
                if args.cluster:
                    cluster = cl.Clustering(training_data.compound, seed=args.random)
                    cluster.cluster_training(model)
                    trained_model.cluster = cluster
                    if args.dataout:
                        cluster_file = args.dataout + ".cluster"
                        with open(cluster_file, 'wb') as fid:
                            pickle.dump(cluster, fid)
    else:
        trained_model = False
    return trained_model


def test_model(trained_model, args, seed=False, proxy=False, pred=False):
    if trained_model is False:
        assert args.datain, "If the model hasn't been trained, input data\
                             must be specified"
        model_file = args.datain + ".model"
        trained_model = False
        cluster = False
        with open(model_file, 'rb') as fid:
            trained_model = pickle.load(fid)
        if args.cluster:
            cluster_file = args.datain + ".cluster"
            with open(cluster_file, 'rb') as fid:
                cluster = pickle.load(fid)
    else:
        if args.test_input:
            (user, experimental, chemofeatures, fingerprint) = check_features(args)
            testing = rtest.Read(args.test_input, user=user, id_name=_id)
            testing.predictors = False
            if (args.distance is True) or (args.cluster is True) or (args.impute is True):
                '''These functions all require a distance matrix, which is best collected
                using the fingerprint data'''
                fingerprint = True
            verbose_print(args.verbose, "Adding PubChem Features")
            testing_data = add_pubchem_features(testing, args, user=user,
                                                proxy=proxy,
                                                fingerprint=fingerprint,
                                                experimental=experimental,
                                                chemofeatures=chemofeatures,
                                                id_name=_id,
                                                chunks=_chunks)
            distance = collect_distance_matrix(testing_data)
            verbose_print(args.verbose, "Extracting features")
            testing_data = extract_features(testing_data, args, user=user,
                                            fingerprint=fingerprint,
                                            experimental=experimental,
                                            chemofeatures=chemofeatures,
                                            remove_static=False)
            test = btest.Process(testing_data, trained_model.model.features,
                                 trained_model.model.train)
            test.impute_values(distance)
    if trained_model.cluster or args.cluster:
        verbose_print(args.verbose,
                      "Rejecting based on clustering of training")
        cluster = trained_model.cluster
        testmatrix = []
        testcompounds = []
        for i in range(len(test.test[:, 0])):
            testmatrix.append(test.test[i, :])
            testcompounds.append(test.test_names[i])
        testmatrix = np.array(testmatrix)
        predicted = cluster.cluster_testing(testmatrix)
        for pred in predicted:
            print(test.test_names[pred])
    prediction = trained_model.model.clf.predict_proba(test.test)
    for i, name in enumerate(test.test_names):
        print(name, prediction[i, ])


def define_pred(args):
    '''
    Users specify the prediction value in the input file
    :param pred: prediction value
    '''
    assert args.pred, "Prediction value must be defined"
    return args.pred


def check_arguments(args):
    '''
    Checks to make sure train or test or both
    modes
    :param args.train
    :param args.test
    '''
    parser = argparse.ArgumentParser()
    if args.cluster:
        warnings.warning('Cluster module is still not fully functional')
    if not (args.train or args.test):
        parser.error('No action requested, add --train or --test')
    if (args.test) and not (args.test_input):
        parser.error("If testing, must specify test data, use -t/--test_input\
                     <<DATAFILE>>")


def main():
    args = parse_arguments()
    check_arguments(args)
    seed = define_random_seed(args)
    proxy = define_proxy(args)
    pred = define_pred(args)
    if args.train:
        trained_model = train_model(args, seed, proxy, pred)
    if args.test:
        test_model(trained_model, args, seed, proxy, pred)

if __name__ == '__main__':
    main()
    exit()

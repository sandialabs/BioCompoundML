.. BioCompoundML documentation master file, created by
   sphinx-quickstart on Wed Aug  3 21:40:38 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to BioCompoundML's documentation!
=========================================

BioCompoundML provides a chemoinformatic tool for predicting chemical properties using machine learning.

Purpose
=======

Quantitative Structure Property/Activity Relationships (QSPRs and QSARs) often attempt to determine as exactly as possible the exact value of a given chemical property. Tools for predicting these properties are incredibly useful, but are often limited in one of two regards:
	1) They lack generality - they cannot be readily rebuilt for an arbitrary set of chemical properties.
	2) They can require highly sophisticated and expensive computation.

If however, our question is discrete, for instance, binary classification above or below a certain threshold (e.g., melting point at room temperature), then the problem can be reframed as a machine learning classification problem and solved rather quickly.

We came to this problem with a particular interest in mind, one we feel is common in cheminformatics. The rapid screening of a large number of compounds for multiple chemical properties can oftentimes be handled efficiently and effectively using a classification paradigm.

Common BioCompoundML Workflow
================================

Input file
----------
BioCompoundML starts with the training of a model using random forest classification. As such it requires an initial training file. This file provides a list of compounds and with a measured value. This can take a variety of forms, however the easiest one uses a tab-delimited file with the name of the compound, a PubChem identifier and a measured value. ::

	#Name	RON	PubChem
	1-Butene	98.8	7844
	1-Ethyl-3-Methylcyclopentane	57.6	19502
	1-Heptene	54.5	11610
	1-Hexene	76.4	11597
	1-Isopropyl-4-methylcyclohexane	67.3	7459
	1-Methyl-1-ethylcyclohexane	68.7	35411

There are a few important things to recognize. 1) The header line is essential and must start with #. 2) Name and PubChem are important and must be capitalized as they are here. 3) The format of the file is tab-delimited text. 4) You can output tab-delimited text from Excel, but do it with caution, export from MSOffice products can have unexpected effects.

BioCompoundML uses the NCBI PubChem API heavily. There are ways of handling CAS numbers, but PubChem ID (CID) is the easiest and most direct. Providing CAS requires a separate call to NCBI to retrieve CIDs.

The user must specify the name of the feature being trained, in the above case ``RON``. If a split-value isn't provided, then BioCompoundML splits on the median.

Feature Collection
------------------
User-provided features
^^^^^^^^^^^^^^^^^^^^^^
The next step in the workflow is to collect Cheminformatic features. There are a variety of these. One is to simply use 'user' provided features. Below is an example of the original training file including OH Rate Constant. ::

	#Name	RON	PubChem	OH_Rate_Constant
	Methyl acetate	120	6584	0.2598
	O-Xylene	120	7237	6.5119
	Ethyl acetate	118	8857	1.7038
	Ethyl buanoate	115.4	7762	3.3339
	Propylbenzene	111	7668	7.31

When additional columns are included in the training file and 'user' is selected as a parameter (see :doc:`script` for examples of the parameters available), this feature is added to the model. This is particularly useful when you have private or licensed values. It is important to remember that if you wish to use these features for prediction, you will need to provide them in both the training and testing datasets.

PubChem features
^^^^^^^^^^^^^^^^
In addition to user-provided features, BioCompoundML also collects features directly from PubChem. These include computationally predicted/experimentally measured features, 881 SMILES fingerprints and Structural Data Files (SDFs). 

The choice of which of these features to collect are specified in the ``--fingerprints``, ``--experimental`` and ``--chemoinformatics`` parameters for ``bcml.py``

Fingerprints
""""""""""""
Fingerprints directly collects a CACTVS string and converts this to 881 binary SMILES features (e.g., C > 4 or C(-C)(-C)(=C)), a full list can be found at ftp://ftp.ncbi.nlm.nih.gov/pubchem/specifications/pubchem_fingerprints.txt

Experimental and Computationally-Predicted features
"""""""""""""""""""""""""""""""""""""""""""""""""""
The experimental/estimated features that are used by the feature extraction package include experimentally measured properties (e.g., melting point, boiling point, vapor pressure); inferred structural features (e.g., rotatable bond count, heavy atom count); chemical properties (e.g., molecular weight, formal charge) and inferred chemical properties (e.g., XLogP3 – which estimates the Octanol-Water partitioning – a property directly related to hydrophobicity). The file ``Chemoinformatics/feature_list.txt`` contains the features that are collected from PubChem. Additional features can be added or removed from this file, depending on the importance of the features. The standard file ships with the following features selected ::

	Density
	Vapor Density
	Boiling Point
	Hydrogen Bond Donor Count
	Rotatable Bond Count
	XLogP3
	Flash Point
	Formal Charge
	Undefined Atom Stereocenter Count
	Auto-Ignition
	Molecular Weight
	Hydrogen Bond Acceptor Count
	XLogP3-AA
	LogP
	Defined Atom Stereocenter Count
	Complexity
	Vapor Pressure
	Covalently-Bonded Unit Count
	Isotope Atom Count
	Undefined Bond Stereocenter Count
	Heavy Atom Count
	Exact Mass
	Monoisotopic Mass
	Topological Polar Surface Area
	Melting point
	Defined Bond Stereocenter Count

Additional chemoinformatic features
"""""""""""""""""""""""""""""""""""
Additionally, this package also has the capacity to retrieve substance data files (SDFs) from NCBI. These may be useful in downstream QSPR/QSAR feature extraction.BioCompoundML includes a copy of PaDEL-Discriptor, a molecular descriptor calculator, that takes as input an SDF file and provides thousands of individual QSPR and QSAR descriptors for each compound. By default, BioCompoundML calculates 1444 of these descriptors (1D/2D descriptors). This software is provided with its open source Apache 2.0 License.

Imputing Missing Data
---------------------
Imputing missing values is achieved using a two-step approach. The first step is to perform K-Nearest Neighbors (KNN) imputation. This process takes a distance matrix and imputes missing values using the KNN. The distance matrix in this tool is calculated using the Jaccard Distance/Tanimoto Score, using the 881 NCBI fingerprint variables. This allows the distance matrix to be collected separately from value imputation. This matrix is used to identify the nearest neighbors. The default for BioCompoundML is k=5. The distance matrix is then used to assign a weight to each value for the nearest neighbors and return a weighted average, such that nearer neighbors (in this case compounds) are more heavily weighted. This approach is generalizable and has shown consistent success as an approach to missing data. In cases where features were too sparse to fully resolve using KNN, we used the mean value for the feature as a minimum information imputer.

Feature Reduction
-----------------
The Boruta algorithm was chosen for selecting features for classification. This algorithm generates a set of shadow random features, duplicating and then shuffling the variables. The result of this is a set of Z-score distributions for each feature. Each original feature is compared to the maximum Z-score for the list of shadow features. Features that fail to score significantly better than distribution of shadow features (using standard t-tests) are then excluded from the model. This step can dramatically reduce the complexity of the model - eliminating needless and uninformative features (see https://github.com/danielhomola/boruta_py for examples of this).

Random Forest Classification
----------------------------
The main function of BioCompoundML is to run the Random Forest Classifier. This function ties up the scikit-learn RandomForestClassifier. The default parameters are n_estimators=512, oob_score=True and n_jobs=-1, which specify that the number of estimators be high, the out of bag samples is used to estimate the general error and it is run on as many cpus as possible.

Cross-Validation
----------------
By default BioCompoundML runs 50% leave-out cross-validation 100 times. This allows the calculation of the mean and standard deviations for accuracy, precision, recall and Receiver Operator Characteristic Area Under the Curve.

Feature Weightings
------------------
BioCompoundML also weights individual features that were used to build the model. If Boruta was selected using ``--selection``, this only includes the reduced features ::

	Complexity 0.1778
	XLogP3-AA 0.1751
	Rotatable Bond Count 0.1317
	Monoisotopic Mass 0.0771
	Molecular Weight 0.0671


Testing
-------
Using the ``--test`` command, BioCompoundML, takes a second file after the ``--test_input`` parameter in nearly the same format as the training input (minus the training feature). If user input was provided for training, the same input will need to be provided in this file ::

	#Name	PubChem
	isoamyl acetate	31276
	myrcene	31253
	eucalyptol	2758
	3-carene	26049

The output looks this ::

	isoamyl acetate	[0.033, 0.967]
	myrcene	[0.201, 0.799]
	eucalyptol	[0.084, 0.916]
	3-carene	[0.246, 0.754]

This specifies the compound name, its probability of classification below the threshold and its probability of classification above the threshold.

Indices and tables
==================

* :doc:`index`
* :doc:`script`
* :doc:`module`


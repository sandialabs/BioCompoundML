[![Build Status](https://travis-ci.org/sandialabs/BioCompoundML.svg?branch=master)](https://travis-ci.org/sandialabs/BioCompoundML)

# BioCompoundML

Rapidly screen a large number of compounds for fuel and chemical properties using machine learning. It's quick -- build in minutes, screen in seconds. It's clean -- cluster, predict, report and validate in a single interface. And it directly connects to the PubChem API and a variety of Quantitative Structure (Property and Activity) Relationship predictors (QSPR/QSAR). 

## Documentation

See documentation at http://sandialabs.github.io/BioCompoundML/

## Build

The most difficult part of the build is getting scikit-learn up and running and beautiful-soup. It is best that you use an existing tool, like conda or canopy or another scientific python distribution. If not, it may take some effort to get BioCompoundML running on your machine. Ultimately, you will need numpy, scipy, scikit-learn, matplotlib and beautiful-soup. If you have those the rest of the setup should be fairly painless.

```bash
git clone https://github.com/sandialabs/BioCompoundML.git
pip install -r requirements.txt
python setup.py install
```

### Dependencies
-------------
BioCompoundML is tested to work under Python 2.7.

* numpy==1.10.4
* scikit-learn==0.17.1
* scipy==0.17.0
* beautiful-soup==4.4.1
* matplotlib==1.4.0

## Deeper dive
============


## License

BSD - 3-Clause Copyright 2016 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.

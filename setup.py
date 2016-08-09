from setuptools import find_packages, setup

setup(name="BioCompoundML",
      version="1.0.0a1",
      description="A random forest classifier for biological compound screening",
      author=["Corey M. Hudson",
              "Leanne Whitmore"
              ],
      author_email=["coreymhudson@gmail.com",
                    "lwhitmo@sandia.gov"
                    ],
      platforms=["linux",
                 "osx"
                 ],
      license="BSD 3 clause",
      url="https://github.com/sandialabs/BioCompoundML",
      classifiers=['Development Status :: 3 - Alpha',
                   'Intended Audience :: Developers, AI Researchers, Cheminformaticists, Bioinformaticists',
		   'Topic :: Cheminformatics :: Machine Learning',
   		   'License :: BSD 3 clause',
                   'Programming Language :: Python :: 2',
                   'Programming Language :: Python :: 2.6',
                   'Programming Language :: Python :: 2.7',
		   ],
      keywords='cheminformatics, chemoinformatics, machine learning',
      test_suite="tests",
      packages=find_packages(),
      )

"""
This module collects and processes user data
"""


from __future__ import print_function
import sys


def dictitems(dict):
    if sys.version_info[0]>=3:
        return dict.items()
    else:
        return dictitems(dict)

def verbose_print(verbose, line):
    if verbose:
        print(line)


class TestingUpdate(object):
    def _update_compounds(self):
        for id, compound in dictitems(self.compound.copy()):
            binhash = compound['userhash']
            for key, value in dictitems(binhash.copy()):
                if key not in self.features:
                    self.compound[id]['userhash'].pop(key, None)

    def update(self):
        self._update_compounds()
        del self.features

    def __init__(self, compounds, features):
        self.compound = compounds.compound
        self.features = features


class Update(object):
    def _variable_features(self):
        total = dict()
        variable = dict()
        for id, compound in dictitems(self.compound):
            userhash = compound['userhash']
            for key, value in dictitems(userhash):
                try:
                    val = float(value)
                    if key not in total:
                        total[key] = val
                    elif val != total[key]:
                        variable[key] = val
                except:
                    pass
        self.variable = variable.keys()
        print_string = "Total user features " + str(len(total.keys()))
        verbose_print(self.verbose, print_string)
        print_string = "Variable user features " + str(len(variable.keys()))
        verbose_print(self.verbose, print_string)

    def _update_compounds(self):
        for id, compound in dictitems(self.compound.copy()):
            userhash = compound['userhash']
            for key, value in dictitems(userhash.copy()):
                if key not in self.variable:
                    self.compound[id]['userhash'].pop(key, None)

    def update(self):
        self._variable_features()
        if self.remove is True:
            self._update_compounds()
        del self.variable

    def __init__(self, compounds, verbose=False, remove_static=True):
        self.compound = compounds.compound
        self.verbose = verbose
        self.remove = remove_static

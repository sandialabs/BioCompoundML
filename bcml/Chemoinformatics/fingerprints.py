"""
This module processes fingerprint data from NCBI
"""


def verbose_print(verbose, line):
    if verbose:
        print line


class TestingUpdate(object):
    def _update_compounds(self):
        for id, compound in self.compound.copy().iteritems():
            binhash = compound['binhash']
            for key, value in binhash.copy().iteritems():
                if key not in self.features:
                    self.compound[id]['binhash'].pop(key, None)

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
        for id, compound in self.compound.iteritems():
            binhash = compound['binhash']
            for key, value in binhash.iteritems():
                if key not in total:
                    total[key] = value
                elif value != total[key]:
                    variable[key] = value
        self.variable = variable.keys()
        print_string = "Total fingerprint features " + str(len(total.keys()))
        verbose_print(self.verbose, print_string)
        print_string = "Variable fingerprint features " + str(len(variable.keys()))
        verbose_print(self.verbose, print_string)

    def _update_compounds(self):
        for id, compound in self.compound.copy().iteritems():
            binhash = compound['binhash']
            for key, value in binhash.copy().iteritems():
                if key not in self.variable:
                    self.compound[id]['binhash'].pop(key, None)

    def update(self):
        self._variable_features()
        if self.remove is True:
            self._update_compounds()
        del self.variable

    def __init__(self, compounds, verbose=False, remove_static=True):
        self.compound = compounds.compound
        self.verbose = verbose
        self.remove = remove_static

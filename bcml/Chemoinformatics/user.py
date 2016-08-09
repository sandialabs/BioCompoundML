"""
This module collects and processes user data
"""


def verbose_print(verbose, line):
    if verbose:
        print line


class TestingUpdate(object):
    def _update_compounds(self):
        for id, compound in self.compound.copy().iteritems():
            binhash = compound['userhash']
            for key, value in binhash.copy().iteritems():
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
        for id, compound in self.compound.iteritems():
            userhash = compound['userhash']
            for key, value in userhash.iteritems():
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
        for id, compound in self.compound.copy().iteritems():
            userhash = compound['userhash']
            for key, value in userhash.copy().iteritems():
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

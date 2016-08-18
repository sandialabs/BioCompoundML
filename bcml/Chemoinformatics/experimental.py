"""
This module parses experimental/computed data from NCBI
"""

from __future__ import print_function
import sys
import re
import os
_path = os.path.abspath(__file__)
_directory = os.path.dirname(_path)
_features_list = _directory + '/feature_list.txt'


def dictitems(dict):
    if sys.version_info[0]>=3:
        return dict.items()
    else:
        return dict.iteritems()

def verbose_print(verbose, line):
    if verbose:
        print(line)


def _add_feature_list():
    features = set()
    with open(_features_list) as filep:
        for line in filep:
            features.add(line.strip())
    return features


def _format(string):
    '''This function attempts to get at the
    ideosyncracies associated with formatting
    pubchem xml files. This is the likely
    location of errors in features. It is always
    worth double checking the values here. One
    of the main hacks deep in the code is that we
    select the first best experimental value, which
    is, by convention the best, but may be in
    different units or different pressures or at
    a different temperature than other variables.
    Errors here most commonly resolve to this
    problem
    '''
    string = string.replace(u'\xb0', ' deg ')
    string = string.replace('Deg ', 'deg ')
    string = string.replace(',', '')
    '''Convert F to C'''
    if string.find('def F'):
        F_string = re.compile("(-?\d+) deg F")
        match = F_string.match(string)
        if match:
            string = (float(match.group(1))-32.0)/1.8
            string = str(string)
    '''Convert amoung pressures'''
    if ' kPa' in string:
        hg_string = re.compile("(-?\d+) kPa")
        match = hg_string.match(string)
        if match:
            string = float(match.group(1))*7.50061561303
            string = str(string)
    if string.startswith('kPa'):
        (first, number) = string.split(':')
        string = float(number)*7.50061561303
        string = str(string)
    '''Handle scientific notation'''
    if 'X' in string:
        string = string.replace('X10+', 'e')
        string = string.replace('X10', 'e')
    elif '-' in string[1:]:
        range_string = re.compile("(\d+)-(\d+)")
        match = range_string.match(string)
        if match:
            string = (float(match.group(1))+float(match.group(2)))/2.0
            string = str(string)
    string_out = re.findall(r"-?\d+(?:\.\d*)?(?:[eE][+\-]?\d+)?", string)
    try:
        string = "{:.8f}".format(float(string))
    except:
        pass
    try:
        return float(string_out[0])
    except:
        return False


class TestingUpdate(object):
    def update(self):
        features = self.features
        for id, compound in dictitems(self.compound.copy()):
            if 'xml' in compound.keys():
                for key, value in dictitems(compound['xml'].copy()):
                    if key in features:
                        self.compound[id]['xml'][key] = _format(value)
                    else:
                        self.compound[id]['xml'].pop(key, None)
                self.compound[id]['experimentalhash'] = self.compound[id].pop('xml')
        del self.features

    def __init__(self, compounds, features):
        self.compound = compounds.compound
        self.features = features


class Update(object):
    def update(self):
        self._variable_features()
        if self.remove is True:
            self._update_compounds()
        del self.variable
        del self.total

    def _update_compounds(self):
        '''Remove the invariable compounds from the feature set'''
        for id, compound in dictitems(self.compound.copy()):
            binhash = compound['experimentalhash']
            for key, value in dictitems(binhash.copy()):
                if key not in self.variable:
                    self.compound[id]['experimentalhash'].pop(key, None)

    def _variable_features(self):
        '''Start by deciding which features to include, this is by default
        the 26 features stored in feature_list.txt'''
        features = _add_feature_list()
        self.total = {}
        self.variable = {}
        for _id, compound in dictitems(self.compound.copy()):
            if 'xml' in compound.keys():
                for key, value in dictitems(compound['xml'].copy()):
                    if key in features:
                        self.compound[_id]['xml'][key] = _format(value)
                        if key not in self.total:
                            self.total[key] = self.compound[_id]['xml'][key]
                        elif self.compound[_id]['xml'][key] != self.total[key]:
                            self.variable[key] = 1
                    else:
                        self.compound[_id]['xml'].pop(key, None)
                '''Remove the original xml data from the features and redefine
                it as experimentalhash'''
                self.compound[_id]['experimentalhash'] = self.compound[_id].pop('xml')
        print_string = "Total experimental features " + str(len(self.total.keys()))
        verbose_print(self.verbose, print_string)
        print_string = "Variable experimental features " + str(len(self.variable.keys()))
        verbose_print(self.verbose, print_string)

    def __init__(self, compounds, verbose=False, remove_static=True):
        self.compound = compounds.compound  # Import compounds
        self.verbose = verbose
        self.remove = remove_static  # If true then remove redundant features

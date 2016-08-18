'''
This is the workinng class for Chemofeature collection
'''


from __future__ import print_function
from tempfile import mkdtemp
from subprocess import call
import os
import sys
from shutil import rmtree


_path = os.path.abspath(__file__)
_directory = os.path.dirname(_path)
_fp_file = _directory + '/fingerprints.txt'
_db_directory = _directory + "/db/pubchem"
_testing_directory = _directory + "/db/testing/"
_padel_descriptor = _directory + "/padel_descriptor.out"


def dictitems(dict):
    if sys.version_info[0]>=3:
        return dict.items()
    else:
        return dictitems(dict)

def verbose_print(verbose, line):
    if verbose:
        print(line)


class TestingUpdate(object):
    def _parsePadel(self):
        features = self.features
        with open(_padel_descriptor, "r") as filebuffer:
            header = filebuffer.readline()
            keys = header.strip().split(',')
            for line in filebuffer:
                ld = {}
                linelist = line.strip().split(',')
                for count, var in enumerate(ld):
                    if keys[count] in features:
                        ld[keys[count]] = var.replace('"', '')
                self.compound[linelist[0].replace('"', '')]['padelhash'] = ld
        del self.features

    def update(self, padel=False):
        if padel is False:
            [self.compound[k].pop('sdf') for k, v in dictitems(self.compound.copy())]
        else:
            cid_list = self.compound.keys()
            for cid in cid_list:
                filename = _testing_directory+str(cid)+".sdf"
                if not os.path.isfile(filename) and 'sdf' in self.compound[cid].keys():
                    with open(filename, "w") as filebuffer:
                        filebuffer.write(self.compound[cid]['sdf'])
            [self.compound[k].pop('sdf') for k, v in dictitems(self.compound.copy())]
            call(["java", "-jar", _directory+"/PaDEL-Descriptor.jar", "-dir",
                 _testing_directory, "-2d", "-3d", "-file", _padel_descriptor,
                 "-threads", "-1"])
            for cid in cid_list:
                filename = _testing_directory + str(cid) + ".sdf"
                os.remove(filename)
            self._parsePadel()

    def __init__(self, compounds, features):
        self.compound = compounds.compound
        self.features = features


class Update(object):
    def _update_compounds(self):
        '''Remove invariable features'''
        for id, compound in dictitems(self.compound.copy()):
            for key in self.compound[id]['padelhash'].copy().keys():
                if key not in self.variable:
                    self.compound[id]['padelhash'].pop(key, None)

    def _parsePadel(self):
        '''Directly parse PaDEL-Descriptor data'''
        self.total = dict()
        self.variable = dict()
        with open(_padel_descriptor, "r") as filebuffer:
            header = filebuffer.readline()
            keys = header.strip().split(',')
            for line in filebuffer:
                linedict = {}
                linelist = line.strip().split(',')
                for count, var in enumerate(linelist):
                    var = var.replace('"', '')
                    if not var:
                        var = 0
                    linedict[keys[count]] = var
                    key = keys[count]
                    var = str(var).replace('"', '')
                    if key not in self.total:
                        self.total[key] = var
                    elif var != self.total[key]:
                        self.variable[key] = var
                self.compound[linelist[0].replace('"', '')]['padelhash'] = linedict
        print_string = "Total PaDEL-Descriptors " + str(len(self.total.keys()))
        verbose_print(self.verbose, print_string)
        print_string = "Variable PaDEL-Descriptors features " + str(len(self.variable.keys()))
        verbose_print(self.verbose, print_string)

    def _variable_compounds(self):
        if self.clean is True:  # Create a temporary directory for .sdfs
            dirpath = mkdtemp()
        else:
            dirpath = _db_directory
        if self.padel is False:
            '''Get rid of the sdf features'''
            [self.compound[k].pop('sdf') for k, v in dictitems(self.compound.copy())]
        else:
            '''PaDEL updating is more complicated than the other features, since
            it requires us to not only extract sdfs from PubChem, but to read them
            into a directory and run the java script'''
            cid_list = self.compound.keys()  # Use these to write sdfs
            for cid in cid_list:
                filename = dirpath + '/' + str(cid) + ".sdf"
                if not os.path.isfile(filename) and 'sdf' in self.compound[cid].keys():
                    with open(filename, "w") as filebuffer:
                        filebuffer.write(self.compound[cid]['sdf'])
            [self.compound[k].pop('sdf') for k, v in dictitems(self.compound.copy())]
            '''Specifies to run PaDEL-Descriptor, with max threads, using -2d and -3d
            descriptors'''
            call(["java", "-jar", _directory + "/PaDEL-Descriptor.jar", "-dir",
                 dirpath, "-2d", "-3d", "-file", _padel_descriptor,
                 "-threads", "-1"])
            self._parsePadel()
        if self.clean is True:
            rmtree(dirpath)  # Clean .sdfs from filesystem
            os.remove(_padel_descriptor)  # Clean _descriptors from filesystem

    def update(self, padel=True, clean=True):
        self.padel = padel  # Run PaDEL-Descriptor
        self.clean = clean  # Clean working directory
        self._variable_compounds()
        if self.remove is True:
            self._update_compounds()
        del self.variable
        del self.total

    def __init__(self, compounds, verbose=False, remove_static=True):
        self.compound = compounds.compound  # Import compounds
        self.verbose = verbose
        self.remove = remove_static  # If true then remove redundant features

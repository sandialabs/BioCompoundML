"""
This process takes in the testing dataset and outputs
a data structure that includes the name of the molecule
and the CAS number

Attributes:
    input_file (str): This is the training file that
    is read by the output
    Instance (class): This is a private class which
    structures each instance
    Model (class): This is a public class with the
    total structure of the set
"""


class Read(object):
    """This file reads a testing file"""
    def _read_file(self):
        header = {}
        compounds = []
        with open(self.input_file) as fb:
            line = fb.readline()
            if line.startswith('#'):
                line = line.replace("#", "")
                head = line.strip().split('\t')
                for count, item in enumerate(head):
                    header[count] = item
            else:
                header[0] = 'Name'
                header[2] = 'CAS'
            for line in fb:
                larray = line.strip().split('\t')
                compound = {}
                if self.user is True:
                    compound['userhash'] = {}
                for count, item in enumerate(larray):
                    if header[count] == 'Name':
                        compound['Name'] = item.rstrip()
                    if self.user is True:
                        if header[count] == self.id_name:
                            compound[self.id_name] = item.rstrip()
                        compound['userhash'][header[count]] = item.rstrip()
                    else:
                        compound[header[count]] = item
                compounds.append(compound)
        return compounds

    def __init__(self, input_file, user=False, id_name=False):
        self.input_file = input_file
        self.user = user
        self.id_name = id_name
        self.compounds = self._read_file()

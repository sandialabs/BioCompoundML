"""
This process takes in the training dataset and outputs
a data structure that includes the name of the molecule,
the predictor, and the CAS number

Attributes:
    input_file (str): This is the training file that
    is read by the output
    Instance (class): This is a private class which
    structures each instance
    Model (class): This is a public class with the
    total structure of the set
"""


class Read(object):
    """This file reads a training file"""
    def _read_file(self):
        header = {}
        compounds = []
        predictors = []
        predictor = False
        with open(self.input_file) as fb:
            line = fb.readline()
            if line.startswith('#'):
                line = line.replace("#", "")
                head = line.strip().split('\t')
                for count, item in enumerate(head):
                    header[count] = item.rstrip()
            else:
                header[0] = 'Name'
                header[1] = 'Predictor'
                self.predictor = header[1]
                header[2] = 'CAS'
            for line in fb:
                larray = line.strip().split('\t')
                compound = {}
                if self.user:
                    compound['userhash'] = {}
                for count, item in enumerate(larray):
                    if header[count] == self.predictor:
                        predictor = item
                    elif self.user is True:
                        compound['userhash'][header[count]] = item.rstrip()
                    compound[header[count]] = item
                compounds.append(compound)
                predictors.append(predictor)
        return (compounds, predictors)

    def __init__(self, input_file, predictor=False, user=False, id_name=False):
        self.input_file = input_file
        self.predictor = predictor
        self.user = user
        self.id_name = id_name
        (self.compounds, self.predictors) = self._read_file()

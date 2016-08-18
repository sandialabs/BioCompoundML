"""

This contains the unit tests for the read_training module.

"""


from __future__ import print_function
import unittest
from Chemoinformatics import experimental as exp
from collections import defaultdict
from copy import deepcopy


class Object(object):
    pass


class ExperimentalTest(unittest.TestCase):

    def setUp(self):
        """Create an instance of the Update class"""
        print("Initializing test")
        compound = defaultdict(dict)
        compound = {'7844': {'xml': {u'Density': u'0.577 at 25 deg C/4 deg C', u'Heat of Combustion': u'-2719.1 KJ/mol at constant pressure and 25 deg C', u'Odor Threshold': u'69 ppb', u'Vapor Density': u'1.93 (Air= 1)', u'Decomposition': u'When heated to decomposition it emits acrid smoke and fumes.', u'Boiling Point': u'-6.47 deg C at 760 mm Hg', u'Hydrogen Bond Donor Count': u'0', u'Rotatable Bond Count': u'1', u'XLogP3': u'2.4', u'Melting Point': u'-185.3 deg C (FP)', u'Odor': u'Slightly aromatic odor', u'Flash Point': u'Flammable gas', u'Viscosity': u'7.76X10-3 mPa sec of saturated vapor at 298.15K; 0.186 mPa sec of saturated liquid at 266 K', u'Undefined Atom Stereocenter Count': u'0', u'Solubility': u'Sol in alcohol, ether, <a class="pubchem-internal-link CID-241" href="https://pubchem.ncbi.nlm.nih.gov/compound/benzene">benzene</a>', u'Formal Charge': u'0', u'Library': [u'Main library', u'Replicate library', u'Replicate library'], u'Auto-Ignition': u'725 DEG F (385 DEG C)', u'Vapour pressure': u'kPa at 21\xb0C: 464', u'Molecular Weight': u'56.10632', u'Hydrogen Bond Acceptor Count': u'0', u'LogP': u'log Kow= 2.40', u'Defined Atom Stereocenter Count': u'0', u'Surface Tension': u'0.0121 dyn/cm at 298.15 K', u'Complexity': u'14', u'Vapor Pressure': u'2.253X10+3 mm Hg at 25 deg C', u'Covalently-Bonded Unit Count': u'1', u'Isotope Atom Count': u'0', u'Spectral Properties': u'MAX ABSORPTION (GAS): 162 NM SHOULDER (LOG E= 4.0); 175 NM (LOG E= 4.2); 187 NM (LOG E= 4.1)', u'Heavy Atom Count': u'4', u'Physical Description': u'ODOURLESS COLOURLESS COMPRESSED LIQUEFIED GAS.', u'Color': u'Colorless', u'Exact Mass': u'56.0626', u'Monoisotopic Mass': u'56.0626', u'Heat of Vaporization': u'20.22 KJ/mol at 25 deg C', u'Stability': u'Stable', u'Topological Polar Surface Area': u'0', u'Undefined Bond Stereocenter Count': u'0', u'Defined Bond Stereocenter Count': u'0'}}, '19502': {'xml': {u'Library': u'Main library', u'Hydrogen Bond Donor Count': u'0', u'Rotatable Bond Count': u'1', u'XLogP3-AA': u'3.8', u'Heavy Atom Count': u'8', u'm/z Top Peak': u'55', u'Undefined Atom Stereocenter Count': u'3', u'Formal Charge': u'0', u'Molecular Weight': u'112.21264', u'Hydrogen Bond Acceptor Count': u'0', u'Total Peaks': u'74', u'NIST Number': u'60953', u'Defined Atom Stereocenter Count': u'0', u'm/z 2nd Highest': u'83', u'Complexity': u'66.4', u'Covalently-Bonded Unit Count': u'1', u'Isotope Atom Count': u'0', u'm/z 3rd Highest': u'41', u'Exact Mass': u'112.125201', u'Monoisotopic Mass': u'112.125201', u'Topological Polar Surface Area': u'0', u'Undefined Bond Stereocenter Count': u'0', u'Defined Bond Stereocenter Count': u'0'}}, '11610': {'xml': {u'Density': u'0.6970 @ 20 deg C/4 deg C', u'Vapor Density': u'0.7 (AIR= 1)', u'Decomposition': u'When heated to decomposition it emits acrid smoke and fumes.', u'Boiling Point': u'93.6 deg C', u'Hydrogen Bond Donor Count': u'0', u'Rotatable Bond Count': u'4', u'XLogP3': u'4', u'Melting Point': u'-119.7 deg C', u'Flash Point': u'LESS THAN 32 DEG F (CLOSED CUP)', u'Viscosity': u'0.5 sq mm/s', u'Undefined Atom Stereocenter Count': u'0', u'Solubility': u'Soluble in <a class="pubchem-internal-link CID-702" href="https://pubchem.ncbi.nlm.nih.gov/compound/ethanol">ethanol</a> and ether; slightly soluble in <a class="pubchem-internal-link CID-5943" href="https://pubchem.ncbi.nlm.nih.gov/compound/carbon%20tetrachloride">carbon tetrachloride</a>', u'Formal Charge': u'0', u'Library': [u'Main library', u'Replicate library', u'Replicate library', u'Replicate library'], u'Auto-Ignition': u'500 DEG F', u'Molecular Weight': u'98.18606', u'Hydrogen Bond Acceptor Count': u'0', u'LogP': u'log Kow= 3.99', u'Defined Atom Stereocenter Count': u'0', u'Complexity': u'37.3', u'Vapor Pressure': u'59.3 mm Hg @ 25 deg C', u'Covalently-Bonded Unit Count': u'1', u'Isotope Atom Count': u'0', u'Spectral Properties': u'Index of refraction: 1.3994 @ 20 deg C', u'Heavy Atom Count': u'7', u'Physical Description': u'<a class="pubchem-internal-link CID-11610" href="https://pubchem.ncbi.nlm.nih.gov/compound/N-heptene">N-heptene</a> is a colorless liquid. Flash point 25\xb0F. Insoluble in <a class="pubchem-internal-link CID-962" href="https://pubchem.ncbi.nlm.nih.gov/compound/water">water</a> and much less dense than <a class="pubchem-internal-link CID-962" href="https://pubchem.ncbi.nlm.nih.gov/compound/water">water</a>. Vapors heavier than air. Used to make other chemicals.', u'Color': u'Colorless liquid', u'Exact Mass': u'98.10955', u'Monoisotopic Mass': u'98.10955', u'Topological Polar Surface Area': u'0', u'Undefined Bond Stereocenter Count': u'0', u'Defined Bond Stereocenter Count': u'0'}}}
        self.test_data = Object()
        self.test_data.compound = deepcopy(compound)
        self.original = Object()
        self.original.compound = deepcopy(compound)

    def tearDown(self):
        """Delete data structure"""
        print("Clearing out file")
        del self.test_data

    def test_fingerprint_functions(self):
        """Test the values associated with the test file"""
        print("Run Fingerprint functions")
        self.test_data = exp.Update(self.test_data, remove_static=True, verbose=True)
        self.test_data.update()
        self.assertTrue('7844' in self.test_data.compound)
        self.assertTrue('experimentalhash' in self.test_data.compound['7844'])
        self.assertTrue('Density' in self.test_data.compound['7844']['experimentalhash'])
        self.assertEquals(0.577, self.test_data.compound['7844']['experimentalhash']['Density'])
        variable_length = len(self.test_data.compound['7844']['experimentalhash'].keys())
        self.assertEquals(variable_length, 16)
        self.original = exp.Update(self.original, remove_static=False, verbose=True)
        self.original.update()
        total_length = len(self.original.compound['7844']['experimentalhash'].keys())
        self.assertEquals(total_length, 25)

if __name__ == '__main__':
    unittest.main()

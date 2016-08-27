"""

This contains the unit tests for the read_training module.

"""


from __future__ import print_function
import unittest
from Chemoinformatics import chemofeatures as cf
try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict
from copy import deepcopy


class Object(object):
    pass


class ChemoFeaturesTest(unittest.TestCase):

    def setUp(self):
        """Create an instance of the Update class"""
        print("Initializing test")
        compound = OrderedDict({'7844': {'sdf': '7844\n  -OEChem-07251617462D\n\n 12 11  0     0  0  0  0  0  0999 V2000\n    2.8660   -0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.7320    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.5981   -0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.4675   -0.7249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    3.2646   -0.7249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.3100    0.7869    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.4631    0.5600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.6900   -0.2869    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    3.7320    0.8700    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.5981   -0.8700    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.1350    0.0600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  1  3  1  0  0  0  0\n  1  5  1  0  0  0  0\n  1  6  1  0  0  0  0\n  2  7  1  0  0  0  0\n  2  8  1  0  0  0  0\n  2  9  1  0  0  0  0\n  3  4  2  0  0  0  0\n  3 10  1  0  0  0  0\n  4 11  1  0  0  0  0\n  4 12  1  0  0  0  0\nM  END\n> <PUBCHEM_COMPOUND_CID>\n7844\n\n> <PUBCHEM_COMPOUND_CANONICALIZED>\n1\n\n> <PUBCHEM_CACTVS_COMPLEXITY>\n14\n\n> <PUBCHEM_CACTVS_HBOND_ACCEPTOR>\n0\n\n> <PUBCHEM_CACTVS_HBOND_DONOR>\n0\n\n> <PUBCHEM_CACTVS_ROTATABLE_BOND>\n1\n\n> <PUBCHEM_CACTVS_SUBSKEYS>\nAAADccBgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAAAAAAACACAAAACAAAAAACAACBCAAAAAAAgAAAIAAAAAAgAAAAAAQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA==\n\n> <PUBCHEM_IUPAC_OPENEYE_NAME>\nbut-1-ene\n\n> <PUBCHEM_IUPAC_CAS_NAME>\n1-butene\n\n> <PUBCHEM_IUPAC_NAME>\nbut-1-ene\n\n> <PUBCHEM_IUPAC_SYSTEMATIC_NAME>\nbut-1-ene\n\n> <PUBCHEM_IUPAC_TRADITIONAL_NAME>\nbut-1-ene\n\n> <PUBCHEM_IUPAC_INCHI>\nInChI=1S/C4H8/c1-3-4-2/h3H,1,4H2,2H3\n\n> <PUBCHEM_IUPAC_INCHIKEY>\nVXNZUUAINFGPBY-UHFFFAOYSA-N\n\n> <PUBCHEM_XLOGP3>\n2.4\n\n> <PUBCHEM_EXACT_MASS>\n56.0626\n\n> <PUBCHEM_MOLECULAR_FORMULA>\nC4H8\n\n> <PUBCHEM_MOLECULAR_WEIGHT>\n56.10632\n\n> <PUBCHEM_OPENEYE_CAN_SMILES>\nCCC=C\n\n> <PUBCHEM_OPENEYE_ISO_SMILES>\nCCC=C\n\n> <PUBCHEM_CACTVS_TPSA>\n0\n\n> <PUBCHEM_MONOISOTOPIC_WEIGHT>\n56.0626\n\n> <PUBCHEM_TOTAL_CHARGE>\n0\n\n> <PUBCHEM_HEAVY_ATOM_COUNT>\n4\n\n> <PUBCHEM_ATOM_DEF_STEREO_COUNT>\n0\n\n> <PUBCHEM_ATOM_UDEF_STEREO_COUNT>\n0\n\n> <PUBCHEM_BOND_DEF_STEREO_COUNT>\n0\n\n> <PUBCHEM_BOND_UDEF_STEREO_COUNT>\n0\n\n> <PUBCHEM_ISOTOPIC_ATOM_COUNT>\n0\n\n> <PUBCHEM_COMPONENT_COUNT>\n1\n\n> <PUBCHEM_CACTVS_TAUTO_COUNT>\n1\n\n> <PUBCHEM_COORDINATE_TYPE>\n1\n5\n255\n\n$$$$'}, '19502': {'sdf': '19502\n  -OEChem-07251617462D\n\n 24 24  0     1  0  0  0  0  0999 V2000\n    3.0878    0.4239    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0\n    2.5878   -1.1149    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0\n    2.2788   -0.1639    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.8968   -0.1639    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.5878   -1.1149    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.0878    1.4239    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000   -1.9239    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.9538    1.9239    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.8451    0.8098    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.7482   -0.9819    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.9688    0.3731    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.7124   -0.4160    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.6733   -0.5096    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.2068    0.3731    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    3.5230   -1.7315    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.1942   -1.2438    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.8757    2.0065    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.4772    1.3163    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.4984   -1.5595    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.6356   -2.4255    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.5016   -2.2884    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.2638    1.3870    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.4908    2.2339    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    3.6438    2.4609    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  3  1  0  0  0  0\n  1  4  1  0  0  0  0\n  1  6  1  0  0  0  0\n  1  9  1  0  0  0  0\n  2  3  1  0  0  0  0\n  2  5  1  0  0  0  0\n  2  7  1  0  0  0  0\n  2 10  1  0  0  0  0\n  3 11  1  0  0  0  0\n  3 12  1  0  0  0  0\n  4  5  1  0  0  0  0\n  4 13  1  0  0  0  0\n  4 14  1  0  0  0  0\n  5 15  1  0  0  0  0\n  5 16  1  0  0  0  0\n  6  8  1  0  0  0  0\n  6 17  1  0  0  0  0\n  6 18  1  0  0  0  0\n  7 19  1  0  0  0  0\n  7 20  1  0  0  0  0\n  7 21  1  0  0  0  0\n  8 22  1  0  0  0  0\n  8 23  1  0  0  0  0\n  8 24  1  0  0  0  0\nM  END\n> <PUBCHEM_COMPOUND_CID>\n19502\n\n> <PUBCHEM_COMPOUND_CANONICALIZED>\n1\n\n> <PUBCHEM_CACTVS_COMPLEXITY>\n66.4\n\n> <PUBCHEM_CACTVS_HBOND_ACCEPTOR>\n0\n\n> <PUBCHEM_CACTVS_HBOND_DONOR>\n0\n\n> <PUBCHEM_CACTVS_ROTATABLE_BOND>\n1\n\n> <PUBCHEM_CACTVS_SUBSKEYS>\nAAADceBwAAAAAAAAAAAAAAAAAAAAAYAAAAAAAAAAAAAAAAAAAAAAGAAAAAAADQCAAAACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIAAAAAAAAAAAAAAAEAgEAOAAAAAAAAAAAAAAAAAAAAAQAAAAAAAA==\n\n> <PUBCHEM_IUPAC_OPENEYE_NAME>\n1-ethyl-3-methyl-cyclopentane\n\n> <PUBCHEM_IUPAC_CAS_NAME>\n1-ethyl-3-methylcyclopentane\n\n> <PUBCHEM_IUPAC_NAME>\n1-ethyl-3-methylcyclopentane\n\n> <PUBCHEM_IUPAC_SYSTEMATIC_NAME>\n1-ethyl-3-methyl-cyclopentane\n\n> <PUBCHEM_IUPAC_TRADITIONAL_NAME>\n1-ethyl-3-methyl-cyclopentane\n\n> <PUBCHEM_IUPAC_INCHI>\nInChI=1S/C8H16/c1-3-8-5-4-7(2)6-8/h7-8H,3-6H2,1-2H3\n\n> <PUBCHEM_IUPAC_INCHIKEY>\nPQXAPVOKLYINEI-UHFFFAOYSA-N\n\n> <PUBCHEM_XLOGP3_AA>\n3.8\n\n> <PUBCHEM_EXACT_MASS>\n112.125201\n\n> <PUBCHEM_MOLECULAR_FORMULA>\nC8H16\n\n> <PUBCHEM_MOLECULAR_WEIGHT>\n112.21264\n\n> <PUBCHEM_OPENEYE_CAN_SMILES>\nCCC1CCC(C1)C\n\n> <PUBCHEM_OPENEYE_ISO_SMILES>\nCCC1CCC(C1)C\n\n> <PUBCHEM_CACTVS_TPSA>\n0\n\n> <PUBCHEM_MONOISOTOPIC_WEIGHT>\n112.125201\n\n> <PUBCHEM_TOTAL_CHARGE>\n0\n\n> <PUBCHEM_HEAVY_ATOM_COUNT>\n8\n\n> <PUBCHEM_ATOM_DEF_STEREO_COUNT>\n0\n\n> <PUBCHEM_ATOM_UDEF_STEREO_COUNT>\n3\n\n> <PUBCHEM_BOND_DEF_STEREO_COUNT>\n0\n\n> <PUBCHEM_BOND_UDEF_STEREO_COUNT>\n0\n\n> <PUBCHEM_ISOTOPIC_ATOM_COUNT>\n0\n\n> <PUBCHEM_COMPONENT_COUNT>\n1\n\n> <PUBCHEM_CACTVS_TAUTO_COUNT>\n1\n\n> <PUBCHEM_COORDINATE_TYPE>\n1\n5\n255\n\n> <PUBCHEM_BONDANNOTATIONS>\n1  9  3\n2  10  3\n4  13  3\n\n$$$$'}, '11610': {'sdf': '11610\n  -OEChem-07251617462D\n\n 21 20  0     0  0  0  0  0  0999 V2000\n    5.4641    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.5981   -0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.3301   -0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.7320    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    7.1962    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.8660   -0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    5.8626    0.7249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.0656    0.7249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.1996   -0.7249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.9966   -0.7249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.9316   -0.7249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    6.7287   -0.7249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.1306    0.7249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    3.3335    0.7249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    7.5062   -0.2869    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    7.7331    0.5600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    6.8862    0.7869    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.8660   -0.8700    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000    0.8700    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.4631   -0.0600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  1  3  1  0  0  0  0\n  1  8  1  0  0  0  0\n  1  9  1  0  0  0  0\n  2  4  1  0  0  0  0\n  2 10  1  0  0  0  0\n  2 11  1  0  0  0  0\n  3  5  1  0  0  0  0\n  3 12  1  0  0  0  0\n  3 13  1  0  0  0  0\n  4  6  1  0  0  0  0\n  4 14  1  0  0  0  0\n  4 15  1  0  0  0  0\n  5 16  1  0  0  0  0\n  5 17  1  0  0  0  0\n  5 18  1  0  0  0  0\n  6  7  2  0  0  0  0\n  6 19  1  0  0  0  0\n  7 20  1  0  0  0  0\n  7 21  1  0  0  0  0\nM  END\n> <PUBCHEM_COMPOUND_CID>\n11610\n\n> <PUBCHEM_COMPOUND_CANONICALIZED>\n1\n\n> <PUBCHEM_CACTVS_COMPLEXITY>\n37.3\n\n> <PUBCHEM_CACTVS_HBOND_ACCEPTOR>\n0\n\n> <PUBCHEM_CACTVS_HBOND_DONOR>\n0\n\n> <PUBCHEM_CACTVS_ROTATABLE_BOND>\n4\n\n> <PUBCHEM_CACTVS_SUBSKEYS>\nAAADccBgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAAAAAAACACAAAACAAAAAACAACBCAAAAAAAgAAAIAAAAAAgAAAIAAQAAAAAAgAAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA==\n\n> <PUBCHEM_IUPAC_OPENEYE_NAME>\nhept-1-ene\n\n> <PUBCHEM_IUPAC_CAS_NAME>\n1-heptene\n\n> <PUBCHEM_IUPAC_NAME>\nhept-1-ene\n\n> <PUBCHEM_IUPAC_SYSTEMATIC_NAME>\nhept-1-ene\n\n> <PUBCHEM_IUPAC_TRADITIONAL_NAME>\nhept-1-ene\n\n> <PUBCHEM_IUPAC_INCHI>\nInChI=1S/C7H14/c1-3-5-7-6-4-2/h3H,1,4-7H2,2H3\n\n> <PUBCHEM_IUPAC_INCHIKEY>\nZGEGCLOFRBLKSE-UHFFFAOYSA-N\n\n> <PUBCHEM_XLOGP3>\n4\n\n> <PUBCHEM_EXACT_MASS>\n98.10955\n\n> <PUBCHEM_MOLECULAR_FORMULA>\nC7H14\n\n> <PUBCHEM_MOLECULAR_WEIGHT>\n98.18606\n\n> <PUBCHEM_OPENEYE_CAN_SMILES>\nCCCCCC=C\n\n> <PUBCHEM_OPENEYE_ISO_SMILES>\nCCCCCC=C\n\n> <PUBCHEM_CACTVS_TPSA>\n0\n\n> <PUBCHEM_MONOISOTOPIC_WEIGHT>\n98.10955\n\n> <PUBCHEM_TOTAL_CHARGE>\n0\n\n> <PUBCHEM_HEAVY_ATOM_COUNT>\n7\n\n> <PUBCHEM_ATOM_DEF_STEREO_COUNT>\n0\n\n> <PUBCHEM_ATOM_UDEF_STEREO_COUNT>\n0\n\n> <PUBCHEM_BOND_DEF_STEREO_COUNT>\n0\n\n> <PUBCHEM_BOND_UDEF_STEREO_COUNT>\n0\n\n> <PUBCHEM_ISOTOPIC_ATOM_COUNT>\n0\n\n> <PUBCHEM_COMPONENT_COUNT>\n1\n\n> <PUBCHEM_CACTVS_TAUTO_COUNT>\n1\n\n> <PUBCHEM_COORDINATE_TYPE>\n1\n5\n255\n\n$$$$'}})
        self.test_data = Object()
        self.test_data.compound = deepcopy(compound)
        self.original = Object()
        self.original.compound = deepcopy(compound)

    def tearDown(self):
        """Delete data structure"""
        print("Clearing out file")
        del self.test_data
        del self.original

    def test_chemoinfo_functions(self):
        """Test the values associated with the test file"""
        print("Run Chemoinformatics functions")
        self.test_data = cf.Update(self.test_data, remove_static=True, verbose=True)
        self.test_data.update()
        self.assertTrue('7844' in self.test_data.compound)
        self.assertTrue('padelhash' in self.test_data.compound['7844'])
        self.assertTrue('MATS3v' in self.test_data.compound['7844']['padelhash'])
        self.assertAlmostEquals(-0.05263157, float(self.test_data.compound['7844']['padelhash']['MATS3v']))
        variable_length = len(self.test_data.compound['7844']['padelhash'].keys())
        self.assertEquals(variable_length, 843)
        self.original = cf.Update(self.original, remove_static=False, verbose=True)
        self.original.update()
        total_length = len(self.original.compound['7844']['padelhash'].keys())
        self.assertEquals(total_length, 1876)

if __name__ == '__main__':
    unittest.main()

from PubChemUtils import pubchempy_utils as pcp
import unittest


class PubChemPyUtilException(Exception):

    """Exception type for the Pubchem Class"""

    pass


class PubChemPyUtilTests(unittest.TestCase):
    """Test the Pubchem class."""

    def setUp(self):
        """Initialize before every test."""
        c = [{'userhash': {'PubChem': '7844', 'Name': '1-Butene'}, 'Name': '1-Butene', 'RON': '98.8', 'PubChem': '7844'},
             {'userhash': {'PubChem': '19502', 'Name': '1-Ethyl-3-Methylcyclopentane'}, 'Name': '1-Ethyl-3-Methylcyclopentane', 'RON': '57.6', 'PubChem': '19502'},
             {'userhash': {'PubChem': '11610', 'Name': '1-Heptene'}, 'Name': '1-Heptene', 'RON': '54.5', 'PubChem': '11610'},
             {'userhash': {'PubChem': '11597', 'Name': '1-Hexene'}, 'Name': '1-Hexene', 'RON': '76.4', 'PubChem': '11597'},
             {'userhash': {'PubChem': '7459', 'Name': '1-Isopropyl-4-methylcyclohexane'}, 'Name': '1-Isopropyl-4-methylcyclohexane', 'RON': '67.3', 'PubChem': '7459'},
             {'userhash': {'PubChem': '35411', 'Name': '1-Methyl-1-ethylcyclohexane'}, 'Name': '1-Methyl-1-ethylcyclohexane', 'RON': '68.7', 'PubChem': '35411'},
             {'userhash': {'PubChem': '136729', 'Name': '1-Methyl-2-ethylcyclopentane'}, 'Name': '1-Methyl-2-ethylcyclopentane', 'RON': '57.6', 'PubChem': '136729'},
             {'userhash': {'PubChem': '107252', 'Name': '1-Methyl-2-propylcyclohexane'}, 'Name': '1-Methyl-2-propylcyclohexane', 'RON': '29.9', 'PubChem': '107252'},
             {'userhash': {'PubChem': '8125', 'Name': '1-Octene'}, 'Name': '1-Octene', 'RON': '28.7', 'PubChem': '8125'},
             {'userhash': {'PubChem': '8004', 'Name': '1-Pentene'}, 'Name': '1-Pentene', 'RON': '87.9', 'PubChem': '8004'},
             {'userhash': {'PubChem': '11549', 'Name': '1,1-Dimethylcyclohexane'}, 'Name': '1,1-Dimethylcyclohexane', 'RON': '87.3', 'PubChem': '11549'}]
        self._collect = pcp.Collect(c, verbose=True)

    def tearDown(self):
        """Clean up after each test."""
        del self._collect

    def test_set_proxy(self):
        """Ensure that a proxy can be established"""
        self._collect.set_proxy('http://www.example.com:80')

    def test_add_user(self):
        """Ensure that defined userhash can be created"""
        self._collect.add_user(user=True)

    def test_add_fingerprint(self):
        """Ensure that fingerprints can be added"""
        """Check that they can be added one at a time"""
        self._collect.add_fingerprint(fingerprint=True, chunks=1)
        self.assertTrue('binhash' in self._collect.compound['7844'])
        self.assertTrue('C-C-C=C' in self._collect.compound['7844']['binhash'])
        """Check that they can be added in bulk"""
        self._collect.add_fingerprint(fingerprint=True)
        self.assertTrue('binhash' in self._collect.compound['7844'])
        self.assertTrue('C-C-C=C' in self._collect.compound['7844']['binhash'])
        self.assertEqual('1', self._collect.compound['7844']['binhash']['C-C-C=C'])

    def test_add_sdf(self):
        """Ensure that SDFs can be added"""
        self._collect.add_sdf(sdf=True, chunks=False)
        self.assertTrue('sdf' in self._collect.compound['7844'])
        self.assertEqual('7844', self._collect.compound['7844']['sdf'].split('\n')[0])

    '''def test_add_xml(self):
        """Ensure that XMLs can be added"""
        self._collect.add_xml(xml=True)
        self.assertTrue('xml' in self._collect.compound['7844'])
        self.assertEqual(u'-6.47 deg C at 760 mm Hg', self._collect.compound['7844']['xml']['Boiling Point'])

    def test_pubchem(self):
        """Ensure that features were extracted."""
        self._convert.pubchem(5090)
        self.assertEqual('314.06128', self._convert.values['Exact Mass'])

    def test_pubchem_xrefs(self):
        """Ensure that links were extracted"""
        self._convert.pubchem_xrefs(5090)
        _found = 0
        _test_url = 'http://www.chemspider.com/Chemical-Structure.4911.html'
        if _test_url in set(self._convert.urls):
            _found = 1
        self.assertEqual(_found, 1)

    def test_pubchem_inchi(self):
        """Ensure that you can get cid from InChi"""
        inchi = 'InChI=1S/C3H8/c1-3-2/h3H2,1-2H3'
        cid = self._convert.pubchem_inchi(inchi)
        self.assertEqual(cid, '6334')

    def test_pubchem_synonym(self):
        """Ensure that cas number conversion works"""
        self._convert.synonym('162011-90-7')
        self.assertEqual(self._convert.cid, '5090')


    def test_pubchem_cactvs(self):
        _test_cactvs = 'AAADccB4OABAAAAAAAAAAAAAAAAAAQAAAAAwYAAAAAAAAAABQAAAG'\
                       'gQAAAAADACg2AKyCYAABAqIAiDSCHBCAAAgCBAIiBkAAMgIJDKgNR'\
                       'CCMAAkwAEIqQeIyKCOEAAAAAAAAAAgAAAAAAAAAAAAAAAAAA=='
        _test_binstring = "11000000011110000011100000000000010000000000000000"\
                          "00000000000000000000000000000000000000000000000000"\
                          "00000000000000000000000000000000000000000001000000"\
                          "00000000000000000000000000001100000110000000000000"\
                          "00000000000000000000000000000000000000000000000000"\
                          "00000101000000000000000000000000011010000001000000"\
                          "00000000000000000000000000000000110000000000101000"\
                          "00110110000000001010110010000010011000000000000000"\
                          "00000100000010101000100000000010001000001101001000"\
                          "00100001110000010000100000000000000000001000000000"\
                          "10000001000000001000100010000001100100000000000000"\
                          "00110010000000100000100100001100101010000000110101"\
                          "00010000100000100011000000000000001001001100000000"\
                          "00000100001000101010010000011110001000110010001010"\
                          "00001000111000010000000000000000000000000000000000"\
                          "00000000000000000000000000001000000000000000000000"\
                          "00000000000000000000000000000000000000000000000000"\
                          "0000000000000000000000000000000"
        _test_binstring = ''.join(_test_binstring)
        self._convert.pubchem_cactvs(5090)
        self.assertEqual(self._convert.binstring, _test_binstring)
        self.assertEqual(self._convert.cactvs, _test_cactvs)
    '''

if __name__ == '__main__':
    unittest.main()

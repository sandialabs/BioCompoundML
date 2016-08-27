"""

This contains the unit tests for the main bcml module.

"""


from __future__ import print_function
import unittest
import subprocess
import os


_dir = os.path.dirname(__file__)
_executable = os.path.abspath(os.path.join(_dir, os.pardir, 'bcml.py'))
_data_min = os.path.abspath(os.path.join(_dir, 'data', 'RON_min.txt'))
_data_mini = os.path.abspath(os.path.join(_dir, 'data', 'RON_mini.txt'))


class BCMLTests(unittest.TestCase):

    def setUp(self):
        """Create an instance of main"""
        print("Initializing tests")
        pass

    def tearDown(self):
        """Delete data structure"""
        print("Clearing out test suite")
        pass

    def test_parse_argument(self):
        print("Testing input parameters")
        print("Testing empty run")
        args = ['python', _executable]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 2)
        print("Testing read of incompatible modules")
        args = ['python', _executable, '-i data/null.txt', '--train', '--test']
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 2)
        args = ['python', _executable, '-i data/null.txt', '-d data/null.txt', '--pred', 'RON']
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 2)
        print("Testing random seeds")
        args = ['python', _executable, '-i', _data_min, '--train', '--pred', 'RON', '--user', '--random', '12345']
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)
        args = ['python', _executable, '-i', _data_min,  '--train', '--random', '-12345', '--pred', 'RON']
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 1)
        print("Testing proxy")
        args = ['python', _executable, '-i', _data_min, '--train', '--proxy', 'https://wwwproxy.example.com', '--pred', 'RON', '--user', '--random', '12345']
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)
        args = ['python', _executable, '-i', _data_min, '--train', '--proxy', '500', '--pred', 'RON', '--user', '--random', '12345']
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 1)
        print("Testing user feature")
        args = ['python', _executable, '-i', _data_min, '--train', '--pred', 'RON', '--user', '--random', '12345']
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)
        print("Testing fingerprint feature")
        args = ['python', _executable, '-i', _data_min, '--train', '--pred', 'RON', '--fingerprint', '--random', '12345']
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)
        print("Testing experimental feature")
        args = ['python', _executable, '-i', _data_min, '--train', '--pred', 'RON', '--verbose', '--experimental', '--random', '12345', '--split_value', '85.0']
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)
        print("Testing chemofeatures feature")
        args = ['python', _executable, '-i', _data_min, '--train', '--pred', 'RON', '--chemofeatures', '--random', '12345']
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)
        print("Testing distance feature")
        args = ['python', _executable, '-i', _data_min, '--train', '--pred', 'RON', '--fingerprint', '--distance', '--random', '12345']
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)
        print("Testing split feature")
        args = ['python', _executable, '-i', _data_min, '--train', '--pred', 'RON', '--fingerprint', '--random', '12345', '--split_value', '85.0']
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)
        print("Testing imputation")
        args = ['python', _executable, '-i', _data_min, '--train', '--pred', 'RON', '--fingerprint', '--random', '12345', '--impute']
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)
        print("Testing feature selection")
        args = ['python', _executable, '-i', _data_min, '--train', '--pred', 'RON', '--fingerprint', '--experimental', '--random', '12345', '--selection', '--split_value', '94.4']
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        print(stdoutdata, stderrdata)
        self.assertEqual(process.returncode, 0)
        print("Testing Cross Validation")
        args = ['python', _executable, '-i', _data_min, '--train', '--pred', 'RON', '--experimental', '--split_value', '94.4', '--random', '12345', '--impute', '--cv']
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdoutdata, stderrdata = process.communicate()
        self.assertEqual(process.returncode, 0)

if __name__ == '__main__':
    unittest.main()

from secnet import *
import unittest 

### lone file test 
"""
python3 -m tests.test_secnet
"""
###


class SecNetClass(unittest.TestCase):

    def test__SecNet__generate__case1(self):
        sn = SecNet_sample1()
        assert True

if __name__ == '__main__':
    unittest.main()
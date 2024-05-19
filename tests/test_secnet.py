from secnet import *
import unittest 

### lone file test 
"""
python3 -m tests.test_secnet
"""
###

class IsoRingedChainClass(unittest.TestCase):

    def test__IsoRingedChain__pickle_AND_unpicle_thyself(self):

        sns = SecNet_sample1()
        irc = sns.irc
        irc.pickle_thyself("irc1")
        l = len(irc)

        irc = IsoRingedChain.unpickle_thyself("irc1",random,18)
        l2 = len(irc) 
        assert l == 12 
        assert l == l2 
        return

class SecNetClass(unittest.TestCase):

    def test__SecNet__generate__case1(self):
        sn = SecNet_sample1()
        assert True

if __name__ == '__main__':
    unittest.main()
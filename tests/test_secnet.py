from secnet import *
import unittest 
import time

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

    def test__SecNet__generate__case2(self):
        ts = time.time()
        sn = SecNet_sample_TDirNv1()
        ts2 = time.time() - ts

        assert ts2 <= 30.0 # seconds 
        return

    def test__Secnet__unpickle_thyself(self):
        fp = "codename__ASS_SHIT"

        if os.path.isdir(fp):
            sn = SecNet.unpickle_thyself(fp,\
            DEFAULT_SINGLETON_RANGE,random,9)
            assert type(sn) == SecNet
        
        print("<SECNET> pickle function may not work!")

    def test__SecNet__subgraph_for_TDir(self):
        ##print("-- unpickling")
        sn = SecNet.unpickle_thyself("codename__ASS_SHIT",\
                DEFAULT_SINGLETON_RANGE,random,9)
        ##print("UNPICKLED")

        irc = sn.irc 

        td1 = TDir(13,17,\
                "I",radius=4,velocity=1)

        sg1 = sn.subgraph_for_TDir(td1) 

        q = set(sg1.sp[17].min_paths.keys())
        assert q == {5, 8, 9, 13, 17,\
                25, 28, 34} 

        assert sg1.d == defaultdict(set,\
                {17: {8, 5},\
                34: {25, 9}, 5: {17},\
                8: {17, 28}, 9: {34},\
                28: {8, 13}, 13: {25, 28},\
                25: {34, 13}})

        assert sg1.sn == {34, 8, 9, 13, 25, 28}

        qx = sg1.sp[13]

        qxr = qx.min_paths 
        for (k,v) in qxr.items():
                v_ = v[0]
                assert sum(v_.pweights) <= 4
        return 

if __name__ == '__main__':
    unittest.main()
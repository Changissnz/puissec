from cracker import *
import unittest

### lone file test 
"""
python3 -m tests.test_cracker
"""
###
class BackgroundInfoClass(unittest.TestCase):

    def test__BackgroundInfo__generate_Pr_case1(self):
        s = SecNet_sample1(SecSeq_sample_2(7,31)) 
        for s_ in s.irc.irl:
                s_.explode_contents()

        irc = s.irc
        random.seed(773)

        #lm = BackgroundInfo.generate_background_leak(irc,random)
        bi = BackgroundInfo.generate_instance(irc,s.srm)
        assert type(bi) == BackgroundInfo
        print(bi.dm)

    def test__BackgroundInfo__generate_instance_case1(self):

        fp = "codename__ASS_SHIT"

        if not os.path.isdir(fp):
            pickled_SecNet_sample_Q()

        sn = SecNet.unpickle_thyself(fp,\
                DEFAULT_SINGLETON_RANGE,random,9)

        irc = sn.irc 
        srm = sn.srm

        #random.seed(773)
        random.seed(777)

        #lm = BackgroundInfo.generate_background_leak(irc,random)

        bi = BackgroundInfo.generate_instance(irc,srm)
        assert type(bi) == BackgroundInfo
        print(bi.dm)


    

if __name__ == '__main__':
    unittest.main()
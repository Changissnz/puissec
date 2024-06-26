from cracker import *
import unittest,time 

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

    def test__BackgroundInfo__generate_instance_case2(self):
        sn = SecNet_sample_TDir1v1()

        irc = sn.irc 
        srm = sn.srm

        bi = BackgroundInfo.generate_instance(irc,srm)

        ks = set(bi.opm[0].keys())
        assert ks == {5, 9, 3, 2, 8}

        dv = {5:12,9:8,3:13,2:42,8:110}

        for k,v in bi.opm[0].items():
                #print("K: ",k, " V: ",len(v))
                assert len(v) == dv[k]
        assert bi.dec_map == {0: 2}


    def test__BackgroundInfo__naive_IRC2HypStruct_map(self):

        sn = SecNet.unpickle_thyself("codename__ASS_SHIT",\
                DEFAULT_SINGLETON_RANGE,random,9)
        irc = sn.irc 

        srm = sn.srm
        random.seed(777)

        ircm = BackgroundInfo.naive_IRC2HypStruct_map(irc,\
                full_hypseq=True,naive_split=2)

        ircm2 = BackgroundInfo.naive_IRC2HypStruct_map(irc,\
                full_hypseq=False,naive_split=2)

        assert set(ircm[0].keys()) == {5,9,3,2,8}
        assert set(ircm2[0].keys()) == set(ircm[0].keys())

        assert len(ircm[0][5]) == 12
        assert len(ircm2[0][5]) == 1
        return

class OrderOfCracknClass(unittest.TestCase):

    def test__OrderOfCrackn__order_by_depchain_map__case1(self):
        sn1 = SecNet.unpickle_thyself("codename__ASS_SHIT",\
                DEFAULT_SINGLETON_RANGE,random,9)

        irc = sn1.irc 
        srm = sn1.srm

        bi = BackgroundInfo.generate_instance(irc,srm)

        # testing 
        q = time.time()
        ooc = OrderOfCrackn()
        soln = ooc.order_by_depchain_map(bi.dm)
        qs = time.time() - q
        ##print("process time: ",qs)

        ans = [{1}, {0, 4}, {7}, {6}, {8}, {3}, {2}, {5}]
        assert soln == ans 
        return

    def test__OrderOfCrackn__order_by_depchain_map__case2(self):
        random.seed(243423)
        np.random.seed(324252)

        ss,sndg = SecSeq_sample_5(num_secs=3,max_nconn_ratio=0.0,\
        depconn_ratio=1.0,min_components=1,drange_max=3)

        m1 = ss.sec_instances_to_supermap('l')
        m2 = ss.sec_instances_to_supermap('d')
        m3 = ss.sec_instances_to_supermap('c')
        srm = SRefMap(m1,m2,m3,'c')

        irc = IsoRingedChain(ss,DEFAULT_SINGLETON_RANGE,\
        random,24153)
        bi = BackgroundInfo.generate_instance(irc,srm)

        ooc = OrderOfCrackn()
        soln = ooc.order_by_depchain_map(bi.dm)

        assert soln == [{2}, {0, 1}]
        return

    def test__OrderOfCrackn__order_by_depchain_map__case3(self):

        random.seed(243423)
        np.random.seed(324252)

        ss,sndg = SecSeq_sample_5(num_secs=3,\
        singleton_range=DEFAULT_SINGLETON_RANGE,\
        num_conn=5000,min_components=1,\
        max_nconn_ratio=0.0,drange_max=3,\
        depconn_ratio=1.0,conn_types=[1])

        m1 = ss.sec_instances_to_supermap('l')
        m2 = ss.sec_instances_to_supermap('d')
        m3 = ss.sec_instances_to_supermap('c')
        srm = SRefMap(m1,m2,m3,'c')

        irc = IsoRingedChain(ss,DEFAULT_SINGLETON_RANGE,\
        random,24153)
        bi = BackgroundInfo.generate_instance(irc,srm)

        ooc = OrderOfCrackn()
        soln = ooc.order_by_depchain_map(bi.dm)

        assert bi.dm == {0: [{0}], \
                1: [{1}, {0, 2}], 2: [{2}, {0}]}
        assert soln == [{0}, {2}, {1}]

if __name__ == '__main__':
    unittest.main()
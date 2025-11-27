from agents.cracker import *
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

    def test__BackgroundInfo__partially_naive_IsoRing2HypStruct_map_case1(self):

        sn = SecNet_sample_approxhyp()

        ## declare a <BackgroundInfo> 
        ir = sn.irc[0]

        b1,b2 = BackgroundInfo.partially_naive_IsoRing2HypStruct_map(\
                ir,1.0,0.3,0.7,random)

        ks = list(b2.values())

        c = Counter()
        for k in ks:
            c[k[0]] += 1

        cx = Counter({1:4,0:2})
        assert cx == c

        ############### CASE: entire <IsoRingedChain> map 

        ph1,ph2 = BackgroundInfo.partially_naive_IRC2HypStruct_map(sn.irc,\
                1.0, [0.,1.],[0.,1.],random)

        assert len(ph1) == 1
        assert len(ph2) == 1 
        assert len(ph1[0]) == len(ph2[0])
        #assert False,"PH1\n{}\PH2\n{}".format(ph1,ph2)
        assert len(ph1[0]) == 6 
        return

    def test__BackgroundInfo__naive_populate_IRC2HypStruct_map_case1(self):
        sn = SecNet_sample_approxhyp()

        irc = sn.irc 
        bound_length = 1.0 
        rd_range = [0.,0.5]
        ra_range = [0.3,0.7]
        rnd_struct = random 

        irc2h_map = BackgroundInfo.partially_naive_IRC2HypStruct_map(\
                irc,bound_length,rd_range,ra_range,rnd_struct)

        dx = deepcopy(irc2h_map[0]) 
        populate_ratio_range = [0.5,0.5]
        m2 = BackgroundInfo.naive_populate_IRC2HypStruct_map(irc,dx,populate_ratio_range,\
                rnd_struct=random,scramble=False)
        m2_ = m2[0]

        populate_ratio_range = [1.0,1.0]
        dx2 = deepcopy(irc2h_map[0]) 
        m3 = BackgroundInfo.naive_populate_IRC2HypStruct_map(irc,dx2,populate_ratio_range,\
                rnd_struct=random,scramble=False)
        m3_ = m3[0]

        ir = sn.irc[0] 
        for s in ir.sec_cache:
            d = s.dim()
            l = len(s.opm)

            ##print("{}: {}".format(l / 2,len(m2_[d])))

            # case 1 
            assert abs(len(m2_[d]) - l / 2) <= 1.0
            # case 2 
            assert l == len(m3_[d])

    # TODO: re-write
    """
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
    """

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

    # NOTE: tests checks for correctness of optional argument
    #       `malpermute_degree` in <BackgroundInfo.generate_instance>
    def test__BackgroundInfo__generate_instance_case3(self):

        s = SecNet_sample1(SecSeq_sample_2(7,31)) 
        for s_ in s.irc.irl:
            s_.explode_contents()

        irc = s.irc
        random.seed(773)

        bi = BackgroundInfo.generate_instance(irc,s.srm,True,0.0,0.3)
        assert type(bi) == BackgroundInfo
        q1 = bi.dm
        q1_ = bi.cdm

        random.seed(773) 
        bi = BackgroundInfo.generate_instance(irc,s.srm,True,0.5,1.0)
        assert type(bi) == BackgroundInfo
        q2 = bi.dm
        q2_ = bi.cdm

        assert q1 != q2
        assert q1_ != q2_ 

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

        ######################################################

class OrderOfCracknClass(unittest.TestCase):

    """
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
    """

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
 
class CrackerClass(unittest.TestCase):

    # NOTE: test partially checks for correctness; only asserts 
    #       output of method is a <Cracker> 
    def test__Cracker__generate_instance_by_engineered_BI__case1(self):

        sn = SecNet_sample_TDirNv1()

        irc = sn.irc
        srm = sn.srm
        rnd_struct = random

        dmapargs = ['cd',max,[0,1]]
        bi_args = [False,0.23,0.5,False,dmapargs] 
        irc2hs_args = ["naive",True,1,5] 

        crackling_sz = 6
        radar_radius = 5
        energy = 1000.0 

        cr = Cracker.generate_instance_by_engineered_BI(irc,srm,\
        rnd_struct,bi_args,irc2hs_args,crackling_sz,\
        radar_radius,energy)

        assert type(cr) == Cracker 

if __name__ == '__main__':
    unittest.main()
from psec_env import * 
import unittest

### lone file test 
"""
python3 -m tests.test_crackling
"""
###

class HypInferClass(unittest.TestCase):

    def test__HypInfer__infer_by_LeakInfo_case1(self):

        sn = SecNet_sample_approxhyp()

        irc = sn.irc
        srm = sn.srm
        bi = BackgroundInfo.generate_instance(irc,srm)

        ph1,ph2 = BackgroundInfo.partially_naive_IRC2HypStruct_map(sn.irc,\
                1.0, [0.,1.],[0.,1.],random)

        for x in sn.irc.irl:
            x.set_isorep(0)

        crck = Cracker(ph1,bi,6) 

        se = SecEnv(sn,crck,vb=0)
        se.preprocess()

        for _ in range(4):
            se.run(1.0)

        dim = sn.irc.irl[0].sec.dim()
        L = sn.irc.ircld.fetch_Leak(0,dim)

        li = L.leakm.d[0]

        assert len(crck.cracklings) == 0
        assert len(crck.spent) == 1 
        crackling = crck.spent[0] 

        hs = crackling.hs
        hsx = HypInfer.infer_by_LeakInfo(hs,li) 

        bsx = hsx.suspected_subbounds[0]
        stat = matrix_methods.point_in_bounds(bsx,irc.irl[0].sec.seq) 
        assert stat 
        return

if __name__ == '__main__':
    unittest.main()
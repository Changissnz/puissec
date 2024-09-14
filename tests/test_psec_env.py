from psec_env import *

import unittest,time 
import sys


### lone file test 
"""
python3 -m tests.test_psec_env 
"""
###

class SecEnvClass(unittest.TestCase):

    def test__SecEnv__init__case1(self):
        senv = SecEnv_sample_1()
        assert type(senv) == SecEnv
        return

    def test__SecEnv__instantiate_cracker_target__case1(self):

        senv = SecEnv_sample_1()
        senv.instantiate_cracker_target()

        cx = senv.crck.cracklings[0]
        td = cx.td

        q = td.td
        assert q.target_node == 0
        assert q.location == 1

        p = [1, 4]
        pw = [1]
        npath_actual = NodePath.preload(p,pw) 
        assert q.node_path == npath_actual
        assert senv.sn.node_loc_assignment[0] == q.node_path.p[-1]
        assert q.node_path.p[-1] == 4

    def test__SecEnv__instantiate_cracker_target__case2(self):
        se = SecEnv_sample_1(sn3=SecNet_sample_TDirNv1())
        se.instantiate_cracker_target()
        assert len(se.crck.cracklings) == 1

        # demonstrating the scores of the 
        # <Crackling>'s <TDirector> 
        c = se.crck.cracklings[0]
        assert c.loc() == 58
        c.td_next(1.0)
        assert c.loc() == 58

    ############################

    def test__SecEnv__InstantiateCrackerANDSecNet__case1(self):
        se = SecEnv_sample_2() 

        ##print("-- instantiate cracker target")
        se.instantiate_cracker_target()
        ##print("-- instantiate IRC-td")
        se.instantiate_td_for_IRC(5,1.0)

        radar_ans = {0:0,1:1}
        for q in se.sn.irc.irl:
            #print("SEC {}".format(q.sec.idn_tag))
            #print("RADIUS {}".format(q.td.td.radius))
            #print("$$$$$$$$$$$$$$")
            #print("* path")
            #print(str(q.td.td.node_path))    
            #print("SNGC")
            #q.td.resource_sg.display(0)
            lx = q.td.check_radar()
            ##print(lx)
            ##print("-----##")
            lx_ans = radar_ans[q.sec.idn_tag]
            assert len(lx) == lx_ans
            assert len(q.td.resource_sg.ring_locs) == 2
            assert len(q.td.resource_sg.crackling_locs) == 1 
            assert q.td.td.node_path.cost() == 0 

        ## for <Cracker>
        ##print("\t\t----------- CRACKER DATA------------")
        td_ = se.crck.cracklings[0].td
        assert td_.td.node_path.cost() == 2

    def test__SecEnv__cproc__case1(self):
        se = SecEnv_sample_1(sn3=SecNet_sample_TDirNv1())
        se.instantiate_td_for_IRC(5,1.0)

        se.instantiate_cracker_target()
        #assert len(se.crck.cracklings) == 1
        # demonstrating the scores of the 
        # <Crackling>'s <TDirector> 
        c = se.crck.cracklings[0]
        random.seed(1543)

        ir = se.sn.irc.fetch_IsoRing(18)
        ####
        """
        q = set(ir.td.resource_sg.d.keys())
        print("ISORING")
        for q_ in q: 
            s = ir.td.targetnode_analysis(q_,None)
            print("{}:{}".format(q_,s))
        """
        ####
        se.cproc()
        ####
        """
        print("[2] CRACKLING @ ",c.cidn)
        for q_ in q: 
            s = c.td.targetnode_analysis(q_,random)
            print("{}:\n{}\n".format(q_,s))
        """
        ####
        q67 = se.coloc()
        assert q67 == {'0,18,58'},"GOT {}".format(q67)

    def test__SecEnv__InstantiateCrackerANDSecNet__case2(self):

        se = SecEnv_sample_1(sn3=None)#SecNet_sample_TDirNv1())
        se.instantiate_td_for_IRC(5,1.0)
        se.instantiate_cracker_target()

        #assert len(se.crck.cracklings) == 1
        # demonstrating the scores of the 
        # <Crackling>'s <TDirector> 
        c = se.crck.cracklings[0]
        iring = se.sn.irc.irl[0]

        q = set(c.td.resource_sg.d.keys())
        # targetnode_analysis
        s = c.td.targetnode_analysis(np.min)
        return

    def test__SecEnv__load_IRCLD_into_SecNet__case1(self):

        se = SecEnv_sample_1(sn3=None)#SecNet_sample_TDirNv1())

        for x in se.sn.irc.irl:
            x.explode_contents() 

        se.load_IRCLD_into_SecNet()
        qx = se.sn.irc.ircld
        rx = set(qx.d[0].keys())

        assert rx == set(se.sn.irc.fetch_IsoRing(0).secdim_seq())
        assert set(qx.d.keys()) == {0}

    ### TODO: incorrect
    """
    def test__SecEnv__leak_by_str_idn__case1(self):

        se = SecEnv_sample_1(sn3=None)#SecNet_sample_TDirNv1())

        for x in se.sn.irc.irl:
            x.explode_contents() 

        se.load_IRCLD_into_SecNet()
        ##
        se.instantiate_td_for_IRC(5,1.0)
        se.instantiate_cracker_target()
        ##

        sidn = "0,0,32"
        lbsi = se.leak_by_str_idn(sidn)

        lbsi_ = se.crck.fetch_crackling(0).hs

        ans1 = np.array([\
            [0.13057814,0.17970186],\
            [0.,0.5],\
            [0.,0.5]])

        ans2 = np.array([\
            [0.,0.5],\
            [0.,0.5],\
            [0.,0.5]])

        assert matrix_methods.equal_iterables(\
            lbsi.suspected_subbounds[0],ans1), "got\n{}\nwant\n{}\n".format(\
                lbsi.suspected_subbounds[0],ans1) 

        assert matrix_methods.equal_iterables(\
            lbsi_.suspected_subbounds[0],ans2) 
    """
    
    # NOTE: print-test 
    def test__SecEnv__run__case1(self):

        orig_stdout = sys.stdout
        f = open('out_se.txt', 'w')
        sys.stdout = f

        se = SecEnv_sample_1(sn3=None)#SecNet_sample_TDirNv1())
        se.verbose = 2 

        for x in se.sn.irc.irl:
            x.explode_contents()

        se.load_IRCLD_into_SecNet()
        se.bi_update_to_cracker_hyp()
        ##
        se.instantiate_cracker_target()
        se.instantiate_td_for_IRC(5,1.0)
        
        for i in range(10):
            se.run(1.0)

        f.close()
        sys.stdout = orig_stdout

    def test__SecEnv__run__case2__SingleNodeInterdictionTest(self):

        sn = SecNet_sample_approxhyp()

        irc = sn.irc
        srm = sn.srm
        bi = BackgroundInfo.generate_instance(irc,srm)

        ###

        ph1,ph2 = BackgroundInfo.partially_naive_IRC2HypStruct_map(sn.irc,\
                1.0, [0.,1.],[0.,1.],random)

        for x in sn.irc.irl:
            x.set_isorep(0)

        crck = Cracker(ph1,bi,6) 

        se = SecEnv(sn,crck,vb=2)
        se.preprocess() 

        #se.load_IRCLD_into_SecNet()
        #se.instantiate_cracker_target()
        #se.instantiate_td_for_IRC(5,1.0)

        for _ in range(2):
            se.run(1.0)

        for _ in range(3):
            se.run(1.0)
            assert len(se.icrack) == 1

    def test__SecEnv__run__case3__demo(self):

        se = SecEnv_sample_1(sn3=None)#SecNet_sample_TDirNv1())
        se.verbose = 2 

        for x in se.sn.irc.irl:
            x.explode_contents()
        se.preprocess()
        for i in range(4):
            se.run(1.0)

    def test__SecEnv__run__case4__fstat(self):
        sn = SecNet_sample_TDirNv1()

        irc = sn.irc
        srm = sn.srm
        rnd_struct = random

        dmapargs = ['cd',max,[0,1]]
        bi_args = [False,0.0,None,False,dmapargs] 
        irc2hs_args = ["naive",True,1,5] 

        crackling_sz = 6
        radar_radius = 5
        energy = 1000.0 

        random.seed(224)
        crck = Cracker.generate_instance_by_engineered_BI(irc,srm,\
            rnd_struct,bi_args,irc2hs_args,crackling_sz,\
            radar_radius,energy)

        # open info mode
        se = SecEnv(sn,crck,rnd_struct=random,\
                ct_ratio=5000,vb=1,mode_open_info = (0,2))

        se.preprocess() 
        for _ in range(2):#30): 
            se.run(1.0)

        assert len(se.crck.cracklings) == 0
        assert len(se.cbs) == 0
        return



if __name__ == '__main__':
    unittest.main()
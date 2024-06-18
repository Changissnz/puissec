from psec_env import *
import unittest,time 

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
        """
        q = set(ir.td.resource_sg.d.keys())
        print("ISORING")
        for q_ in q: 
            s = ir.td.targetnode_analysis(q_,None)
            print("{}:{}".format(q_,s))
        """
        se.cproc()
        """
        print("[2] CRACKLING @ ",c.cidn)
        for q_ in q: 
            s = c.td.targetnode_analysis(q_,random)
            print("{}:\n{}\n".format(q_,s))
        """
        q67 = se.coloc()
        print("QQ: ",q67)
        assert q67 == {'0,18,58'}

if __name__ == '__main__':
    unittest.main()
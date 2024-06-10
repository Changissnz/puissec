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
        assert c.loc() == 11
        c.td_next(1.0)
        assert c.loc() == 10

if __name__ == '__main__':
    unittest.main()
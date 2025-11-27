from agents.secnet import *
import unittest 

### lone file test 
"""
python3 -m tests.test_tdir
"""
###

class TDirClass(unittest.TestCase):

    def test__TDir__load_path(self):
        sn = SecNet_sample_TDir1v1()
        q = sn.irc.irl[0].sec.idn_tag

        target_node = sn.node_loc_assignment[q]
        entry_points = sn.entry_points

        p1 = [16, 26, 21, 20, 5, 7, 19, 14, 23, 4, 10]
        w1 = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

        p2 = [25, 10]
        w2 = [1]

        p3 = [11, 28, 22, 27, 1, 13, 25, 10]
        w3 = [1, 1, 1, 1, 1, 1, 1]

        p4 = [23, 4, 10]
        w4 = [1, 1]

        npath1 = NodePath.preload(p1,w1) 
        npath2 = NodePath.preload(p2,w2)
        npath3 = NodePath.preload(p3,w3)
        npath4 = NodePath.preload(p4,w4) 

        tdsl = []
        d = {16:npath1,25:npath2,11:npath3,23:npath4}
        for ep in entry_points:
            qx = sn.sgc.sp[ep]
            xp = qx.min_paths[target_node]     
            td = TDir(ep,q,"C",30,1)
            tdsl.append(td)

            sg = sn.subgraph_for_TDir(td)
            td.load_path(sg)
            assert td.node_path == d[ep]

        td_ = tdsl[0]
        seq = [(16, True),(26, True),(26, True),(21, True),\
        (21, True),(20, True),(20, True),(5, True),\
        (5, True),(7, True)]
        j = None

        for i in range(30):
                q = td_.scaled__next__(0.5)
                if i < len(seq):
                        assert q == seq[i]
                if not q[1]:
                        j = i 
                        break 
        assert j == 20
        return

if __name__ == '__main__':
    unittest.main()
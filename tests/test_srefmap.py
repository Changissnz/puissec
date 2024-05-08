from secnet import * 
import unittest

### lone file test 
"""
python3 -m tests.test_srefmap
"""
###

class SRefMapClass(unittest.TestCase):

    def test__SRefMap__fc_proc__best_nodedec_map__case1(self):
        srm = SRefMap_sample1() 
        nm = srm.fc_proc__best_nodedec_map([0,1])
        assert len(nm) == 12

    def test__SRefMap__prmap_for_nodedec__case1(self):
        sn = SecNet_sample1()

        srm = sn.srm
        srm.reprocess('d')

        opm = srm.opmn[0]
        prmap = srm.prmap_for_nodedec(0,6,0,'greedy-d')

        assert opm[0] > prmap[0]
        assert prmap[6] > opm[6] 


if __name__ == '__main__':
    unittest.main() 
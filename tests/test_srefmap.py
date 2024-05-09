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
        nm = srm.fc_proc__best_nodedec_map(min,[0,1])
        assert len(nm) == 12

    def test__SRefMap__prmap_for_nodedec__case1(self):
        sn = SecNet_sample1()

        srm = sn.srm
        srm.reprocess('d')

        opm = srm.opmn[0]
        prmap = srm.prmap_for_nodedec(0,6,0,'greedy-d')

        assert opm[0] > prmap[0]
        assert prmap[6] > opm[6] 

        prmap = srm.prmap_for_nodedec(0,0,0,'greedy-d')
        assert opm[0] == prmap[0]
        assert prmap[6] == opm[6] 

    def test__SRefMap__collect_prism_points__DecMap__case1(self):
        ## use <SRefMap> on sample 3
        sn = SecNet_sample1() 
        srm = sn.srm

        dm = srm.collect_prism_points__DecMap('c',max,[0])  
        dm2 = srm.collect_prism_points__DecMap('d',max,[0])  
        dm3 = srm.collect_prism_points__DecMap('cd',max,[0])  

        dm4 = srm.collect_prism_points__DecMap('c',max,[1])
        dm5 = srm.collect_prism_points__DecMap('d',max,[1])
        dm6 = srm.collect_prism_points__DecMap('cd',max,[1])

        assert dm != dm2 
        assert dm == dm3 
        assert dm4 != dm5
        assert dm4 == dm6
        assert dm2 == dm5
        return
    
    # NOTE: dummy test; checks for crash-free exec.
    def test__SRefMap__collect_prism_points__PrMap__case1(self):
        sn = SecNet_sample1() 
        srm = sn.srm
        
        dx = srm.collect_prism_points__PrMap('c',"greedy-lone",0)
        dx2 = srm.collect_prism_points__PrMap('c',"greedy-d",0)
        dx3 = srm.collect_prism_points__PrMap('c',"greedy-c",0)
        assert True 

if __name__ == '__main__':
    unittest.main() 
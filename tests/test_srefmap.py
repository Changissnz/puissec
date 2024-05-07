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
        nm = srm.fc_proc__best_nodedec_map("greedy-lone",[0,1])
        assert len(nm) == 12

if __name__ == '__main__':
    unittest.main() 
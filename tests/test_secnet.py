from secnet import *
import unittest 

### lone file test 
"""
python3 -m tests.test_secnet
"""
###


class SecNetClass(unittest.TestCase):

    def test__SecNet__generate__case1(self):
        ss = SecSeq_sample_1(1)
        print(len(ss))

        sec_node_count = 12
        nsec_node_count = 23
        num_entry = 4
        rnd_struct = random
        sn = SecNet.generate(ss,sec_node_count,\
                nsec_node_count,num_entry,\
                rnd_struct,"spine frame",772) 

        assert True 

if __name__ == '__main__':
    unittest.main()
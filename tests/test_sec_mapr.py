from secnet_gen import * 
from iso_ring import * 
import unittest

### lone file test 
"""
python3 -m tests.test_sec_mapr
"""
###

class SecMaprFunctionsClass(unittest.TestCase):

    def test__SecMaprFunctions__metrics_on_node_in_depmap__case1(self):
        s = Sec_list_sample2()
        sndg = SecNetDepGen(s,random,1,0.8,[1,4])
        sndg.assign_conn(1500)
        #sq = sndg.sq
        s = sndg.sq 
        """
        for i in range(len(s)):
                for j in range(len(s)):
                        print("FOR S={}.{}",i,j)
                        metrcs = metrics_on_node_in_depmap(s[2].dm,j)
                        print(metrcs)
        """
        outp = metrics_on_node_in_depmap(s[1].dm,11)
        assert outp == (1, {4})


if __name__ == '__main__':
    unittest.main()
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

    def test__SecMaprFunctions__depchain_for_Sec__case1(self):

        ss = SecSeq_sample_1()
        sm = ss.sec_instances_to_supermap('d')

        ans = [{6}, {11, 5}]
        dc1 = depchain_for_Sec(sm,6)
        assert dc1 == ans 

        ans2 = [{10}, {0, 2, 8, 9, 11}]
        dc2 = depchain_for_Sec(sm,10)
        assert dc2 == ans2
        return

    def test__SecMaprFunctions__connected_subsets_of_codepmap__case1(self):
        ss = SecSeq_sample_1(1)
        sm = ss.sec_instances_to_supermap('c')
        cs = connected_subsets_of_codepmap(sm)
        assert len(cs) == 1

        ss2 = SecSeq_sample_1(5)
        sm = ss2.sec_instances_to_supermap('c')
        cs = connected_subsets_of_codepmap(sm)
        assert len(cs) == 5

if __name__ == '__main__':
    unittest.main()
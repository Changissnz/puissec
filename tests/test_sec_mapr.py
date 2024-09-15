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
        print("XXX")
        for i in range(len(s)):
                for j in range(len(s)):
                        print("FOR S={}.{}",i,j)
                        metrcs = metrics_on_node_in_depmap(s[2].dm,j)
                        print(metrcs)
        """

        outp = metrics_on_node_in_depmap(s[2].dm,7)
        assert outp == (1, {7})

    def test__SecMaprFunctions__depchain_for_Sec__case1(self):
        ss = SecSeq_sample_1()
        sm = ss.sec_instances_to_supermap('d')

        ans1 = [{1}, {9, 3}]
        dc1 = depchain_for_Sec(sm,1)
        assert dc1 == ans1
        return

    def test__SecMaprFunctions__connected_subsets_of_codepmap__case1(self):
        ss = SecSeq_sample_1(1)
        sm = ss.sec_instances_to_supermap('c')
        cs = connected_subsets_of_codepmap(sm)
        assert len(cs) == 2

        ss2 = SecSeq_sample_1(5)
        sm = ss2.sec_instances_to_supermap('c')
        cs = connected_subsets_of_codepmap(sm)
        assert len(cs) == 5

    def test__SecMaprFunctions__connected_subsets_of_codepmap__case2(self):
        s = Sec_list_sample2()
        sndg = SecNetDepGen(s,random,1,0.8,[1,4])
        sndg.assign_conn(5000,[2,3])
        ss = SecSeq(sndg.sq)

        sm = ss.sec_instances_to_supermap('c')
        cs = connected_subsets_of_codepmap(sm)
        assert len(cs) == 1

    def test__SecMaprFunctions__permute_setseq__case1(self):
        random.seed(51214)
        sseq = [{0,1},{3,4,5},{6,8,7,10}]
        rnd_struct = random
        sseq2  = permute_setseq(sseq,rnd_struct,20)
        ans = ([{0, 5}, {1, 4, 7}, {8, 10, 3, 6}], 20)
        assert sseq2 == ans



class PDMapIterClass(unittest.TestCase):

    def test__PDMapIter__next__case1(self):
        pdec = defaultdict(set,{8:{2,3,4},1:{0,1},2:{1,2,3},0:{4,3}})

        ans1 = defaultdict(None, {0: 3, 1: 0, 2: 1, 8: 2})
        ans6 = defaultdict(None, {0: 3, 1: 0, 2: 2, 8: 4})
        ans36 = defaultdict(None, {0: 4, 1: 1, 2: 3, 8: 4})
        ansd = {1:ans1,6:ans6,36:ans36}

        pmi = PDMapIter(pdec)
        while not pmi.reached_end():
                q = next(pmi)
                if pmi.sz in ansd: 
                        ans = ansd[pmi.sz]
                        assert ans == q 
        assert pmi.sz == 36 

if __name__ == '__main__':
    unittest.main()
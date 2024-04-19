from struct_samples import * 
import unittest

class SecClass(unittest.TestCase):

    def test__Sec__next__case1(self):
        singleton_range = [-1.,2.]
        sec = Sec_sample_1()
        tindex = -1

        for i in range(10):
            q = next(sec)

            if type(q[0]) == type(None):
                tindex = i
                break 

        assert tindex == 4 
        for k,v in sec.obfsr.bpoints.items():
            assert v.shape[0] == 4
            assert len(np.unique(v,axis=0)) == 4

            assert np.min(v) >= singleton_range[0]
            assert np.max(v) <= singleton_range[1]

        counter_ans = [Counter({0: 2, 1: 1, 2: 1}),\
                Counter({0: 3, 3: 1}),\
                Counter({2: 4}),\
                Counter({0: 4})]
        assert sec.obfsr.dpm.cnt == counter_ans

        vx = sec.lone_pr_vec_for_bloom()
        assert vx[0] == [0.2, 0.2, 0.4, 0.2]
        assert vx[1] == defaultdict(float,\
            {'0,1': 0.0})            
        return

    def test__Sec__next__case2(self):
        sec = Sec_sample_2() 
        sec.process_one_bloomiso()
        sec2 = sec.generate_next_Sec()
        assert len(sec2[0].seq) == 2
        assert sec2[0].obfsr.sz == 8
        assert not sec2[0].obfsr.tstat

        sec2[0].process_one_bloomiso() 
        sec3 = sec2[0].generate_next_Sec()
        assert sec3[0].obfsr.obf.oseeds.shape == (1296, 8)

if __name__ == '__main__':
    unittest.main()

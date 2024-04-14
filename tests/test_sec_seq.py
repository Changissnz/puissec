from sec_seq import * 
import unittest

def Sec_sample_1():
    sequence = np.array([0.4,1.8,2.0])
    singleton_range = [-1.,2.]

    o1 = matrix_methods.vector_to_string([0.2,1.5,0.8],float)
    o2 = matrix_methods.vector_to_string([-0.5,1.1,0.4],float)
    o3 = matrix_methods.vector_to_string([1.,1.,0.5],float)
    o4 = matrix_methods.vector_to_string(sequence,float)
    optima_pr_map = defaultdict(float)
    optima_pr_map[o1] = 0.2
    optima_pr_map[o2] = 0.2
    optima_pr_map[o3] = 0.2
    optima_pr_map[o4] = 0.4

    dep_map = defaultdict(float)
    codep_map = defaultdict(float)

    sec = Sec(sequence,singleton_range,optima_pr_map,\
        dep_map,codep_map)
    return sec

def Sec_sample_2():
    random.seed(3) 
    sequence = np.array([0.4,1.8,2])
    singleton_range = [-1.,2.]
    sec = Sec_sample_1()

    # add more points
    std_bound = np.array([deepcopy(singleton_range) \
        for _ in range(3)])
    ps = generate_pointgen_stdrand(std_bound,20,np.random)
    for p in ps:
        pstr = matrix_methods.vector_to_string(p,float) 
        sec.opm[pstr] = random.random()
    optima_pr_map = map_freq2pr(sec.opm)
    sec.opm = optima_pr_map

    # construct the dep. and co-dep. maps 
    bis = aprng_gauge.BatchIncrStruct(len(optima_pr_map),True,True,2)
    dep_map = defaultdict(float)
    codep_map = defaultdict(float)
    
    stat = True
    while stat: 
        q = next(bis) 
        stat = not type(q) == type(None) 
        if not stat: continue 
        vs = matrix_methods.vector_to_string(q,int)
        is_in_dep = 1 if random.random() >= 0.5 else 0 

        if is_in_dep:
            dep_map[vs] = random.random()
        else:
            codep_map[vs] = random.random() 

    sec = Sec(sequence,singleton_range,optima_pr_map,\
        dep_map,codep_map)

    return sec 

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

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
        return

if __name__ == '__main__':
    unittest.main()

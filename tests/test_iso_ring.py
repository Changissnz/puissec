from iso_ring import *
import unittest

# 
"""
python3 -m unittest tests.test_iso_ring
"""
#

def IsoRing_sample_1():
    random.seed(12)
    np.random.seed(12)

    singleton_range = [0.,1.] 
    dimension = 5
    num_optima = 2
    countermeasure = (0.6,0.5) 
    secs = []

    for i in range(5): 
            sec = Sec.generate_bare_instance(singleton_range,dimension,num_optima,\
            countermeasure,rnd_struct=np.random)
            secs.append(sec)

    sndg = SecNetDepGen(secs,random,2,0.75,[1,4])

    sndg.assign_conn()

    superbound = np.ones((5,2)) * np.array([0.,1.])
    spacing_ratio_range = [0.,0.2]
    outlier_pr_ratio = 0.4#1.0
    num_bounds = 8
    sb = deepcopy(superbound)

    obf = BoundedObjFunc.generate_BoundedObjFunc(\
            superbound,spacing_ratio_range,\
            outlier_pr_ratio,num_bounds,3) 

    return IsoRing(sndg.sq[0],obf,sb)

class BoundedObjFuncClass(unittest.TestCase):

    def test__BoundedObjFunc__output(self):
        # in global bounds (0,5) x 3 
        bounds_seq = [\
            [[0.,1],\
            [0.,1],\
            [0.,1]],\
            [[1,1.5],\
            [1,2.],\
            [1,5]],\
            [[1.5,5.],\
            [2,5],\
            [1,5]]]
        bounds_seq = [np.array(b) for b in bounds_seq]

        corr_objf = []
        corr_objf.extend([ObjFunc("e.d.",8),ObjFunc("s.m.",8),\
                ObjFunc("e.d.",10)])

        default_objf = ObjFunc("e.d.",12)
        do2 = deepcopy(default_objf)
        bof = BoundedObjFunc(bounds_seq,corr_objf,\
            default_objf)
        q = deepcopy(corr_objf) 


        # values 

        ## case 1 
        v1 = np.array([0.,0.5,0.75])
        ref1 = np.array([0.,0.25,0.])

        d1 = bof.output(ref1,v1) 
        ans1 = matrix_methods.euclidean_point_distance(ref1,v1)
        assert np.round(ans1,5) == d1 

        ## case 2 
        v2 = np.array([1.25,1.5,3])
        d2 = bof.output(ref1,v2) 
        ans2 = q[1].output(ref1,v2) 
        assert np.round(ans2,5) == d2 

        ## case 3
        v3 = np.array([1.25,3.,4.])
        d3 = bof.output(ref1,v3) 
        ans3 = do2.output(ref1,v3)
        assert np.round(ans3,5) == d3 

class IsoRingClass(unittest.TestCase):

    def test__IsoRing__explode_contents_case1(self):
        ir = IsoRing_sample_1()
        ir.explode_contents()
        assert len(ir.sec_cache) == 5

if __name__ == '__main__':
    unittest.main()
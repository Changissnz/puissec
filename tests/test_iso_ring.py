from iso_ring import *
import unittest

# 
"""
python3 -m unittest tests.test_iso_ring
"""
#

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

    def test__IsoRing__unpickle_thyself_case1(self):
        rx = IsoRing.unpickle_thyself('isosave')
        assert type(rx) == IsoRing
        ans_seq = np.array([0.083  , 0.67198, 0.80659, 0.98274, 0.63566])
        assert matrix_methods.equal_iterables(rx.sec.seq,ans_seq,5)

if __name__ == '__main__':
    unittest.main()
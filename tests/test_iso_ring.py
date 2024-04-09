from iso_ring import *
import unittest

# 
"""
python3 -m unittest tests.test_iso_ring
"""
#

class IsoRingClass(unittest.TestCase):

    def test__IsoRing_init(self):
        # typically a proper-bounds
        contents = np.array([[-1.0,0.0],\
                [0.0,1.0],\
                [1.0,20.0]])
        local_optima = None

        # args. for declaring `obj_func`
        arg1 = "s.m." # "e.d"  
        arg2 = 28 # random integer seed
        obj_func = ObjFunc(arg1,arg2)

        ir1 = IsoRing(contents,rnd_seed=34,\
            local_optima=None,obj_func=obj_func,\
                is_entry_point=False)

        # look at generated local_optima
        for x in ir1.local_optima:
            print(x)
            print()
            x_ = np.array(x) 
            assert morebs2.matrix_methods.point_in_bounds(contents,x_)
        return

class IsoChainedRingClass(unittest.TestCase):

    def test__IsoChainedRing_init(self):
        assert True 
        return


if __name__ == '__main__':
    unittest.main()
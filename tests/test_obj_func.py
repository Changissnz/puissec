from obj_func import *
import unittest

def point_pair_case_1():
    p1 = np.array([0.,0.,0.])
    p2 = np.array([5,2.5,5])
    return p1,p2

def template_obj_func_test(arg1,arg2,p1,p2):
    qx = ObjFunc(arg1,arg2)
    x1 = []
    for i in range(10):
        x1.append(qx.output(p1,p2))    
    return x1 

class ObjFuncClass(unittest.TestCase):

    def test__ObjFunc_output__case1(self):
        p1,p2 = point_pair_case_1()
        
        x1 = template_obj_func_test("e.d.",8,p1,p2)
        assert x1 == [7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5]

        x2 = template_obj_func_test("s.m.",8,p1,p2)
        assert x2 == [-22.5, 60.0, 52.5, -22.5, -22.5, -22.5, -22.5, -22.5, -22.5, 60.0]

        x21 = template_obj_func_test("s.m.",18,p1,p2)
        assert x21 == [45.0, -22.5, 60.0, -7.5, -75.0, 7.5, -60.0, 22.5, -45.0, 37.5]
        assert x21 != x2

        x3 = template_obj_func_test("r.n.",8,p1,p2)
        assert x3 == [7.76124, 7.34399, 7.4746, 7.51323, 7.63004, 7.58102, 7.30296, 7.53593, 7.32129, 7.71564]

        x31 = template_obj_func_test("r.n.",18,p1,p2)
        assert x31 == [7.64221, 7.53183, 7.54785, 7.48878, 7.3995, 7.56931, 7.49549, 7.58475, 7.4761, 7.5896]
        assert x31 != x3

if __name__ == '__main__':
    unittest.main()
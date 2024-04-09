from imod import * 
import unittest

### lone file test 
"""
python3 -m tests.test_imod
"""
###

random.seed(7) 

def imod_test_sample_input_1():
    drange = np.array([[-1,1.0],\
                    [2.,6.]])
    return drange 

class DefaultLPSModClass(unittest.TestCase):

    def test__DefaultLPSMod__load_DefaultLPSMod_function(self):
        drange = imod_test_sample_input_1()
        lpsm = DefaultLPSMod(None,drange,4,lambda x: int(x * 2))

        fx = DefaultLPSMod.load_DefaultLPSMod_function(lpsm,\
            np.add)
        assert True
        return

class Index2DModClass(unittest.TestCase):

    def test__Index2DMod__next__case1(self):
        drange = imod_test_sample_input_1()
        lpsm = DefaultLPSMod(None,drange,4,lambda x: int(x * 2))
        fx = DefaultLPSMod.load_DefaultLPSMod_function(lpsm,\
            np.add)
            
        bis = aprng_gauge.BatchIncrStruct(5,True,True,2)
        dim2d = np.array([[0,3],[0,6]])

        i2dm = Index2DMod(dim2d,bis,f=fx)
        tindex = None 

        for i in range(10000):
            if i2dm.fin_stat:
                print('finished at {}'.format(i))
                tindex = i
                break
            q = next(i2dm)
            #print("next: ", q)
        assert tindex == 26 

if __name__ == '__main__':
    unittest.main()
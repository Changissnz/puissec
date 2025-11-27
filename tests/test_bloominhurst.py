from base.bloominhurst import * 
import unittest

### lone file test 
"""
python3 -m tests.test_bloominhurst
"""
###

def OptimaBloomFunc_sample_args_1():
    oseeds = np.array([[0.75,3,2.5],\
            [1.75,0.3,0.25],\
            [0.8,2.9,2.8],\
            [1.7,1.3,2.5],\
            [0.75,1.3,2.5]])

    drange = np.array([[0,1.5],\
            [2,4.],\
            [2,3.]])
    return oseeds,drange

class OptimaBloomFuncClass(unittest.TestCase):

    def test__OptimaBloomFunc__next__case1(self):
        oseeds,drange = OptimaBloomFunc_sample_args_1()

        selector_func = None 
        q = np.add
        q2 = np.subtract 

        bloom_func = AltBaseFunc([q,q2],random)
        bloom_func = AltBaseFunc.load_AltBaseFunc_function(bloom_func)

        d = 5
        split_sz = 6
        splitsz_delta = DEFAULT_OPTIMA_BLOOM_SZ_DELTA_FUNC

        random.seed(1000)
        obf = OptimaBloomFunc(oseeds,drange,selector_func,\
            bloom_func,d,split_sz,splitsz_delta)

        fin_index = None
        for i in range(10):
            if obf.finished_stat:
                #print("finished at {}".format(i))
                fin_index = i
                break 
            next(obf)

        assert fin_index == 3
        return


    def test__OptimaBloomFunc__next__case2(self):
        oseeds,drange = OptimaBloomFunc_sample_args_1()

        selector_func = None 
        q = np.add
        q2 = np.subtract 

        bloom_func = AltBaseFunc([q,q2],random)
        bloom_func = AltBaseFunc.load_AltBaseFunc_function(bloom_func)

        d = 3
        split_sz = 6
        splitsz_delta = DEFAULT_OPTIMA_BLOOM_SZ_DELTA_FUNC

        random.seed(1000)
        obf = OptimaBloomFunc(oseeds,drange,selector_func,\
            bloom_func,d,split_sz,splitsz_delta)

        fin_index = None
        for i in range(10):
            if obf.finished_stat:
                #print("finished at {}".format(i))
                fin_index = i
                break 
            next(obf)

        assert fin_index == 5
        return


if __name__ == '__main__':
    unittest.main()

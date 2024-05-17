from leakette import * 
import unittest

### lone file test 
"""
python3 -m tests.test_leakette
"""
###
class LeaketteClass(unittest.TestCase):

    def test__Leakette__functiontest__case1(self):

        ir1 = IsoRing_sample_1()
        random.seed(1127) 
        degree = [1.,0.,(0.0,1.0)]
        fx = subbound_for_decimal

        q = leakf__type_MV(ir1,random,\
                degree,fx)

        ans = np.array([[0.91875, 0.91875],\
            [0.90071, 0.90071],\
            [0.03342, 0.03342],\
            [0.95695, 0.95695],\
            [0.13721, 0.13721]])

        assert matrix_methods.equal_iterables(ans,q)

        degree = [0.6,0.1,(0.0,1.0)] 
        q = leakf__type_MV(ir1,random,\
                degree,fx)

        ans = np.array([[0.81875,1.],\
                [np.nan,np.nan],\
            [0.,0.13342],\
            [np.nan,np.nan],\
            [0.03721,0.23721]])

        assert (ans[0] == q[0]).all()
        #assert np.sum(ans[1]) is np.nan
        assert (ans[2] == q[2]).all()
        #assert np.sum(ans[3]) is np.nan
        assert sum(abs((ans[4] - q[4]))) <= 10 ** -6

    def test__Leakette__functiontest__case2(self):

        ir1 = IsoRing_sample_1()
        random.seed(1127) 
        degree = 0.6
        fx = idn_decimal

        q = leakf__type_MV(ir1,random,\
                degree,fx)
        ans = np.array([0.91875,np.nan,\
                0.03342,np.nan,\
                0.13721])
        for (a,q_) in zip(ans,q):
            if np.isnan(q_):
                assert np.isnan(a) 
            else: 
                assert abs(a - q_) <= 10 ** -6

    def test__Leakette__functiontest__case3(self):

        ir1 = IsoRing_sample_1()
        random.seed(1127) 
        degree = (0.,0.)
        fx = choose_multiple

        q = leakf__type_MV(ir1,random,\
                degree,fx)
        assert np.all(np.isnan(q))

        degree = (1.,1.)
        q = leakf__type_MV(ir1,random,\
                degree,fx)
        assert matrix_methods.equal_iterables(q,ir1.sec.seq)

        degree = (1.0,0.0)
        q = leakf__type_MV(ir1,random,\
                degree,fx)
        ans = np.array([1.e-05,1.e-05,\
                1.e-05,1.e-05,\
                1.e-05])
        assert matrix_methods.equal_iterables(ans,q)

if __name__ == '__main__':
    unittest.main()
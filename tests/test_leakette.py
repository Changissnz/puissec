from secnet import * 
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
        return

    def test__leakf__type_MV__case1(self):
        random.seed(2004)
        np.random.seed(2004)

        ss = SecSeq_sample_4(num_secs=1,\
                singleton_range=DEFAULT_SINGLETON_RANGE,\
                num_conn=1,min_components=1,max_nconn_ratio = 0.3,\
                drange_max=1)

        sc = ss[0] 
        bound = [0.,1.]
        irc = IsoRingedChain(sc,bound,random,71)

        q = irc[0]
        degree = (1.,0.1,(0.,1.),2)
        fx = subbound_for_decimal_with_hop
        qx = leakf__type_MV(q,random,degree,fx)

        qx2 = np.array([[ 0.48775,  0.68775],\
        [ 0.70444,  0.90444],\
        [ 0.37034,  0.57034],\
        [ 0.3857 ,  0.5857 ],\
        [ 0.65874,  0.85874],\
        [ 0.69192,  0.89192],\
        [ 0.48632,  0.68632],\
        [ 0.90797,  1.03599],\
        [ 0.4632 ,  0.6632 ],\
        [-0.01441,  0.15677]])

        stat = matrix_methods.equal_iterables(qx,qx2)
        assert stat 


if __name__ == '__main__':
    unittest.main()
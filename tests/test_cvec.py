from base.cvec import * 
import unittest

### lone file test 
"""
python3 -m tests.test_cvec
"""
###

class CVecISelectorClass(unittest.TestCase):

    def test__CVecISelector__output__case1(self):
        ref_range = [0.9,1.0]
        prev_range = [0.0,0.2]

        cvis = CVecISelector(ref_range,prev_range)

        clvec = [80,20,10,5]

        ansd = {80:\
            ([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,\
            12, 13, 14, 15], [72, 73, 74, 75, 76, 77, 78, 79]),\
            20: ([0, 1, 2, 3], [18, 19]),\
            10: ([0, 1], [9]),\
            5: ([0], [4])}

        for c in clvec: 
            outp = cvis.output(c)
            assert ansd[c] == outp 
        return 

class CVecClass(unittest.TestCase):

    def test__CVec__cmp__case1(self):
        ciseq = default_cvec_iselector_seq()
        cvec = CVec(cvis=ciseq)

        x = [6,5,4,3,2,1,0.5,0.25,0.125]
        bvec = cmp_seq_with_cvec(cvec,x) 
        assert np.all(bvec == np.array([ True,  True,  True,  True,  True,  True,  True])) 

        x = [5,5]
        bvec = cmp_seq_with_cvec(cvec,x) 
        assert np.all(bvec == np.array([False, False, False, False, False, False, False]))

        x = [5,4]
        bvec = cmp_seq_with_cvec(cvec,x) 
        assert np.all(bvec == np.array([ True, False,  True, False, False, False, False]))

        x = [-1,-2,-1,-2,0,0,1,1,-1,-2,-2,-3]
        bvec = cmp_seq_with_cvec(cvec,x) 
        assert np.all(bvec == np.array([False,  True, False, False, False, False, False]))

    def test__CVec__scan_kmult_search(self):
        # case 1 
        V_m = np.array([33,36.,6,40,12,16])
        V_f = np.array([11,9.,2,8.,1,4])

        bs,me = CVec__scan__kmult_search(V_m,V_f,depth=5)
        assert round(abs(bs - 3.82),5) == 0.0

        # case 2 
        V_m = np.array([51.,36.,21.])
        V_f = np.array([3.,4.,7.]) 
        bs,me = CVec__scan__kmult_search(V_m,V_f,depth=5)
        assert round(abs(bs - 5.99999999),5) == 0.0




if __name__ == '__main__':
    unittest.main()
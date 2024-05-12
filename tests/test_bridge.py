from bridge import * 
import unittest

### lone file test 
"""
python3 -m tests.test_bridge
"""
###
class CBridgeClass(unittest.TestCase):

    """
    demonstration w/ wrong hypothesis
    """
    def test__CBridge__next__case1(self):

        ir = IsoRing_sample_1()
        ir.explode_contents()

        # declare a hypothesis 
        sb1 = np.ones((5,2)) * np.array([0.,0.5]) 
        suspected_subbounds = [sb1]
        hs = HypStruct(0,1,suspected_subbounds,sb_pr=np.array([1.0]))

        # declare a Crackling
        c = Crackling()
        c.load_HypStruct(hs)

        cb = CBridge(c,ir,hs,ssih=5)

        stat = True 
        i = 0 
        while stat:
                qc = next(cb)
                ##print("q: ",qc)
                stat = type(qc) != type(None)
                stat = stat and (i < 100)
                i += 1 

    """
    demonstration w/ the perfect hypothesis
    """
    def test__CBridge__next__case2(self):

        ir = IsoRing_sample_1()
        ir.explode_contents()

        hs = one_correct_HypStruct_for_IsoRing(ir)

        c = Crackling()
        c.load_HypStruct(hs)

        cb = CBridge(c,ir,hs,ssih=5)
        qc = next(cb)

        assert matrix_methods.equal_iterables(qc,ir.sec.seq,5)
        return

    def test__CBridge__next__case3(self):
        return -1 

if __name__ == '__main__':
    unittest.main()
  
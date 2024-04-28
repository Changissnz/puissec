from bridge import * 
import unittest

class CBridgeClass(unittest.TestCase):

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

if __name__ == '__main__':
    unittest.main()
  
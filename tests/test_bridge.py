from secnet import * 
import unittest

def CBridge_case_X(guess_index=3):


    bound = np.zeros((5,2),dtype=float)
    bound[:,1] = 2.0
    spacing_ratio_range = [0.5,0.6]
    outlier_pr_ratio = 0.3
    num_bounds = 10
    rs_seed = 1225 
    
    bof = BoundedObjFunc.generate_BoundedObjFunc(\
        deepcopy(bound),spacing_ratio_range,\
        outlier_pr_ratio,num_bounds,rs_seed)

        # make the <Sec>
    np.random.seed(1229)
    singleton_range = [0.,2.]
    dimension = 5
    num_optima = 4
    optima_countermeasure = (0.6,0.8)
    s1 = Sec.generate_bare_instance(singleton_range,\
        dimension,num_optima,optima_countermeasure,\
        rnd_struct=np.random)
    s1.idn_tag = 2
    
    
    bs = np.array([s1.seq - 0.2,s1.seq + 0.2]).T
    
    ir = IsoRing(s1,bof,bound,None)
    
    # make a HypStruct and a Crackling
    suspected_subbounds = [deepcopy(bs)]
    hs = HypStruct(s1.idn_tag,guess_index,\
        suspected_subbounds,sb_pr=np.array([1.0]))
        
    c = Crackling(cidn=1,cvsz=100)
    c.load_HypStruct(hs)
    
    cb = CBridge(c,ir,hs,ssih=2,cidn=1,\
        batch_size=2000,verbose=True)
    return cb 

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
    ### TODO: check!
    """
    def test__CBridge__next__case2(self):

        ir = IsoRing_sample_1()
        ir.explode_contents()

        hs = one_correct_HypStruct_for_IsoRing(ir)

        c = Crackling()
        c.load_HypStruct(hs)

        cb = CBridge(c,ir,hs,ssih=5)
        qc = next(cb)
        print("assert case2")
        #assert matrix_methods.equal_iterables(qc,ir.sec.seq,5),\
        #        "got {}\nwant {}".format(qc,ir.sec.seq)
        cv = c.cvec
        assert len(c.cvec.v) == 1
        assert len(c.cvec.input_samples) == 1
        assert c.cvec.v[0] == 0.0
        return
    """

    ## TODO: REMOVE!
    """
    def test__CBridge__next__case3(self):
        ir = IsoRing_sample_1()
        ir.explode_contents()
        sec = ir.sec.deepcopy(new_idn_tag=ir.sec.idn_tag,transfer_obfsr=False)

        seq_idn = sec.idn_tag
        ti = ir.sec.seq_index()
        bs = np.array([0.91875, 0.90071, 0.03342, 0.95695, 0.13721])
        q = bs - 0.1
        q2 = bs + 0.1
        sbs = [np.array([q,q2]).T]
        hs = HypStruct(seq_idn,ti,sbs,sb_pr=np.array([1.0]))

        c = Crackling()
        c.load_HypStruct(hs)

        cb = CBridge(c,ir,hs,ssih=2)

        for i in range(32): 
            qc = next(cb) 

        ans = defaultdict(float, {'0.91875,0.90071,0.03342,0.95695,0.13721': 0.4})
        assert c.cracked_dict == ans, "want {}\ngot {}".format(c.cracked_dict,\
                ans)
    """

    ## NOTE: remove 
    """
    def test__CBridge__next__case5(self):

        # generate the <IsoRing>
        bound = np.zeros((6,2),dtype=float)
        bound[:,1] = 2.0
        spacing_ratio_range = [0.,0.05]
        outlier_pr_ratio = 1.0
        num_bounds = 2
        rs_seed = 1224 

        bof = BoundedObjFunc.generate_BoundedObjFunc(\
        deepcopy(bound),spacing_ratio_range,\
        outlier_pr_ratio,num_bounds,rs_seed)

        np.random.seed(1224)
        singleton_range = [0.,2.]
        dimension = 6
        num_optima = 3
        optima_countermeasure = (0.3,0.2)
        s1 = Sec.generate_bare_instance(singleton_range,\
        dimension,num_optima,optima_countermeasure,\
        rnd_struct=np.random)
        s1.idn_tag = 2

        ir = IsoRing(s1,bof,bound,None)

        # make a HypStruct and a Crackling
        suspected_subbounds = [deepcopy(bound)]
        hs = HypStruct(s1.idn_tag,1,suspected_subbounds,sb_pr=np.array([1.0]))
        print("HYPE")
        print(str(hs))
        print()

        c = Crackling(cidn=1,cvsz=100)
        c.load_HypStruct(hs)

        cb = CBridge(c,ir,hs,ssih=3,cidn=1,\
        batch_size=2000,verbose=False)#True)
        i = 0 

        while i < 2000:
                ##print("I: ",i)
                cbx = next(cb)
                if type(cbx) == type(None):
                        print("NO MORE")
                        break 
                ##print(cbx)
                i += 1 

        assert i == 732, "got {} want 732".format(i)
        assert len(c.flagged_pts) == 260 
    """

    ############################### 
    
    ## NOTE: remove 
    """
    def test__CBridge__next__case4(self):
        ss = SecSeq_sample_2(40,200)
        bound = [0.,1.]
        irc = IsoRingedChain(ss,bound,random,71)

        ir = irc[33]
        ir.explode_contents()
        bds = np.ones((5,2)) * np.array([0.,1.])
        ir.ofunc = BoundedObjFunc.one_simple_BoundedObjFunc(\
                bds,12)

        B = np.array([[0.60,0.75],\
                [0.75,0.9],\
                [0.5,0.8],\
                [0.72,0.82],\
                [0.55,0.75]])

        ti = irc[33].sec.seq_index()
        hs = HypStruct(33,ti,\
                [B],sb_pr=np.array([1.]))

        ssih = 9
        c = Crackling()
        c.load_HypStruct(hs)
        cb = CBridge(c,ir,hs,ssih=ssih)

        for i in range(200): 
                q = next(cb) 
                if type(q) == type(None):
                        #print("early T=",i)
                        break 
                #print("Q: ",q)
        assert len(c.flagged_pts) == 46, "got {} want 46".format(\
            len(c.flagged_pts))
    """

    ## NOTE: remove 
    """
    def test__CBridge__next__case6(self):
        cb = CBridge_case_X(guess_index=3)
        
        i = 0 
        while i < 35:
            print("I: ",i)
            cbx = next(cb)
            if type(cbx) == type(None):
                print("NO MORE")
                break 
            i += 1 

        s = cb.isoring.sec
        prm = s.optima_points_to_index_pr_map()
        crackling = cb.crackling 
        assert len(crackling.cracked_dict) == 1, "want {} got {}".format(1,\
            len(crackling.cracked_dict)) 
        assert set(crackling.cracked_dict.values()) == \
                {prm[3]}, "got values\n{}".format(set(crackling.cracked_dict.values()))
        return 
    """

if __name__ == '__main__':
    unittest.main()
  
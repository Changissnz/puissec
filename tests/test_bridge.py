from secnet import * 
import unittest

# TODO: unused,delete?
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
        suspected_subbounds,sb_pr=np.array([1.0]),\
        hs_vec=np.array([2]))
        
    c = Crackling(cidn=1,cvsz=100)
    c.load_HypStruct(hs)
    
    cb = CBridge(c,ir,cidn=1,\
        batch_size=2000,verbose=True)
    return cb 
###############################################

def CBridge_case_ApproximateHypStruct_num0(is_approx_hypstruct:bool,\
    hop_size=5):

    random.seed(2004)
    np.random.seed(2004)

    ss = SecSeq_sample_4(num_secs=1,\
            singleton_range=DEFAULT_SINGLETON_RANGE,\
            num_conn=1,min_components=1,max_nconn_ratio = 0.3,\
            drange_max=1)

    sc = ss[0] 
    print(sc[0]) 

    bound = [0.,1.]
    irc = IsoRingedChain(sc,bound,random,71)
    for x in irc.irl:
        x.explode_contents()
    q = irc[0]

    # [10, 6, 3, 2, 8, 9]
    seci = 2
    bound_length = np.array([0.4,0.8,0.5])
    location_ratio = np.array([0.1,0.2,0.6])

    if is_approx_hypstruct:
        hs = one_approximate_HypStruct_for_IsoRing(q,seci,\
            bound_length,location_ratio,hop_size)
    else: 
        q.set_isorep(seci)
        hs = one_correct_HypStruct_for_IsoRing(q)

    hs_seq = [hop_size for _ in hs.hs_vec]
    hs.hs_vec = hs_seq 
    # run a <CRridge> cracking session
    cr = Crackling(cidn=2,cvsz=200)
    cr.load_HypStruct(hs) 

    cb  = CBridge(cr,q,cidn=12,batch_size=100,verbose=True)
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
        hs = HypStruct(0,1,suspected_subbounds,sb_pr=np.array([1.0]),\
            hs_vec=np.array([5]))

        # declare a Crackling
        c = Crackling()
        c.load_HypStruct(hs)

        cb = CBridge(c,ir)

        stat = True 
        i = 0 
        while stat:
                qc = next(cb)
                ##print("q: ",qc)
                stat = type(qc) != type(None)
                stat = stat and (i < 100)
                i += 1
        assert not c.astat

    def test__CBridge__next__ApproximateHypStruct__num0__case1(self):

        cb = CBridge_case_ApproximateHypStruct_num0(True) 
        cr = cb.crackling 
        q = cb.isoring 
        i = 0 
        il = 300
        while i < il and not cr.astat: 
            print("i: ",i)
            next(cb)
            i += 1

        assert i == 8

        l = len(cr.flagged_pts)

        print("STAT")
        print("- c : ",cr.cstat)
        print("- i : ",cr.istat)
        print("- a: ", cr.astat)

        print("FLAGGED POINTS")
        print(l)

        assert not cr.cstat and not cr.istat
        assert cr.astat 
        assert l == 4
        print(cr.flagged_pts)

        assert len(cr.cracked_dict) == 1 
        sx = '0.00515,0.62506,0.15195'
        vx = cr.cracked_dict[sx]
        assert vx == q.sec.opm[sx]

        qr1 = matrix_methods.string_to_vector(sx,float)
        qr1 = np.round(qr1,5) 
        assert (q.sec.seq == qr1).all() 

    def test__CBridge__next__ApproximateHypStruct__num0__case2(self):
        cb = CBridge_case_ApproximateHypStruct_num0(True,hop_size=10)
        cr = cb.crackling

        i = 0 
        il = 3000
        while i < il and not cr.astat: 
            print("i: ",i)
            next(cb)
            i += 1

        assert i == 127
        assert len(cr.flagged_pts) == 91
        assert len(cr.cracked_dict) == 1

    def test__CBridge__next__ApproximateHypStruct__num0__case3(self):
        cb = CBridge_case_ApproximateHypStruct_num0(False,hop_size=7)
        cr = cb.crackling

        i = 0 
        il = 3000
        while i < il and not cr.astat: 
            print("i: ",i)
            next(cb)
            i += 1

        assert i == 1
        assert len(cr.cracked_dict) == 1

if __name__ == '__main__':
    unittest.main()
  
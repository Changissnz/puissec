from leakette import * 
from morebs2 import relevance_functions
from hype import * 
from tdir import * 

class Crackling:

    def __init__(self,cmp_deg=1.5,cidn=None):
        self.cmp_deg = cmp_deg 
        self.cidn = cidn
        self.hs = None 
        self.td = None
        self.rss = None
        self.set_cvec()
        self.astat = False
        self.cracked_dict = defaultdict(float)
        self.flagged_pts = []

    def load_HypStruct(self,hs):
        assert type(hs) == HypStruct
        self.hs = hs 

    def target_of_HypStruct(self):
        if type(self.hs) == type(None): 
            return None
        return self.hs.seq_idn

    def load_TDirector(self,td): 
        assert type(td) == TDirector
        self.td = td
        return

    def fetch_tdir(self):
        if type(self.td) == type(None):
            return None
        return self.td.td 

    def loc(self):
        if type(self.td) == type(None):
            return None
        return self.td.loc()  

    """
    """
    def td_next(self,timespan=1.0):

        if self.td.obj_stat == "search for target":
            print("---")
            print("SEARCHING FOR TARGET")
            # check if <TDir> is still active.
            if not self.td.td.active_stat:
                self.td.extloc_search__set_TDir(extf=max,rnd_struct=random)
            print("path is ")
            print(str(self.td.td.node_path))
            print("------")
            self.td.td.scaled__next__(timespan)
            return
        elif self.td.obj_stat == "capture target":
            if not self.td.td.active_stat:
                isoring_loc = self.td.check_radar()
            assert False, "not programmed yet."
        return

    # TODO: test 
    """
    used by <Crackling> instance to decide 
    to accept an entry point as a point
    """
    def accept_TDirector_at_entry_point(self):
        assert type(self.td) == TDirector
        q = self.td.check_radar()
        return len(q) > 0

    def set_cvec(self):
        ciseq = default_cvec_iselector_seq()
        cvec = CVec(cvis=ciseq)
        self.cvec = cvec 

    """
    p := vector, the attempted point
    q := vector, length equal to number of optima,
            each i'th value is the score of p w.r.t.
            the i'th optimum point. 
    astat := bool, ?is `p` equal to any local optima 
                   of the targeted <IsoRing> instance.
    """
    def register_response(self,p,q,astat:bool): 
        assert type(self.hs) != type(None)
        self.astat = astat 
        ##print("PP: {}\nQQ: {}".format(p,q))
        x = q[self.hs.target_index] 
        self.cvec.append(x,p)
        
        qx1 = self.cvec.cmp(measures.zero_div)
        qx2 = self.cvec.cmp(np.less_equal) 
        ##print("QX1: ",qx1)
        ##print("QX2: ",qx2)

        ##s = [qx1[self.hs.target_index],qx2[self.hs.target_index]] 
        
        ## arg #1
        s = [np.sum(qx1) >= len(qx1) / self.cmp_deg,\
            np.sum(qx2) >= len(qx2) / self.cmp_deg]
        
        d = s[0] or s[1] 
        """
        d = metric_2dboolmat_to_bool(np.array([qx1,qx2]),\
            0.4,True)
        """
        if d: 
            self.flagged_pts.append(len(self.cvec) - 1)
        return d

    def register_lo_pr(self,prx):
        assert len(self.cvec.input_samples) > 0
        v = self.cvec.input_samples[-1]
        s = matrix_methods.vector_to_string(v,float)
        self.cracked_dict[s] = prx 
        return

"""
- return: 
a <RChainHead> instance derived from a 
<IsoRing>. 
"""
def IsoRing_and_Crackling_to_base_RChainHead(ir:IsoRing,cracklng:Crackling):
    assert type(ir) == IsoRing

    def outputf1(p):
        ##print("IR registers attempt")
        ##print(p)
        q,stat = ir.register_attempt(p)
        ##
        """
        print("-- register")
        print("scorevec: ",q)
        print("stat: ",stat)
        """
        ##
        """
        if stat:
            return None
        """
        ##print("cracking resp: ")
        d = cracklng.register_response(p,q,stat)
        if stat:
            prx = ir.response_to_prompt(cracklng.hs.target_index)
            ##print("PRXXX: ",prx)
            cracklng.register_lo_pr(prx) 
        ##print("response: ",d)
        return d

    rch = relevance_functions.RChainHead()

    argseq1 = ['nr',outputf1]
    rch.add_node_at(argseq1)
    return rch

"""
"""
def default_base_RSSI(ir:IsoRing,cracklng:Crackling,hs:HypStruct,ssih):
    assert type(ir) == IsoRing
    assert type(hs) == HypStruct

    ##print("\t-- converting to RCH")
    rch = IsoRing_and_Crackling_to_base_RChainHead(ir,cracklng)
    ##print("\t-- CONVERTED")
    resplattingMode = ("relevance zoom",rch)

    mpsb = hs.most_probable_subbound()
    start_point = deepcopy(mpsb[:,0])
    rss = rssi.ResplattingSearchSpaceIterator(mpsb,\
        start_point,SSIHop=ssih,resplattingMode = \
            resplattingMode)
    ##print("DECLARED RSS")
    return rss

"""
class w/ no variables; contains basic functions to 
improve a <HypStruct> instance. 
"""
class HypInfer: 

    def __init__(self):
        return

    @staticmethod
    def infer_by_LeakInfo(hypStruct,leak_info):
        assert type(hypStruct) == HypStruct
        assert type(leak_info) == LeakInfo

        for lk1,lk2 in leak_info.leak_info.items():
            for lk2_ in lk2:
                hypStruct = HypInfer.infer_FX(\
                    hypStruct,lk2_,lk1)
        return hypStruct

    # TODO: test 
    @staticmethod
    def infer_FX(hypStruct,leak_value,leak_idn):
        assert leak_idn in {0,1,2}

        def exec_fn(sx):
            if leak_idn == 0:
                return adjust_bounds__F0(sx,leak_value)
            elif leak_idn == 1:
                return adjust_bounds__F1(sx,leak_value)
            else: 
                return adjust_bounds__F2(sx,leak_value)

        for (i,s) in enumerate(hypStruct.suspected_subbounds):
            s_ = exec_fn(s)
            hypStruct.suspected_subbounds[i] = s_ 
        return hypStruct

def closest_reference_to_bound_start(b,V_f):

    null_vec = np.zeros((len(V_f,)))
    start_vec = deepcopy(null_vec) 

    for (i,v) in enumerate(V_f):
        if not np.isnan(v): 
            qx = measures.zero_div(b[i,0],v,np.nan)

            if np.isnan(qx):
                null_vec[i] = b[i,0]
                continue

            xq = int(round(qx))
            start_vec[i] = xq * v
        else: 
            null_vec[i] = b[i,0]

    return start_vec + null_vec 

def adjust_bounds__F0(b,V_f):

    r = closest_reference_to_bound_start(b,V_f)
    bx = deepcopy(b[:,1])

    ##print("ADJUSTING ",bx) 
    ##print("UNIT ",V_f)

    qi = [vf if (not np.isnan(vf)) else 0.0 for vf in V_f] 
    V_f = np.array(qi) 
    ##print("UNIT2: ",V_f)

    q = CVec__scan__kmult_search(bx,V_f,depth=1)
    ##print("Q0: ",q)
    q = q[0]
    ##print("Q: ",q)

    if type(q) in {type(None),tuple}: 
        print("k-mult search failed.")
        return b 

    ##print("V_F: ",V_f)
    ##print("Q: ",q)
    nb = np.array([r,r+ V_f*q]).T
    print("UPDATED:\n\t",nb)
    return nb 

def adjust_bounds__F1(b,V_f):
    bx = deepcopy(b)

    for (i,bx_) in enumerate(bx):
        if not np.isnan(V_f[i]):
            rx = np.array([V_f[i],V_f[i] + 0.05]).T
            bx[i] = rx
    return bx

def adjust_bounds__F2(b,V_f):
    ##assert matrix_methods.is_proper_bounds_vector(V_f) 
    assert V_f.shape[1] == 2 and V_f.shape[0] > 0 
    b2 = deepcopy(b)

    for (i,bx) in enumerate(V_f):
        stat1 = np.isnan(bx[0]) 
        stat2 = np.isnan(bx[1])
        stat = stat1 or stat2

        if stat: continue
        b2[i] = deepcopy(bx)

    return b2 
from iso_ring import * 
from morebs2 import relevance_functions
from hype import * 

class Crackling:

    def __init__(self,target_indices=None):
        self.hs = None 
        self.rss = None
        self.set_cvec() 
        self.astat = False
        self.cracked_dict = defaultdict(float)

    def load_HypStruct(self,hs):
        assert type(hs) == HypStruct
        self.hs = hs 

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
        
        s = [qx1[self.hs.target_index],qx2[self.hs.target_index]] 
        
        return s[0] or s[1] 
        #####
        """
        print("QXX")
        print(qx1)
        print(qx2)
        print()
        m = np.array([qx1,qx2])
        return metric_2dboolmat_to_bool(m,vx=0.65)
        """
        ##### 

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
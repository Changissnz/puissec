from iso_ring import * 
from morebs2 import relevance_functions
from hype import * 

class Crackling:

    def __init__(self,target_indices=None):
        self.hs = None 
        self.rss = None
        self.set_cvec() 
        self.astat = False

    def load_HypStruct(self,hs):
        assert type(hs) == HypStruct
        self.hs = hs 

    def set_cvec(self):
        ciseq = default_cvec_iselector_seq()
        cvec = CVec(cvis=ciseq)
        self.cvec = cvec 

    """
    only want 
    """
    def register_response(self,p,q,astat:bool): 
        assert type(self.hs) != type(None)
        self.astat = astat 
        x = q[self.hs.target_index] 
        self.cvec.append(x,p)
        
        qx1 = self.cvec.cmp(measures.zero_div)
        qx2 = self.cvec.cmp(np.less_equal) 

        m = np.array([qx1,qx2])
        return metric_2dboolmat_to_bool(m,vx=0.65)

    # TODO: 
    def load_prtable(self):
        return -1 


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
        print(q)
        print(stat)
        """
        ##
        
        if stat:
            return None
        ##print("cracking resp: ")
        d = cracklng.register_response(p,q,stat)
        ##print("\t",d)
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
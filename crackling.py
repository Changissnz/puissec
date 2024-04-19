from iso_ring import * 
from morebs2 import relevance_functions
from hype import * 

"""
- return: 
a <RChainHead> instance derived from a 
<IsoRing>. 
"""
def IsoRing_to_base_RChainHead(ir:IsoRing):
    assert type(ir) == IsoRing

    def outputf1(p):
        q,stat = ir.register_attempt(p)

        if stat:
            return None
        m = ir.cvec_feedback()
        return m

    def outputf2(m):
        ## NOTE 
        # USE W/ CAUTION HERE 
        if type(m) == type(None):
            return False
        return metric_2dboolmat_to_bool(m)

    rch = relevance_functions.RChainHead()

    argseq1 = ['nr',outputf1]
    rch.add_node_at(argseq1)

    argseq2 = ['nr',outputf2]
    rch.add_node_at(argseq2)
    return rch

"""

"""
def default_base_RSSI(ir:IsoRing,hs:HypStruct,ssih):
    assert type(ir) == IsoRing
    assert type(hs) == HypStruct

    rch = IsoRing_to_base_RChainHead(ir)
    resplattingMode = ("relevance zoom",rch)

    mpsb = hs.most_probable_subbound()
    start_point = deepcopy(mpsb[:,0])
    rss = rssi.ResplattingSearchSpaceIterator(mpsb,\
        start_point,SSIHop=ssih,resplattingMode = \
            resplattingMode)

    return rss

class Crackling:

    def __init__(self): 
        self.rss = None

    def declare_new_rssi(self,ir,hs,ssih):
        self.rss = default_base_RSSI(ir,lhs,ssih)
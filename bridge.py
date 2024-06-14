from crackling import * 

class CBridge:

    def __init__(self,crackling,isoring,hs,ssih=5,\
        cidn=None):
        assert type(crackling) == Crackling
        assert type(isoring) == IsoRing 
        assert type(hs) == HypStruct

        self.crackling = crackling 
        self.isoring = isoring
        self.hs = hs 
        self.rssi = None 
        self.load_crackf(ssih)
        self.cidn = cidn
        return

    def load_crackf(self,h=5):
        self.rssi = default_base_RSSI(self.isoring,self.crackling,\
            self.hs,h)

    """
    return:
    - (<Crackling> idn,<IsoRing> idn)
    """
    def agent_idns(self):
        return (self.crackling.cidn,self.isoring.sec.idn_tag)

    def __next__(self):
        p = next(self.rssi)
        return p
        
def one_correct_HypStruct_for_IsoRing(ir:IsoRing):
    ti = ir.sec.seq_index()
    seq = ir.sec.optima_points()[ti]
    b = np.array([seq,seq + 10 ** -4]).T
    sbs = [b] 
    sb_pr = np.array([1.0])
    return HypStruct(ir.sec.idn_tag,ti,sbs,sb_pr)
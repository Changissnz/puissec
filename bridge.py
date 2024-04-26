from crackling import * 

class CBridge:

    def __init__(self,crackling,isoring,hs,ssih=5):
        assert type(crackling) == Crackling
        assert type(isoring) == IsoRing 
        assert type(hs) == HypStruct

        self.crackling = crackling 
        self.isoring = isoring
        self.hs = hs 
        self.rssi = None 

        self.load_crackf(ssih)
        return

    def load_crackf(self,h=5):
        self.rssi = default_base_RSSI(self.isoring,self.crackling,\
            self.hs,h)

    def __next__(self):
        p = next(self.rssi)
        return p 



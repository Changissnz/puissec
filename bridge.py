from crackling import * 

class CBridge:

    def __init__(self,crackling,isoring,hs,ssih=5):
        assert type(crackling) == Crackling
        assert type(isoring) == IsoRing 

        self.crackling = crackling 
        self.isoring = isoring
        self.hs = hs 
        self.crackling.declare_new_rssi(self.isoring,hs,ssih)
        return

    def __next__(self):
        p = next(self.crackling.rssi)
        return p 



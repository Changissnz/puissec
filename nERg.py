from cit import *

class NerG:

    """
    value := float, non-negative value representing cumulative
             energy level of associated agent.
    rf := str|None; reference file used by instance
    """
    def __init__(self,value:float,rf:str=None):
        assert value > 0.0
        self.v = value
        self.open_rf(rf)
        return

    def open_rf(self,rf):
        self.rf = rf
        if type(rf) == type(None):
             return None
        # TODO: 
        assert False

    def default_crackling_size(self):
        return 1

    def default_crackling_velocity(self):
        return 1.0 

    def __add__(self,v):
        assert type(v) in {np.int32,np.float64,\
            np.float32,int,float}

        q = deepcopy(self)
        q.v = q.v + v
        return q

    def __sub__(self,v):
        assert type(v) in {np.int32,np.float64,\
            np.float32,int,float}

        q = deepcopy(self)
        q.v = q.v - v
        return q

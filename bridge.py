from iso_ring import *

class CBridge:

    def __init__(self,isoring):
        self.isoring = isoring
        return

    def msg_isoring(self,p):
        diff,stat = self.isoring.register_attempt(p)
        return diff,stat 



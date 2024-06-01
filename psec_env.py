from secnet import * 

class SecEnv:

    def __init__(self,sn,cracker):
        assert type(sn) == SecNet
        assert type(cracker) == Cracker

        self.sn = sn
        self.cracker = cracker
        return

    def load_crackling(self):

        return -1 
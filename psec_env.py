from secnet import * 

class TDBridge:

    def __init__(self):
        return 

    def pass_info__G(self,c:Crackling,g:SecNet):
        #tdr:TDirector,\):
        assert type(c) == Crackling
        assert type(g) == SecNet

        tdr = c.td.tdir

        G = g.subgraph_for_TDir(tdr)

        c.td.load_graph(G)
        return

"""
environment for the dual-agent 
activity of 1 <SecNet> instance 
and 1 <Cracker> instance. 
"""
class SecEnv:

    def __init__(self,sn,cracker):
        assert type(sn) == SecNet
        assert type(cracker) == Cracker

        self.sn = sn
        self.cracker = cracker
        self.td = TDBridge()
        return

    """
    fetch the subgraph for the <TDir> of 
    <TDirector>.
    """
    def tdbridge_op__pass_G_to_C(self,cidn):
        if cidn not in self.sn.occ_crackl:
            return None


        c = self.sn.occ_crackl[cidn][0] 
        assert type(c.td) != type(None)

        self.td.pass_info__G(c,self.sn)
        return

    def load_crackling(self):

        return -1 
from cracker import * 

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

    def __init__(self,sn,crck):
        assert type(sn) == SecNet
        assert type(crck) == Cracker

        self.sn = sn
        self.crck = crck 
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

    def run_agent(self,is_sn:bool=True):

        if is_sn:
            q = self.sn
        else: 
            q = self.cracker

        return -1

    def run_cracker(self):

        return -1

    def load_crackling(self):
        return -1 

def SecEnv_sample_1():
    sn3 = SecNet_sample_C3()
    bi = BackgroundInfo.generate_instance(\
            sn3.irc,sn3.srm)
    i2hm = BackgroundInfo.naive_IRC2HypStruct_map(sn3.irc,full_hypseq=False,\
            naive_split=2)

    crck = Cracker(i2hm,bi,6)

    return SecEnv(sn3,crck) 
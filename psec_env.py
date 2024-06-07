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

    def __init__(self,sn,crck,rnd_struct=random):
        assert type(sn) == SecNet
        assert type(crck) == Cracker

        self.sn = sn
        self.crck = crck 
        self.rnd_struct = rnd_struct
        self.td = TDBridge()
        return

    # TODO: test or delete. 
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
            q = self.crck
        return -1

    def run_cracker(self):
        return -1

    # TODO: test this
    """
    process that handles the instantiation
    of <Crackling>+<TDirector> for the 
    <Cracker>'s next target (set of <IsoRing>s). 
    """
    def instantiate_cracker_target(self):
        s = self.cracker.next_target()
        d = self.isoringset_dim(s)
        dx = [(k,v) for k,v in d.items()]

        # load the cracklings
        self.cracker.load_cracklings_for_secset(dx)
        return -1

    def update_cracklings_to_SecNet(self):
        stat = True
        for i in range(len(self.cracker.cracklings)):
            stat2 = self.update_crackling_to_SecNet(i)

            if not stat: continue
            stat = stat and stat2
        return stat

    def update_crackling_to_SecNet(self,crackling_index):
        crcklng = self.crck.cracklings[crackling_index]
    
        # default set at node 0 
        self.sn.set_crackling(crcklng,0)

        # find an entry point that crcklng accepts
        tds = []
        cidn = crcklng.cidn
        stat_ = False
        for x in self.sn.entry_points:
            tdirector = self.sn.tdirector_instance_for_crackling_at_entry_point(\
                cidn,x,radius=self.crck.radar_radius)
        
            stat_ = self.crck.accept_TDirector_at_entry_point(cidn,tdirector)
            if stat_: 
                self.crck.load_TDirector(cidn,tdirector)
                break
            tds.append(tdirector)

        if not stat_:
            if len(tds) == 0:
                return False
            qi = rnd_struct.randrange(0,len(tds))
            self.crck.load_TDirector(cidn,tds[qi])
        return True

def SecEnv_sample_1():
    sn3 = SecNet_sample_C3()
    bi = BackgroundInfo.generate_instance(\
            sn3.irc,sn3.srm)
    i2hm = BackgroundInfo.naive_IRC2HypStruct_map(sn3.irc,full_hypseq=False,\
            naive_split=2)

    crck = Cracker(i2hm,bi,6)

    return SecEnv(sn3,crck) 
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

    def pickle_thyself(self,cpath,spath):
        self.crck.pickle_thyself(cpath)
        self.sn.pickle_thyself(spath)

    # TODO: slow process, needs improvements. 
    @staticmethod
    def unpickle_thyself(spath,sargs,cpath,rnd_struct):
        assert len(sargs) == 3
        print("unpickling <Cracker>")
        crck = Cracker.unpickle_thyself(cpath)

        print("unpickling <SecNet>")
        sn = SecNet.unpickle_thyself(spath,\
            sargs[0],sargs[1],sargs[2])
        return SecEnv(sn,crck,rnd_struct)

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

    ############### methods for instantiating and running
    ############### <Crackling> agents. 

    def run(self,timespan=1.0):
        self.run_agent(True)
        self.run_agent(False)
        return

    def run_agent(self,is_sn:bool=True):
        if is_sn:
            q = self.sn
        else: 
            q = self.crck
        return -1

    



    """
    process that handles the instantiation
    of <Crackling>+<TDirector> for the 
    <Cracker>'s next target (set of <IsoRing>s). 
    """
    def instantiate_cracker_target(self):
        s = self.crck.next_target()
        d = self.sn.isoringset_dim(s)
        dx = [(k,v) for k,v in d.items()]

        # load the cracklings
        self.crck.load_cracklings_for_secset(dx)
        self.update_cracklings_to_SecNet()
        return 

    # TODO: 
    def run_cracklings(self):

        # iterate through each crackling in 
        return -1

    # TODO: 
    def run_crackling(self,crackling_index):
        crckl = self.crck.cracklings[crackling_index]
        return -1 

    """
    pass graph info to the Crackling tdirector
    """
    def update_crackling_graph(self,crackling_index):
        crckl = self.crck.cracklings[crackling_index]
        tdr = crckl.fetch_tdir()
        sg = self.sn.subgraph_for_TDir(tdr)
        crckl.td.load_graph(sngc)

    ########################################################

    """
    update the <Crackling> instances declared by 
    <Cracker> to <SecNet> var.
    """
    def update_cracklings_to_SecNet(self):
        stat = True
        for i in range(len(self.crck.cracklings)):
            stat2 = self.update_crackling_to_SecNet(i)

            if not stat: continue
            stat = stat and stat2
        return stat

    """
    pre-requisite method for above method. 
    """
    def update_crackling_to_SecNet(self,crackling_index):
        crcklng = self.crck.cracklings[crackling_index]
    
        # find an entry point that crcklng accepts
        tds = []
        cidn = crcklng.cidn
        stat_ = False
        for x in self.sn.entry_points:
            # set crackling at entry point
            self.sn.set_crackling(crcklng,x)
            tdirector = self.sn.tdirector_instance_for_crackling_at_entry_point(\
                cidn,x,radius=self.crck.radar_radius)
            print("X: ",x)
            print(tdirector.loc())
            print()
            stat_ = self.crck.accept_TDirector_at_entry_point(cidn,tdirector)
            if stat_: 
                tdirector.switch_obj_stat()
                self.crck.load_TDirector(cidn,tdirector)
                break
            tds.append(tdirector)

        if not stat_:
            print("NO ENTRY POINT; choose random.")
            if len(tds) == 0:
                return False
            qi = self.rnd_struct.randrange(0,len(tds))

            print("TDS SELECTION")
            xs = tds[qi]
            print(xs)
            print("location: ",xs.loc())

            self.crck.load_TDirector(cidn,tds[qi])
            self.sn.set_crackling(crcklng,tds[qi].loc()) 
        return True

def SecEnv_sample_1(sn3=None):
    if type(sn3) == type(None):
        sn3 = SecNet_sample_C3()

    print("generating background info")
    bi = BackgroundInfo.generate_instance(\
            sn3.irc,sn3.srm)
    i2hm = BackgroundInfo.naive_IRC2HypStruct_map(sn3.irc,full_hypseq=False,\
            naive_split=2)
    crck = Cracker(i2hm,bi,6)
    return SecEnv(sn3,crck) 
from cracker import * 

class TDBridge:

    def __init__(self):
        return 

    def pass_info__G(self,c:Crackling,g:SecNet):
        assert type(c) == Crackling
        assert type(g) == SecNet
        tdr = c.td.tdir
        G = g.subgraph_for_TDir(tdr)
        c.td.load_graph(G)
        return

# NOTE: no error-handling for arguments.
class Colocks:

    """
    every element of sets `iset`,`cset` is a stringized
    representation of a co-location instance in a <SecEnv>,
    of the form 
        [0] crackling idn
        [1] isoring idn
        [2] node colocation
    """
    def __init__(self,iset,cset):
        assert type(iset) == set
        assert type(cset) == set

        self.iset = iset
        self.cset = cset
        return

    def update_i(self,iset):
        assert type(iset) == set
        self.iset = iset
        return

    def update_c(self,cset):
        assert type(cset) == set
        self.cset = cset 
        return

    def agent_coloc_by_seti(self,aidn,is_cidn,is_iset:bool):
        s = set()
        q = self.iset if is_iset else self.cset 
        c = 0 if is_cidn else 1
        c2 = (c + 1) % 2
        for q_ in q:
            q2 = Colocks.parse_coloc_str(q_)
            if q2[c] == aidn:
                s = s | {q2[c2]}
        return s 

# TODO: check Colocks


    """
    return:
    - 0 for no co-loc with `agent_idn` OR
      1 for co-loc (interdiction) with `agent_idn` OR
      2 for co-loc (cracking) with `agent_idn`. 
    """
    def cstat(self,agent_idn,is_isoring:bool):
        target_index = 0 if not is_isoring else 1 

        for ix in self.iset:
            q = Colocks.parse_coloc_str(ix)
            if q[target_index] == agent_idn:
                return 1

        for ix in self.cset:
            q = Colocks.parse_coloc_str(ix)
            if q[target_index] == agent_idn:
                return 2
        return 0 
    

    @staticmethod
    def parse_coloc_str(s):
        s_ = matrix_methods.string_to_vector(s,int)
        assert len(s_) == 3
        return s_

"""
environment for the dual-agent 
activity of 1 <SecNet> instance 
and 1 <Cracker> instance. 
"""
class SecEnv:

    """
    sn := SecNet
    crck := Cracker
    ct_ratio := (maximum number of cracking attempts by
                <Crackling>) / (time := 1)
    """
    def __init__(self,sn,crck,rnd_struct=random,\
        ct_ratio=5000,vb=0):
        assert type(sn) == SecNet
        assert type(crck) == Cracker
        assert ct_ratio >= 1
        assert type(ct_ratio) == int 
        assert vb >= 0

        self.sn = sn
        self.crck = crck 
        self.rnd_struct = rnd_struct
        self.ct_ratio = ct_ratio
        self.verbose = vb

        self.td = TDBridge() 
        # active <CBridge> instances
        self.cbs = [] 
        self.cbs_idn_ctr = 0

        # coloc info
        self.ci = Colocks(set(),set())
        self.coloc_register()
        self.coloc_leak_update()
        return

    def preprocess(self):
        self.load_IRCLD_into_SecNet()
        self.bi_update_to_cracker_hyp() 
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

    # TODO: test 
    def run(self,timespan=1.0):
        if self.verbose: 
            print("***************** ONE RUN,T={}".format(timespan))
            print()

        self.cproc(timespan)
        self.iproc(timespan)
        self.postmove_update() 

        if self.verbose: 
            print("\n********************************")
            print()

        return

    """
    main process for the <Cracker> of <SecEnv>; it is 
    the decision-chain function that <Cracker> uses.
    """
    def cproc(self,timespan=1.0):
        # check status
        cd = self.crck.cstat()
        vs = set(cd.values())

        # case: load new cracking targets
        if vs == {2}:
            stat = self.instantiate_cracker_target()
            if not stat: return False
            return self.cproc()

        # proceed to run all <Crackling> instances
        if self.verbose: 
            print("------------ CPROC\n")

        for i in range(len(self.crck.cracklings)):
            self.cproc_(i,timespan)
        return True

    def cproc_(self,index,timespan=1.0):
        c = self.crck.cracklings[index]
        s = self.crck.crackling_stat(index) 

        if self.verbose:
            print("-- CSTAT ON {}: {}".format(c.cidn,s))

        # done
        if s == 2: 
            return False

        # continue cracking
        if s == 1:
            if self.verbose: 
                print("-- CRACKING VIA CBRIDGE")

            cidn = c.cidn 
            br = self.fetch_bridge(cidn,True)

            # case: no existing bridge, declare new bridge
            if type(br) == type(None): 
                ai = self.ci.agent_coloc_by_seti(cidn,True,False)
                assert len(ai) == 1
                self.make_bridge(cidn,ai.pop(),ssih=5)

            iterations = int(round(self.ct_ratio * timespan))
            return self.run_CBridge(iterations,cidn)

        # interdiction
        if s == 0:
            if self.verbose:
                print("-- INTERDICTION")
            return True

        # TODO: review dec.
        c.td_next(timespan,set_roam=True,\
            verbose=self.verbose)

    """
    process that handles the instantiation
    of <Crackling>+<TDirector> for the 
    <Cracker>'s next target (set of <IsoRing>s). 
    """
    def instantiate_cracker_target(self):
        s = self.crck.next_target()
        if len(s) == 0: return False

        d = self.sn.isoringset_dim(s)
        dx = [(k,v) for k,v in d.items()]

        # load the cracklings
        self.crck.load_cracklings_for_secset(dx)
        self.update_cracklings_to_SecNet()
        return True

    """
    pass graph info to the Crackling tdirector
    """
    def update_crackling_graph(self,crackling_index):
        crckl = self.crck.cracklings[crackling_index]
        tdr = crckl.fetch_tdir()
        sg = self.sn.subgraph_for_TDir(tdr)
        crckl.td.load_graph(sngc)

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
                self.sn.set_crackling(crcklng,tdirector.loc()) 

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

    # TODO: use or delete
    def load_subgraphs_for_Cracklings(self):
        for q in self.crck.cracklings:
            td = q.fetch_td()
            assert type(td) == TDir 
            sng = self.sn.subgraph_for_TDir(td)
            q.td.load_graph(sng)
        return

    ##################### methods to handle <IsoRing> decisions


    def iproc(self,timespan = 1.0):
        if self.verbose:
            print("------------ IPROC")
        q = len(self.sn.irc.irl)

        for i in range(q):
            self.iproc_(i,timespan)
        return

    # TODO: add more code and test. 
    def iproc_(self,iindex,timespan = 1.0):

        ir = self.sn.irc.irl[iindex]
        cstat = self.ci.cstat(ir.sec.idn_tag,True)
        
        if self.verbose:
            print("I={},L={},STAT={}".format(ir.sec.idn_tag,ir.td.loc(),\
                cstat))
            print("S.G. KEYS")
            print(ir.td.resource_sg.d.keys())
            print("CLOCS")
            print(ir.td.resource_sg.crackling_locs)
        
        # case: is being cracked --> immobilized.
        if cstat == 2:
            if self.verbose: 
                print("CRACKING")
            return

        v = bool(self.verbose)
        ir.default_secproc(timespan,self.rnd_struct,v) 

    ############ TODO: methods to handle <CBridge>s.
    ########################################################

    """
    Instantiates a <CBridge> instance b/t <Crackling> 
    of identifier `cidn` + <Isoring> of identifier `iidn`. 

    arguments:
    - cidn := identifier for <CBridge> 
    - iidn := identifier for <IsoRing>
    - ssih := hop value for <RSSI>.
    """
    def make_bridge(self,cidn,iidn,ssih=5):
        c = self.crck.fetch_crackling(cidn)
        i = self.sn.irc.fetch_IsoRing(iidn)
        hs = c.hs

        cb = CBridge(c,i,hs,ssih,self.cbs_idn_ctr)
        cb.verbose = self.verbose

        self.cbs_idn_ctr += 1
        self.cbs.append(cb) 
        return

    """
    Fetches the <CBridge> holding an agent of 
    identifier `idn`. 

    return:
    - <CBridge>|None
    """
    def fetch_bridge(self,idn,is_cidn:bool=True,\
        allb:bool=False):
        index = 0 if is_cidn else 1 
        q = []
        print("# BRIDGES: ",len(self.cbs))
        for cb in self.cbs:
            idnx = cb.agent_idns()[index]
            if idnx == idn: 
                if not allb: return cb
                q.append(cb)

        if not allb:  
            return None
        return q 

    def fetch_bridge_(self,cidn,iidn): 

        for cb in self.cbs:
            idnx = cb.agent_idns()
            if idnx[0] != cidn:
                continue
            if idnx[1] != iidn:
                continue
            return cb 
        return None

    """
    conducts cracking operation on <CBridge> 
    for `next_rate` iterations.

    return:
    - bool, ?is not early termination?
    """
    def run_CBridge(self,next_rate,iidn):

        def rcb(cb_):
            if self.verbose:
                print("-- CBRIDGE OP: {}".format(cb_.agent_idns()))
                print("\t-/-/-/-/")
            for i in range(next_rate):
                print("ITER=",i)
                try:
                    qc = next(cb_)
                except:
                    return False

                stat = type(qc) != type(None)
                if not stat: return False
            return True 

        cb = self.fetch_bridge(iidn,False,True)
        print("FETCHING BRIDGE")
        print(cb)

        statvec = []
        for cb_ in cb:
            v = rcb(cb_)
            statvec.append(v)
        return statvec 

    #####################################
    ############## <TDirector> instantiation for
    ############## each <IsoRing>|<Crackling>
    ############## agent. 
    #####################################

    def instantiate_td_for_IRC(self,rd,td):
        self.sn.load_TDirectors_for_IRC(rd,td)
        self.load_subgraphs_for_IRC() 
        return

    def load_subgraphs_for_IRC(self):
        for q in self.sn.irc.irl:
            td = q.fetch_td()
            assert type(td) == TDir 
            sng = self.sn.subgraph_for_TDir(td)
            q.td.load_graph(sng)
        return

    #####################################
    ############### post-move update
    #####################################

    """
    updates the locations of <Crackling>+<IsoRing>
    agents onto the variables of <SecNet>, as well
    as each agent's <TDirector> graph. 
    """
    def postmove_update(self):
        self.sn.update_occ_crackl()
        self.sn.update_nla() 

        self.coloc_register()
        self.coloc_leak_update()
        return

    # TODO: 
    def coloc_register(self):

        iset = set()
        cset = set()

        s = self.coloc() 

        while len(s) > 0:
            s_ = s.pop()
            s2 = Colocks.parse_coloc_str(s_)
            stat = self.sn.is_secnode(s2[2])

            if stat:
                cset = cset | {s_}
            else:
                iset = iset | {s_}

        self.ci.update_c(cset)
        self.ci.update_i(iset)

        self.crackling_stat_update()
        return

    """
    def crackling_stat_update(self,is_cset:bool):
        ws = self.ci.cset if is_cset else \
                self.ci.iset

        for ws_ in ws:
            s2 = Colocks.parse_coloc_str(ws_)
            cr = self.crck.fetch_crackling(s2[0]) 
            assert type(cr) != type(None)

            if is_cset:
                cr.cstat = True
            else:
                cr.istat = True
    """


    def crackling_stat_update(self): 

        for c in self.crck.cracklings:
            q = c.td.coloc()
            c.cstat,c.istat = False,False 
            if len(q) > 0:
                if self.sn.is_secnode(c.loc()):
                    c.cstat = True
                else:
                    c.istat = True
        return

    #############

    """
    return:
    - set::(crackling idn,isoring idn,node colocation)
    """
    def coloc(self):
        # iterate through each <Crackling>
        # and check status
        dx = set()
        for c in self.crck.cracklings:
            q = c.td.loc() 
            target = c.target_of_HypStruct()
            stat = target in c.td.resource_sg.ring_locs
            q2 = None 
            if stat:
                q2 = c.td.resource_sg.ring_locs[target]

            if q == q2: 
                a = np.array([c.cidn,target,q2])
                s = matrix_methods.vector_to_string(a,int)
                dx = dx | {s}
        return dx

    # TODO: 
    def coloc_leak_update(self):
        if self.verbose: 
            print("--- PERFORMING LEAKS : ",len(self.ci.iset))

        for interdict in self.ci.iset:
            self.leak_by_str_idn(interdict)
        return


    ################ methods for leaking

    # TODO: test 
    def bi_update_to_cracker_hyp(self):
        bi = self.crck.bi
        assert type(bi) == BackgroundInfo
        bi.apply_leakmap_on_IRC2HypStruct_map(\
            self.crck.hyp_map) 

    # TODO: test 
    def load_IRCLD_into_SecNet(self):
        self.sn.irc.load_default_IRCLD()

    # TODO: test
    def leak_by_str_idn(self,sidn):
        psidn = Colocks.parse_coloc_str(sidn)
        crckling = self.crck.fetch_crackling(psidn[0])

        assert type(crckling) == Crackling
        assert type(crckling.hs) == HypStruct

        # fetch the appropriate leak
        sec_idn = psidn[1]

        ir = self.sn.irc.fetch_IsoRing(sec_idn)
        opt_dim = ir.secdim_seq()[ir.repi]
        L = self.sn.irc.ircld.fetch_Leak(sec_idn,opt_dim)
        outp = L.leak_info(ir)

        # case: no more leaks
        if type(L.prev_li) == type(None):
            return None

        # case: update <HypStruct> by latest <LeakInfo>
            # fetch the <HypStruct>

            # retrieve the `leak_value` and `leak_idn`
        lv = L.prev_li[2]
        fidn = L.prev_li[1] 
        
        hs_ = HypInfer.infer_FX(\
            deepcopy(crckling.hs),lv,fidn)
        crckling.hsi = hs_ 
        return hs_ 

    ################ end-round update 

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

def SecEnv_sample_2():
    sn = SecNet_sample_CSmall()
    return SecEnv_sample_1(sn)
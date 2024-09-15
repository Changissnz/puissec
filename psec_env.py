from cracker import * 

# TODO: unused, possibly delete. 
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

    def __str__(self):
        s1 = str(self.iset) + "\n"
        s1 += str(self.cset) + "\n\n"
        s1 += "----------------------"
        return s1 

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
    vb := int, verbosity in {0,1,2}
    mode_open_info := (0|1|2,0|1|2);  
        [0] corresponds to <SecNet>,
        [1] corresponds to <Cracker>. 

    """
    def __init__(self,sn,crck,rnd_struct=random,\
        ct_ratio=5000,vb=0,mode_open_info = (0,0)):
        assert type(sn) == SecNet
        assert type(crck) == Cracker
        assert ct_ratio >= 1
        assert type(ct_ratio) == int 
        assert vb >= 0
        assert type(mode_open_info) == tuple and len(mode_open_info) == 2
        assert mode_open_info[0] in {0,1,2} and \
            mode_open_info[1] in {0,1,2} 

        self.sn = sn
        self.crck = crck 
        self.rnd_struct = rnd_struct
        self.ct_ratio = ct_ratio
        self.verbose = vb

        self.td = TDBridge() 
        # active <CBridge> instances
        self.cbs = [] 
        self.cbs_idn_ctr = 0

        self.icrack = {}

        # coloc info
        self.ci = Colocks(set(),set())
        self.coloc_register()
        self.coloc_leak_update()
        self.moi = mode_open_info

        self.tcum = 0.0 
        return

    def preprocess(self):
        self.load_IRCLD_into_SecNet()
        self.bi_update_to_cracker_hyp() 

        self.instantiate_cracker_target()
        self.instantiate_td_for_IRC(5,1.0)

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

    # TODO: test 
    """
    passes information between <IsoRing>s and 
    <Crackling>s based on `self.moi`.
    """
    def pass_info(self,is_IsoRing:bool):
        if self.moi == (0,0):
            return

        q = self.moi[0] if is_IsoRing else self.moi[1]
        # pass IsoRing info to each <Crackling>
        if is_IsoRing: 
            for c in self.crck.cracklings:
                # case: target IsoRing not found
                if c.td.obj_stat == "search for target":
                    continue
                i = c.td.td.target_node
                c_ = c.cidn 
                self.pass_info_(i,c_,is_IsoRing,q)
        else:
            for ir in self.sn.irl:
                # case: no cracklings in sight
                if ir.td.obj_stat == "radar null":
                    continue

                # pass info from all visible <IsoRing>
                qs = ir.check_radar(True)

                for q_ in qs:
                    self.pass_info_(q_,ir.sec.idn_tag,is_IsoRing,q)
        return

    """
    for an <IsoRing> identified by `i`, and a 
    <Crackling> identified by `c`, method
    passes information by one of:
    i-->c,
    c-->i,
    according to `is_IsoRing`. The information 
    passed is denoted by `info_type` in {0,1,2}:
    
    arguments:
    - i := str|int, <IsoRing> idn
    - c := str|int, <Crackling> idn
    - is_IsoRing := bool, specifies directionality of info-pass.
    - info_type := int, one of 0|1|2,
        0 := no info passed,
        1 := velocity passed,
        2 := node location passed. 

    return:
    - bool, ?is info by type passed?
    """
    # TODO: test 
    def pass_info_(self,i,c,is_IsoRing:bool,info_type:int):
        assert info_type in {0,1,2}
        if info_type == 0: return True

        # fetch the two agents `i`,`c`. 
        i_,c_ = i,c 
        i = self.fetch_subagent(i,True)
        c = self.fetch_subagent(c,False)

        if type(i) == type(None) or type(c) == type(None): return False

        # case: I passes to C  
        if is_IsoRing:
            v = i.td.td.velocity if info_type == 1 else i.td.td.location
            c.recv_open_info(info_type,i_,v)
            return True
            
        v = c.td.td.velocity if info_type == 1 else c.td.td.location
        i.recv_open_info(info_type,c_,v)
        return True

    ############### methods for instantiating and running
    ############### <Crackling> agents. 

    # TODO: test 
    def run(self,timespan=1.0):
        if self.verbose: 
            print("***************** ONE RUN,T={},T_CUM={}".format(timespan,self.tcum))
            print()

        self.run_(timespan)
        self.postmove_update() 

        if self.verbose: 
            print("\n********************************")
            print()

        self.tcum += timespan
        return

    """
    auxiliary method for the above method `run` that 
    takes into account `open info`. 
    """
    def run_(self,timespan):
        x0,x1 = None,None
        b = False

        # base case: no open info
        if self.moi == (0,0): 
            x0,x1 = self.cproc,self.iproc
            b = True
        elif self.moi[0] > self.moi[1]:
            x0,x1 = self.cproc,self.iproc
        elif self.moi[0] < self.moi[1]:
            x0,x1 = self.iproc,self.cproc
            b = True
        else: 
            q = self.rnd_struct.uniform(0.0,1.0)
            if q >= 0.5:
                x0,x1 = self.iproc,self.cproc 
                b = True
            else: 
                x0,x1 = self.cproc, self.iproc

        # move the first
        x0(timespan)

        # ?pass info?
        self.pass_info(b)

        # move the second
        x1(timespan)
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
            print("NEW TARGET")
            stat = self.instantiate_cracker_target()
            if not stat: return False
            return self.cproc()

        # proceed to run all <Crackling> instances
        if self.verbose: 
            print("------------ CPROC\n")
            print("* active <Cracklings>: ",len(self.crck.cracklings))

        for i in range(len(self.crck.cracklings)):
            self.cproc_(i,timespan)
        return True

    def cproc_(self,index,timespan=1.0):
        c = self.crck.cracklings[index]
        s = self.crck.crackling_stat(index) 

        if self.verbose:
            print("-- CSTAT ON {}: {}".format(c.cidn,s))

        print()
        # done
        if s == 2: 
            print("DONE")
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
                self.make_bridge(cidn,ai.pop())

            iterations = int(round(self.ct_ratio * timespan))
            return self.run_CBridge(iterations,cidn)

        # interdiction
        if s == 0:
            if self.verbose:
                print("-- INTERDICTION")
            return True

        # TODO: review dec.
        vm = c.td_next(timespan,\
            verbose=self.verbose)
        self.crck.energy -= vm 

    """
    process that handles the instantiation
    of <Crackling>+<TDirector> for the 
    <Cracker>'s next target (set of <IsoRing>s). 
    """
    def instantiate_cracker_target(self):
        s = self.crck.next_target()
        if len(s) == 0: return False

        self.crck.initiated = True 
        d = self.sn.isoringset_dim(s)
        dx = [(k,v) for k,v in d.items()]
        print("INSTANTIATING CRACKER TARGET")
        print(dx) 
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
        print("\n\nUPDATING CRACKLINGS ",len(self.crck.cracklings))
        for i in range(len(self.crck.cracklings)):
            stat2 = self.update_crackling_to_SecNet(i)

            if not stat: continue
            stat = stat and stat2  
        return stat

    """
    pre-requisite method for above method. 
    """
    def update_crackling_to_SecNet(self,crackling_index):
        print("\n\nUPDATING")
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
            print("ACCEPTING @ ENTRY POINT: ",stat_)
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
        vl = ir.td_next(timespan,self.rnd_struct,v)
        print("VLLL: ",vl)
        self.sn.energy -= vl 

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
    def make_bridge(self,cidn,iidn):
        c = self.crck.fetch_crackling(cidn)
        i = self.sn.irc.fetch_IsoRing(iidn)
        hs = c.hs
        cb = CBridge(c,i,self.cbs_idn_ctr,\
            #2, True)#self.verbose == 2) 
            self.ct_ratio,self.verbose == 2)

        self.cbs_idn_ctr += 1
        self.cbs.append(cb) 
        return

    """
    Fetches the <CBridge> holding an agent of 
    identifier `idn`. 

    return:
    - list(<CBridge>)|None
    """
    def fetch_bridge(self,idn,is_cidn:bool=True,\
        allb:bool=False):
        index = 0 if is_cidn else 1 
        q = []
        print("# BRIDGES: ",len(self.cbs))
        print("IDN: ", idn, " CIDN: ",is_cidn)
        for cb in self.cbs:
            idnx = cb.agent_idns()
            print("AGENT IDNS: ",idnx)
            idnx = idnx[index]
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
    - list(<bool,int>), 
        [0] ?is not early termination?
        [1] number of iterations
    """
    def run_CBridge(self,next_rate,iidn):
        
        def rcb(cb_):
            if self.verbose:
                print("-- CBRIDGE OP: {}".format(cb_.agent_idns()))
                print("\t-/-/-/-/")
            j = 0
            for i in range(next_rate):
                if self.verbose: 
                    print("ITER=",i)
                                
                ### NOTE: new
                qc = next(cb_) 
                if type(qc) == type(None):
                    if not cb_.crackling.astat:
                        cb_.cfail = True
                        cb_.crackling.fstat = True
                        return False, j 

                j = i + 1
                if cb_.crackling.astat:
                    # add the result to <Cracker.csoln>
                    ks = list(list(cb_.crackling.cracked_dict.items())[0])
                    ks[0] = matrix_methods.string_to_vector(ks[0],float)
                    sol = (cb_.crackling.hs.seq_idn,\
                        cb_.crackling.hs.target_index,\
                        ks[0],ks[1])
                    self.crck.csoln += sol
                    self.crck.oopi += 1
                    return False, j  

            return True,next_rate 

        cb = self.fetch_bridge(iidn,True,True)
        print("FETCHING BRIDGE")
        print(cb)

        statvec = []
        for cb_ in cb:
            v = rcb(cb_)
            statvec.append(v)

            self.crck.energy -= v[1]

        self.remove_spent_CBridge() 
        return statvec

    def remove_spent_CBridge(self):
        cbs2 = []
        for cb in self.cbs:
            if cb.cfail == False:
                cbs2.append(cb)
        self.cbs = cbs2
            

    #####################################
    ############## <TDirector> instantiation for
    ############## each <IsoRing>|<Crackling>
    ############## agent. 
    #####################################

    """
    loads a <TDirector> instance for each <IsoRing>
    """
    def instantiate_td_for_IRC(self,rd,td):
        self.sn.load_TDirectors_for_IRC(rd,td)
        self.load_subgraphs_for_IRC() 
        return

    """
    loads the respective graph-in-sight for each
    <IsoRing> 
    """
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

        # clear open information from last timespan
        self.sn.clear_open_info()
        self.crck.clear_open_info() 

        # clear fstat cracklings
        self.crck.remove_cracklings__fstat()
        return

    """
    registers the co-locations of (<IsoRing>,<Crackling>)
    pairs based on SEC status of node:
    - stat := cracking <--> SEC
    - stat := interdiction <--> NSEC 
    """
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
    updates the (cstat,istat) of each <Crackling>
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

    def coloc_leak_update(self):
        if self.verbose: 
            print("--- PERFORMING LEAKS : ",len(self.ci.iset))

        for interdict in self.ci.iset:
            q = self.leak_by_str_idn(interdict)

            # case: cracked 
            if type(q) == tuple:
                if self.verbose: 
                    print("\t\t** CRACK BY INTERDICTION")
                assert len(q) == 2

                psidn = Colocks.parse_coloc_str(interdict)

                # check to see if already cracked
                ir = self.sn.irc.fetch_IsoRing(psidn[1])
                if ir.sec.idn_tag in self.icrack:
                    d = ir.sec.dim() 
                    if d in self.icrack[ir.sec.idn_tag]:
                        if self.verbose == 2:
                            print("ALREADY CRACKED. REMOVING SPENT.")
                        self.crck.remove_spent_crackling(psidn[0])
                        continue 

                crckling = self.crck.fetch_crackling(psidn[0])
                ti = crckling.hs.target_index
                self.transfer_V(interdict,ti,q[0],q[1])                
        return

    def fetch_subagent(self,agent_idn,is_IsoRing:bool):
        if not is_IsoRing:
            return self.crck.fetch_crackling(agent_idn)
        return self.sn.fetch_IsoRing(agent_idn)

    def transfer_V(self,sidn,opt_index,svec,pr_score):
        psidn = Colocks.parse_coloc_str(sidn)
        sec_idn = psidn[1]
        ir = self.sn.irc.fetch_IsoRing(sec_idn)
        opt_dim = ir.secdim_seq()[ir.repi]
        self.crck.csoln = self.crck.csoln + (sec_idn,opt_index,svec,pr_score)
        if psidn[0] not in self.icrack:
            self.icrack[psidn[0]] = {}
        self.icrack[psidn[0]][opt_dim] = True 

    ################ methods for leaking

    # TODO: test OR delete. 
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

        if self.verbose == 2:
            print("OUTP: ",outp)

        # case: no more leaks
        if type(L.prev_li) == type(None):
            prv = ir.sec.seq_pr() 
            return (ir.sec.seq,prv) 

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

# TODO: complete|delete.
"""
irc_args := (number of <Sec> sources,singleton_range,\
    dimension_range,num_optima_range,optima_countermeasure_range)
sn_args := (sec node count,nsec node count,num entry points,rnd_struct,path-in-mem size,*sngs_args)
"""
def generate_SecEnv(irc_args,sn_args):
    sn = SecNet.generate(irc_args,sn_args)
    return -1 
from collections import defaultdict 
from bridge import * 
from struct_samples import *

class IsoRingedChain:

    def __init__(self,ss:SecSeq,bound,\
        rnd_struct,rs_seed:int):
        self.ss = ss 
        self.irl = []
        self.ircld = None 
        self.load_IsoRings(bound,rnd_struct,rs_seed)
        self.rnd_struct = rnd_struct

    def explode(self,osl:int,num_blooms:int):
        
        for i in self.irl:
            self.explode_IsoRing(i.sec.idn_tag,osl,num_blooms)
        return

    def explode_IsoRing(self,sec_idn,osl:int,num_blooms:int):
        ir = self.fetch_IsoRing(sec_idn)
        assert type(ir) != type(None)
        ir.explode_contents(osl,num_blooms,self.rnd_struct)
        return

    def fetch_IsoRing(self,sec_idn):
        for x in self.irl:
            if x.sec.idn_tag == sec_idn:
                return x
        return None

    def revert_to_SecSeq(self):
        sx = []
        for q in self.irl:
            sx.append(q.sec) 
        ss = SecSeq(sx)
        return ss 

    """
    default args. are 
        - spacing_ratio_range := 
        - outlier_pr_ratio := ()
        - num_bounds := min(100,|sec.opm| * 2)
        - rnd_seed := 113 

    Each <IsoRing> uses the same default <CVec> generator 
    scheme. 
    """
    @staticmethod
    def default_Sec2IsoRing(sec:Sec,bound,rnd_struct,rs_seed):

        if type(rs_seed) != type(None): 
            rnd_struct.seed(rs_seed)
        else: 
            rnd_struct.seed(113)
        
        # generate the BOF
        spacing_ratio_range = rnd_struct.uniform(0.1,0.8)
        spacing_ratio_range = [spacing_ratio_range,rnd_struct.uniform(spacing_ratio_range,0.8)]
        outlier_pr_ratio = rnd_struct.uniform(0.4,0.6) 

        num_bounds = min([100,2 * len(sec.opm)])
        bof = BoundedObjFunc.generate_BoundedObjFunc(\
            deepcopy(bound),spacing_ratio_range,\
            outlier_pr_ratio,num_bounds,rs_seed)
        return IsoRing(sec,bof,\
            bound,cvecl=None)

    def __len__(self): 
        return len(self.irl)

    def __getitem__(self,i):
        assert i <= len(self.irl)
        assert i >= 0 
        return self.irl[i] 

    def load_IsoRings(self,singleton_bound=DEFAULT_SINGLETON_RANGE,\
        rnd_struct=random,rs_seed=8): 

        if len(self.ss.sequence) == 0:
            self.ss = None
            return 

        self.irl = [] 
        if type(rs_seed) != type(None):
            rnd_struct.seed(rs_seed)

        for s in self.ss.sequence:
            # get the dim. of s
            q = s.dim() 
            bound_ = np.ones((q,2)) * singleton_bound

            sd = rnd_struct.randrange(0,500)
            ir = IsoRingedChain.default_Sec2IsoRing(s,\
                bound_,rnd_struct,sd) 
            self.irl.append(ir) 
        self.ss = None 

    # TODO: test
    def load_default_IRCLD(self):
        self.ircld = IRCLeakDict(DEFAULT_LEAKSZ_RANGE,\
            self.rnd_struct)    

        for i in self.irl:
            self.ircld.load_default_Leak_for_IsoRing(i)
        return

    """
    fp := folder path
    """
    def pickle_thyself(self,fp):
        assert fp not in {"","/","\\"}

        clear_existing_dir(fp)
        if not os.path.isdir(fp):
            os.mkdir(fp)

        for x in self.irl:
            fp_ = fp + "/node_{}".format(x.sec.idn_tag)
            x.pickle_thyself(fp_)
    
    @staticmethod
    def unpickle_thyself(fp,rnd_struct,rs_seed):

        q = os.listdir(fp)
        ls = [] 
        for q_ in q:
            if q_ in {"excessory","sngc"}: continue

            fp_ = fp + "/" + q_
            ir = IsoRing.unpickle_thyself(fp_) 
            ls.append(ir) 

        # sort all elements by idn
        ls = sorted(ls,key=lambda x: x.sec.idn_tag) 

        ss = SecSeq([])
        irc = IsoRingedChain(ss,None,rnd_struct,rs_seed)
        irc.irl = ls 
        return irc

"""
A graph structure that serves as an 
environment for activity programmed 
in other structs. 
"""
class SecNet:

    def __init__(self,irc,G,sec_nodeset,\
        node_loc_assignment= None,entry_points=3,\
        bound=DEFAULT_SINGLETON_RANGE,\
        rnd_struct=random,rnd_seed=9,\
        path_size=10,sngc=None,energy=1000.0):
        
        assert len(irc) > 0 and type(irc) == SecSeq
        assert len(G) >= len(irc), "len G {} len IRC {}".format(len(G),len(irc))

        self.ss = irc
        self.irc = IsoRingedChain(irc,bound,rnd_struct,rnd_seed) 
        self.rnd_struct = rnd_struct 
        self.path_size = path_size 
        self.d = G 
        self.sec_nodeset = sec_nodeset
        
        # sec idn -> node location 
        self.node_loc_assignment = node_loc_assignment
        self.assign_entry(entry_points)
        if type(self.node_loc_assignment) == type(None):
            self.assign_loc()
        else:             
            assert type(self.node_loc_assignment) == dict
            assert max(self.node_loc_assignment.keys()) == len(self.irc) - 1 
            assert min(self.node_loc_assignment.keys()) == 0 
            assert len(self.node_loc_assignment) == len(self.irc) 
            assert len(self.node_loc_assignment) == len(set(self.node_loc_assignment.values()))

        self.rnd_struct = rnd_struct
        self.srm = self.load_srm()
        self.ss = None
        self.sgc = None
        # occupied cracklings
        self.occ_crackl = {} # Crackling cidn-> [<Crackling>,node location]

        if type(sngc) == type(None):
            self.preprocess_shortest_paths() 
        else:
            assert type(sngc) == SNGraphContainer
            self.sgc = sngc 
        self.energy = NerG(energy)
        self.terminated = self.is_terminated()
        return

    ######################## location update and
    ######################## status 

    def update_occ_crackl(self):
        ks = []
        for k in self.occ_crackl.keys():
            v = self.occ_crackl[k] 
            if v[0].fstat:
                ks.append(k)
                continue

            qv = (v[0],v[0].loc())
            self.occ_crackl[k] = qv 

            # update TDirector if existing
            if type(v[0].td) == type(None):
                continue

            tdx = v[0].td.td
            sgx = self.subgraph_for_TDir(tdx)
            v[0].td.load_graph(sgx)

        # remove fstat <Crackling>s
        for k in ks:
            del self.occ_crackl[k]
        return 

    ## TODO: refactoradorii above function w/ this. 
    def update_nla(self):
        for k in self.node_loc_assignment.keys():
            s = self.irc.fetch_IsoRing(k) 
            self.node_loc_assignment[k] = s.loc() 
            if type(s.td) == type(None):
                continue
            tdx = s.td.td
            sgx = self.subgraph_for_TDir(tdx)
            s.td.load_graph(sgx)
        return

    def is_secnode(self,n):
        return n in self.sec_nodeset
    
    ##############################################

    @staticmethod
    def alt_instantiate(irc,G,sec_nodeset,\
        node_loc_assignment,srm,entry_points,\
        bound=DEFAULT_SINGLETON_RANGE,\
        rnd_struct=random,rnd_seed=9,sngc=None):

        assert type(irc) == IsoRingedChain
        assert type(srm) == SRefMap 
        ss = irc.revert_to_SecSeq()        
        sn = SecNet(ss,G,sec_nodeset,node_loc_assignment,\
            entry_points=1,bound=bound,\
            rnd_struct=rnd_struct,rnd_seed=rnd_seed,\
            sngc=sngc)
        sn.irc = irc
        sn.entry_points = entry_points
        sn.srm = srm 
        return sn 

    """
    similar to the pickle-pattern for <IsoRingedChain>, except 
    has an extra file called `excessory` that contains a list 
    of objects:
    [0] G, the graph 
    [1] sec nodeset
    [2] dict,sec idn.->node loc
    [3] SRefMap instance
    [4] set of entry points
    """
    def pickle_thyself(self,fp):

        excessory = [deepcopy(self.d),deepcopy(self.sec_nodeset),\
                deepcopy(self.node_loc_assignment),\
                deepcopy(self.srm),deepcopy(self.entry_points)]

        self.irc.pickle_thyself(fp)
        efile = open(fp + "/excessory","wb")
        pickle.dump(excessory,efile)
        efile.close()

        self.sgc.pickle_thyself(fp + "/sngc")
        return

    @staticmethod
    def unpickle_thyself(fp,\
        bound=DEFAULT_SINGLETON_RANGE,rnd_struct=random,rnd_seed=9): 
        assert os.path.exists(fp)

        irc = IsoRingedChain.unpickle_thyself(fp,rnd_struct,rnd_seed)
        efile = open(fp + "/excessory","rb")
        excessory = pickle.load(efile)
        efile.close()

        G = excessory[0]
        sec_nodeset = excessory[1]
        node_loc_assignment = excessory[2]
        srefmap = excessory[3]
        entry_points = excessory[4]

        sngc = None

        if os.path.exists(fp + "/sngc"):
            sngc = SNGraphContainer.unpickle_thyself(fp + "/sngc")
        
        ##print("instantiating w/ rndstruct=",rnd_struct)
        sn = SecNet.alt_instantiate(irc,G,sec_nodeset,\
        node_loc_assignment,srefmap,entry_points,\
        bound=bound,rnd_struct=rnd_struct,\
        rnd_seed=rnd_seed,sngc=sngc)
        return sn 

    ######################## loc-set methods for 
    ######################## <Crackling>,<IsoRing>

    def rc_agent_locs_for_subgraph(self,sgc:SNGraphContainer):
        assert type(sgc) == SNGraphContainer
        k = set(sgc.d.keys())

        nla,oc = {},{}

        for k_,v_ in self.node_loc_assignment.items():
            if v_ in k: 
                nla[k_] = v_
        
        for k_,v_ in self.occ_crackl.items():
            if v_[1] in k:
                q = v_[0].hs.seq_idn
                oc[k_] = (v_[1],q)
        sgc.update_rc_agent_locs(nla,oc)
        return nla,oc 
        
    ######################## graph structure functions 

    def subgraph_for_TDir(self,tdir):
        assert type(self.sgc) != type(None)
        rx = tdir.radius 
        v = tdir.velocity
        ##print("RX: ",rx, v)
        sgc1 = self.sgc.subgraph_by_radius_at_refnode(tdir.\
            location,rx)
        self.rc_agent_locs_for_subgraph(sgc1)
        return sgc1

    def preprocess_shortest_paths(self):
        d2 = defaultdict(set)
        for k in self.d.keys():
            d2[k] = set(self.d[k])

        ring_locs = deepcopy(self.node_loc_assignment)

        ocm = self.occ_crackl_map(set(self.d.keys()))

        self.sgc = SNGraphContainer(d2,deepcopy(self.sec_nodeset),\
            ring_locs,ocm,deepcopy(self.entry_points))
        self.sgc.path_size = self.path_size
        self.sgc.DFSCache_fullproc() 
 
    def to_graphvars(self,dx=None):
        assert type(dx) in {defaultdict,type(None)}
        if type(dx) == type(None): 
            dx = d 

        nla = {} 
        for k in dx.keys():
            nla[k] = self.node_loc_assignment[k] 

        return deepcopy(dx),\
            nla,\
            deepcopy(self.sec_nodeset),\
            deepcopy(self.entry_points)

    def load_srm(self):
        opmn = self.ss.sec_instances_to_supermap("l")
        dms = self.ss.sec_instances_to_supermap("d")
        cdms = self.ss.sec_instances_to_supermap("c")
        return SRefMap(opmn,dms,cdms)

    """
    uses `rnd_struct` to assign each element 
    of `irc` a location on a `sec_nodeset`. 
    """
    def assign_loc(self):
        self.node_loc_assignment = {}

        s = list(deepcopy(self.sec_nodeset))
        i = len(self.irc)
        i2 = 0
        while i > 0:
            if len(s) == 0: 
                s = list(set(self.d.keys()) - self.sec_nodeset)
                if len(s) == 0: raise ValueError 
            
            j = self.rnd_struct.randrange(0,len(s))
            q = s.pop(j)
            self.node_loc_assignment[i2] = q 
            i2 += 1
            i -= 1 
        return 

    def assign_entry(self,entry_points):
        assert type(entry_points) in {int,set}
        q = list(self.d.keys())

        # choose random values
        if type(entry_points) == int:
            assert entry_points > 0 and entry_points <= len(q)
            eps = [] 
            while entry_points > 0: 
                i = self.rnd_struct.randrange(0,len(q))
                eps.append(q.pop(i))
                entry_points -= 1
            self.entry_points = set(eps)
            return

        for ep in entry_points:
            assert ep in q
        self.entry_points = entry_points

    """
    generates an instance of a <SecNet> using 
    a <SecSeq> instance (that may or may not 
    have dep./codep.) along with other variables 
    pertaining to the graph to be generated.
    """
    @staticmethod
    def generate(irc,sec_node_count,\
        nsec_node_count,num_entry,path_size,rnd_struct,*sngs_args):
        q = sec_node_count + nsec_node_count
        assert q >= len(irc), "got {} want {}".format(q,len(irc))
        assert num_entry <= q and num_entry > 0 

        q_ = [i for i in range(q)] 
        rnd_struct.shuffle(q_)

        print("make frame")
        snv = q_[:sec_node_count]
        nsnv = q_[sec_node_count:]
        print("SNV: ", snv)
        print("NSNV: ", nsnv)

        sngs = SecNetGenScheme(*sngs_args)
        snfg = SecNetFrameGen(snv,nsnv,sngs)
        snfg.construct_frame()
        print("after frame")
        print("graph")
        print(snfg.d)
        node_loc_assignment = None
        sn = SecNet(irc,snfg.d,set(snfg.sec_nodevec),\
            node_loc_assignment,entry_points=num_entry,\
            rnd_struct=rnd_struct,path_size=path_size)
        return sn 

    # TODO: test this. 
    """
    NOTE: virtually identical to the above method `generate`, 
          excepts generates a graph from <SecNetFrameGen> 
          using some random struct.

    irc_args := (number of <Sec> sources,singleton_range,\
        dimension_range,num_optima_range,optima_countermeasure_range)
    sn_param_args := (p_s,p_n,rnd_struct,num_entry,
    superbound{Sec2IsoRing conversion}); arguments for <SecNetFrameGen>.
    
    """
    @staticmethod
    def generate__autograph(irc_args,sn_param_args):
 
        sec_list = Sec_list_SEEDTYPE_PythonANDNumpy(irc_args[0],\
            irc_args[1],irc_args[2],irc_args[3],irc_args[4])
        sec_list = SecSeq(sec_list)

        G = SecNetFrameGen.generate_graph__type_prop(len(sec_list),\
            sn_param_args[0],sn_param_args[1],sn_param_args[2])

        assert sn_param_args[3] >= 0.0 and sn_param_args[3] <= 1.0 
        num_entry = int(round(sn_param_args[3] * len(G[0])))
        
        rnd_seed = sn_param_args[2].randrange(0,1000) 

        return SecNet(sec_list,G[0],G[1],entry_points=num_entry,\
            bound=sn_param_args[4],rnd_struct=sn_param_args[2],\
            path_size=10)

    """
    irc_args := (number of <Sec> sources,singleton_range,\
        dimension_range,num_optima_range,optima_countermeasure_range)
    sn_args := (sec node count,nsec node count,num entry points,rnd_struct,path-in-mem size,*sngs_args)
    """
    @staticmethod
    def full_generate(irc_args,sn_args):

        sec_list = Sec_list_SEEDTYPE_PythonANDNumpy(irc_args[0],\
            irc_args[1],irc_args[2],irc_args[3],irc_args[4])
        sec_list = SecSeq(sec_list)
        sn = SecNet.generate(sec_list,sn_args[0],\
            sn_args[1],sn_args[2],sn_args[3],sn_args[4],*sn_args[5])
        return sn 

    #################################################
    ###### methods for routing <Crackling>
    #################################################
    """
    set <Crackling> `c` loaded with a 
    <HypStruct> and <TDir> at location `node`. 
    """
    def set_crackling(self,c,node):
        assert node in self.entry_points
        assert type(c) == Crackling
        assert type(c.hs) == HypStruct
        
        # use the <HypStruct> for `c` to set 
        # <TDir>. 
        self.occ_crackl[c.cidn] = (c,node)
        return

    """
    map info passed to <SNGraphContainer>
    """
    def occ_crackl_map(self,ks):
        d = {} 
        for k,v in self.occ_crackl.items():
            if v[1] in ks:
                target = v[0].target_of_HypStruct()
                d[k] = (v[1],target)
        return d 

    # TODO: test.
    def tdirector_instance_for_crackling_at_entry_point(self,\
        cidn,entry_point,radius=4):

        assert cidn in self.occ_crackl
        assert entry_point in self.entry_points
        assert type(self.occ_crackl[cidn][0].hs) == HypStruct

        target = self.occ_crackl[cidn][0].hs.seq_idn
        td = TDirector(entry_point,target,"C",cidn,radius)

        sngc = self.subgraph_for_TDir(td.td)
        td.load_graph(sngc)
        return td

    def isoringset_dim(self,ir_set):
        dx = {}
        for ir in ir_set:
            irx = self.irc.fetch_IsoRing(ir)
            assert type(irx) != type(None)
            sx = irx.sec_cache[irx.repi]
            dx[ir] = sx.dim()
        return dx 

    #################################################
    ###### methods for routing <IsoRing>
    #################################################

    def load_TDirectors_for_IRC(self,rdict,tdict):
        assert type(rdict) in {int,dict}
        assert type(tdict) in {float,dict}

        def fetch_rt(sec_idn):
            if type(rdict) == dict:
                x = rdict[sec_idn]
            else: 
                x = rdict

            if type(tdict) == dict:
                y = tdict[sec_idn]
            else: 
                y = tdict
            return x,y 

        for i in self.irc.irl:
            q = i.sec.idn_tag
            x,y = fetch_rt(q)
            ##print("RT: ",x,y)
            l = self.node_loc_assignment[q]
            i.default_TDirector_instance(l,x,y)

    def is_terminated(self):
        self.terminated = self.energy.v <= 0.0 
        return self.terminated 

    #################### open info mode
    def clear_open_info(self):
        for ir in self.irc.irl:
            ir.td.td.open_info_var.clear()

############################################################
############################################################

def SecNet_sample1(ss=None):
    if type(ss) != SecSeq:
        ss = SecSeq_sample_1(1)

    #ss = SecSeq_sample_1(1)
    sec_node_count = 12
    nsec_node_count = 23
    num_entry = 4
    rnd_struct = random
    sn = SecNet.generate(ss,sec_node_count,\
            nsec_node_count,num_entry,\
            10,rnd_struct,"spine frame",772) 
    return sn 

def pickled_SecNet_sample_Q(): 
    s = SecNet_sample1(SecSeq_sample_2(9,55)) 

    for s_ in s.irc.irl:
            s_.explode_contents()

    s.pickle_thyself("codename__ASS_SHIT")

"""
sample for testing <TDirector> in 1v1 setting. 
"""
def SecNet_sample_TDir1v1():
    random.seed(100)
    ss = SecSeq_sample_2(num_secs=1,num_conn=1,\
        min_components=1,drange_max=1)

    sec_node_count = 19
    nsec_node_count = 10
    num_entry = 4
    rnd_struct = random
    sn = SecNet.generate(ss,sec_node_count,\
            nsec_node_count,num_entry,\
            10,rnd_struct,"pairing frame",223) 

    for s_ in sn.irc.irl:
        s_.explode_contents()
    return sn 

"""
sample used to demonstrate correctness of 
<TDirector>. 
"""
def SecNet_sample_C3():

    random.seed(14324)
    np.random.seed(14324)
    ss,sndg = SecSeq_sample_4(num_secs=1,singleton_range=DEFAULT_SINGLETON_RANGE,\
        num_conn=5000,min_components=1,max_nconn_ratio=0.4,drange_max=1)

    G = SecNet_graph_sample_C3()
    sec_nodeset = {3,6,9,14}
    node_loc_assignment = {0:4}
    entry_points = {1,5,10,12,17} 
    bound = DEFAULT_SINGLETON_RANGE
    rnd_struct = random
    rnd_seed = 94312

    sn = SecNet(ss,G,sec_nodeset,\
            node_loc_assignment,entry_points,\
            bound,rnd_struct,rnd_seed)
    return sn

# TODO: needs to be used in tests. 
"""
sample for testing <TDirector> on Nv1.
"""
def SecNet_sample_TDirNv1():
    random.seed(100734)
    np.random.seed(371224)
    ss,sndg = SecSeq_sample_4(num_secs=25,singleton_range=DEFAULT_SINGLETON_RANGE,\
        num_conn=100,min_components=3,max_nconn_ratio=0.8,drange_max=4)

    sec_node_count = 25
    nsec_node_count = 34
    num_entry = 4
    rnd_struct = random
    sn = SecNet.generate(ss,sec_node_count,\
            nsec_node_count,num_entry,\
            1,rnd_struct,"pairing frame",223) 

    for s_ in sn.irc.irl:
        ##print("exploding {}".format(s_.sec.idn_tag))
        s_.explode_contents(num_blooms=2)
    return sn 

# TODO: needs to be used in tests. 
def SecNet_sample_CSmall():
    G = SecNet_graph_sample_CSmall()

    random.seed(5721)
    np.random.seed(935223)
    ss,sndg = SecSeq_sample_4(num_secs=2,singleton_range=DEFAULT_SINGLETON_RANGE,\
        num_conn=50,min_components=1,max_nconn_ratio=0.4,drange_max=1)

    sec_nodeset = set()  
    sn = SecNet(ss,G,sec_nodeset,{0:2,1:4},\
        entry_points = 1,bound=DEFAULT_SINGLETON_RANGE,\
        rnd_struct=random,rnd_seed=9,\
        path_size=10,sngc=None,energy=1000.0)
    return sn 

def SecNet_sample_CallSia():
    random.seed(25)
    np.random.seed(25)
    ss,sndg = SecSeq_sample_5(num_secs=3,max_nconn_ratio=0.0,\
    depconn_ratio=1.0,min_components=1,drange_max=3)
    sn = SecNet.generate(ss,6,10,3,4,random,"pairing frame",1015)
    return sn 

def SRefMap_sample1(): 
    sn = SecNet_sample1()
    return sn.srm 
    
###################

def SecNet_sample_approxhyp():

    random.seed(2004)
    np.random.seed(2004)

    ss,_ = SecSeq_sample_4(num_secs=1,\
            singleton_range=DEFAULT_SINGLETON_RANGE,\
            num_conn=1,min_components=1,max_nconn_ratio = 0.3,\
            drange_max=1)

    sn = SecNet.generate(ss,0,\
            1,1,10,random,"spine frame",223) 

    for s_ in sn.irc.irl:
        s_.explode_contents()

    return sn 
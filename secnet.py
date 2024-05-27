from collections import defaultdict 
from bridge import * 
from tdir import * 

class IsoRingedChain:

    def __init__(self,ss:SecSeq,bound,\
        rnd_struct,rs_seed:int):
        self.ss = ss 
        self.irl = []
        self.load_IsoRings(bound,rnd_struct,rs_seed)  

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
        ##print("RND_STRUCT: ",rnd_struct)

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
        ##print("RND STRUCT: ",rnd_struct)
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
            if q_ == "excessory": continue

            fp_ = fp + "/" + q_
            ir = IsoRing.unpickle_thyself(fp_) 
            ls.append(ir) 

        # sort all elements by idn
        ls = sorted(ls,key=lambda x: x.sec.idn_tag) 

        ss = SecSeq([])
        irc = IsoRingedChain(ss,None,rnd_struct,rs_seed)
        irc.irl = ls 
        return irc

# TODO: test this. 
"""
A graph structure that serves as an 
environment for activity programmed 
in other structs. 
"""
class SecNet:

    def __init__(self,irc,G,sec_nodeset,\
        node_loc_assignment= None,entry_points=3,\
        bound=DEFAULT_SINGLETON_RANGE,\
        rnd_struct=random,rnd_seed=9):
        
        assert len(irc) > 0 and type(irc) == SecSeq
        assert len(sec_nodeset) >= len(irc) 
        self.ss = irc
        ##print("RNDSTRUCT#1: ",rnd_struct)
        self.irc = IsoRingedChain(irc,bound,rnd_struct,rnd_seed) 
        self.rnd_struct = rnd_struct 

        self.d = G 
        self.sec_nodeset = sec_nodeset
        # sec index -> node location 
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
        self.occ_crackl = [] 
        self.preprocess_shortest_paths() 
        return

    @staticmethod
    def alt_instantiate(irc,G,sec_nodeset,\
        node_loc_assignment,srm,entry_points,\
        bound=DEFAULT_SINGLETON_RANGE,\
        rnd_struct=random,rnd_seed=9):

        assert type(irc) == IsoRingedChain
        assert type(srm) == SRefMap 

        ss = irc.revert_to_SecSeq()
        
        ##print("INSTANTIATE S.N.")
        sn = SecNet(ss,G,sec_nodeset,node_loc_assignment,\
            entry_points=1,bound=bound,\
            rnd_struct=rnd_struct,rnd_seed=rnd_seed)
        sn.irc = irc
        sn.entry_points = entry_points
        sn.srm = srm 
        return sn 
    ###
    """
SecNet:
    def __init__(self,irc,G,sec_nodeset,\
        node_loc_assignment= None,entry_points=3,\
        bound=DEFAULT_SINGLETON_RANGE,\
        rnd_struct=random,rnd_seed=9):
    """
    ###

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
        """
        self.d = G 
        self.sec_nodeset = sec_nodeset
        # sec index -> node location 
        self.node_loc_assignment = node_loc_assignment
        self.assign_entry(entry_points)
        if type(self.node_loc_assignment) == type(None):
            self.assign_loc()
        """
        excessory = [deepcopy(self.d),deepcopy(self.sec_nodeset),\
                deepcopy(self.node_loc_assignment),\
                deepcopy(self.srm),deepcopy(self.entry_points)]

        self.irc.pickle_thyself(fp)

        efile = open(fp + "/excessory","wb")
        pickle.dump(excessory,efile)
        efile.close()
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

        ##print("instantiating w/ rndstruct=",rnd_struct)
        sn = SecNet.alt_instantiate(irc,G,sec_nodeset,\
        node_loc_assignment,srefmap,entry_points,\
        bound=bound,\
        rnd_struct=rnd_struct,rnd_seed=rnd_seed)
        return sn 

    ######################## graph structure functions 
    # TODO: test 
    def subgraph_for_TDir(self,tdir):
        assert type(self.sgc) != type(None)
        rx = tdir.radius 
        return self.sgc.subgraph_by_radius_at_refnode(tdir.location,\
            rx,deepcopy(self.node_loc_assignment))

    def preprocess_shortest_paths(self):
        d2 = defaultdict(set)
        for k in self.d.keys():
            d2[k] = set(self.d[k])

        self.sgc = StdGraphContainer(d2,deepcopy(self.sec_nodeset))
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
            assert entry_points > 0 and entry_points < len(q)
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
        nsec_node_count,num_entry,rnd_struct,*sngs_args):
        q = sec_node_count + nsec_node_count
        assert q >= len(irc) 
        assert num_entry <= q and num_entry > 0 

        q_ = [i for i in range(q)] 
        rnd_struct.shuffle(q_)

        snv = q_[:sec_node_count]
        nsnv = q_[sec_node_count:]
        sngs = SecNetGenScheme(*sngs_args)
        snfg = SecNetFrameGen(snv,nsnv,sngs)
        snfg.construct_frame()

        node_loc_assignment = None
        sn = SecNet(irc,snfg.d,set(snfg.sec_nodevec),\
            node_loc_assignment,entry_points=num_entry,\
            rnd_struct=rnd_struct)
        return sn 

def SecNet_sample1(ss=SecSeq_sample_1(1)):
    #ss = SecSeq_sample_1(1)
    sec_node_count = 12
    nsec_node_count = 23
    num_entry = 4
    rnd_struct = random
    sn = SecNet.generate(ss,sec_node_count,\
            nsec_node_count,num_entry,\
            rnd_struct,"spine frame",772) 
    return sn 

def pickled_SecNet_sample_Q(): 
    s = SecNet_sample1(SecSeq_sample_2(9,55)) 

    for s_ in s.irc.irl:
            s_.explode_contents()
    """
    irc = s.irc
    irc.pickle_thyself("codename__ASS_SHIT")
    """
    s.pickle_thyself("codename__ASS_SHIT")


def SRefMap_sample1(): 
    sn = SecNet_sample1()
    return sn.srm 
    
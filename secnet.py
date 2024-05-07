from collections import defaultdict 
from bridge import * 

class IsoRingedChain:

    def __init__(self,ss:SecSeq,bound,\
        rnd_struct,rs_seed:int):
        self.ss = ss 
        self.irl = []
        self.load_IsoRings(bound,rnd_struct,rs_seed)  

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

    def load_IsoRings(self,singleton_bound=DEFAULT_SINGLETON_RANGE,\
        rnd_struct=random,rs_seed=8): 

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
    
    def pickle_thyself(self):
        return -1

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
        return

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

def SecNet_sample1():
    ss = SecSeq_sample_1(1)
    sec_node_count = 12
    nsec_node_count = 23
    num_entry = 4
    rnd_struct = random
    sn = SecNet.generate(ss,sec_node_count,\
            nsec_node_count,num_entry,\
            rnd_struct,"spine frame",772) 
    return sn 

def SRefMap_sample1(): 
    sn = SecNet_sample1()
    return sn.srm 
    
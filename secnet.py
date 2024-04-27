from collections import defaultdict 
from bridge import * 

"""
A graph structure that serves as an 
environment for activity programmed 
in other structs. 
"""
class SecNet:

    def __init__(self,irc,G,sec_nodeset,\
        node_loc_assignment= None,
        entry_points=3,
        rnd_struct=random):
        assert len(irc) > 0 and type(irc) == SecSeq
        
        #for i in irc: assert type(i) == IsoRing
        assert len(sec_nodeset) >= len(irc) 

        self.irc = irc 
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
        return

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
    a sequence of <IsoRings> (that may or may not 
    have dep./codep.) along with other variables 
    pertaining to the graph to be generated.
    """
    @staticmethod
    def generate(irc,sec_node_count,\
        nsec_node_count,num_entry,rnd_struct,*sngs_args):
        assert node_count >= len(irc) 
        q = node_count + nsec_node_count
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

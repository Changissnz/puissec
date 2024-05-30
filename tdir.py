from dfs_struct import *

"""
a std. graph is of type dict|defaultdict
"""
class StdGraphContainer:

    def __init__(self,m,sec_nodeset):
        assert type(m) in {dict,defaultdict}
        self.d = m
        self.sn = sec_nodeset 
        self.dfs_cache_map = defaultdict(None)
        self.sp = defaultdict(NodePath)
        self.ring_locs = {}

    def DFSCache_proc(self,n):
        dfsc = DFSCache(n,deepcopy(self.d),\
                search_head_type=1)
        dfsc.exec_DFS()
        dfsc.store_minpaths(num_paths=1)
        self.sp[n] = dfsc
        return 

    def DFSCache_fullproc(self):
        for k in self.d.keys():
            self.DFSCache_proc(k)

    def subgraph_by_radius_at_refnode(self,r,radius,\
        node_loc_assignment):
        ns = {r}

        dfsc_ = self.sp[r] 

        print("SP")

        # fetch the `subgraph`
        for k in dfsc_.min_paths.keys():
            # fetch the path 
            s = dfsc_.min_paths[k]
            if len(s) == 0: continue
            s = np.sum(s[0].pweights)
            if s <= radius:
                ns = ns | {k}

        nla = {} 
        if type(node_loc_assignment) != type(None): 
            for k,v in node_loc_assignment.items():
                if v in ns:
                    nla[k] = v

        snx = ns.intersection(self.sn)
        sg = self.subgraph(ns)

        spx = defaultdict(NodePath)
        for ns_ in ns:
            q = deepcopy(self.sp[ns_])
            k = set(q.min_paths.keys())

            for k_ in k:
                if k_ not in ns:
                    del q.min_paths[k_] 
            spx[ns_] = q 

        sgc = StdGraphContainer(sg,snx)
        sgc.ring_locs = nla 
        sgc.sp = spx
        return sgc 

    def subgraph(self,ns): 
        assert type(ns) == set
        dx = {} 
        for n in ns:
            dx[n] = ns.intersection(self.d[n])
        return defaultdict(set,dx)


"""
Was supposed to be named Traversal Directing [Wang Fong Qhong].
By and through the connection. 
By,through, and for the connection.
"""

class TDir:

    def __init__(self,loc,target_node,\
        vantage_point,radius=4,velocity=1):
        assert vantage_point in {"I","C"}
        self.location = loc
        self.target_node = target_node
        self.vantage_point = vantage_point
        self.radius = radius
        self.velocity = velocity 
        self.node_path = None
        self.index = None
        # time reference to start before travel. 
        self.t = 0.0
        self.t_ = 0.0
        return

    def load_path(self,G):
        assert type(G) == StdGraphContainer
        assert self.location in G.d

        target_loc = self.search_for_target_node(self.target_node)

        if type(target_loc) == type(None):
            print("error: target node not found")
            return

        nodePath = G.d[self.location][target_loc][0]
        self.node_path = nodePath.invert() 
        self.index = 0

    def search_for_target_node(self,G):
        assert type(G) == StdGraphContainer

        if self.target_node not in self.ring_locs:
            return None

        loc = self.ring_locs[self.target_node]
        return loc 

    # TODO: test 
    def __next__(self):
        if self.index >= len(self.node_path):
            return None
        self.index += self.velocity
        q = min([len(self.node_path) - 1,self.index])

        q = self.node_path.p[self.index]
        self.location = q
        return q

    def scaled__next__(self,scale = 1.0):

        tx = self.t + scale
        r = int(tx - self.t_)

        stat = True
        for i in range(r):
            q = self.__next__()
            if type(q) == type(None): 
                stat = not stat
                break 

        self.t = tx 
        self.t_ += r 
        return self.location,stat 

    """
    on results. 
    """
    def reflect(self):

        return -1




"""

"""
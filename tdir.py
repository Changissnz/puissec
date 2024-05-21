
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

        # fetch the `subgraph`
        for k in self.min_paths.keys():
            # fetch the path 
            s = self.min_paths[k]
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
        sg = self.subgraph(ns,snx)

        spx = defaultdict(NodePath)
        for ns_ in ns:
            spx[ns_] = self.sp[ns_]

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
        vantage_point,velocity=1):
        assert vantage_point in {"I","C"}
        self.location = loc
        self.target_node = target_node
        self.vantage_point = vantage_point
        self.velocity = velocity 
        self.node_path = None
        self.index = None
        return

    def load_path(self,G):
        assert type(G) == StdGraphContainer
        assert self.location in G.d
        nodePath = G.d[self.location][self.target_node][0]
        self.node_path = nodePath.invert() 
        self.index = 0

    def __next__(self):
        if self.index >= len(self.node_path):
            return None
        self.index += self.velocity
        q = min([len(self.node_path) - 1,self.index])

        q = self.node_path.p[self.index]
        return q
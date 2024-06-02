from dfs_struct import *

"""
a std. graph is of type dict|defaultdict
"""
class SNGraphContainer:

    def __init__(self,m,sec_nodeset,ring_locs,clocs):
        assert type(m) in {dict,defaultdict}
        self.d = m
        self.sn = sec_nodeset 
        ##self.dfs_cache_map = defaultdict(None)
        self.sp = defaultdict(NodePath)
        # <IsoRing> idn -> location 
        self.ring_locs = ring_locs
        # <Crackling> idn -> (location,target node)
        self.crackling_locs = clocs

    def agent_loc(self,agent_idn,is_isoring:bool):
        if is_isoring:
            if agent_idn not in self.ring_locs: return None
            return self.ring_locs[agent_idn]

        if agent_idn not in self.crackling_locs:
            return None
        return self.crackling_locs[agent_idn][0]

    def agents_at_loc(self,loc,is_isoring:bool): 
        ks = []
        if is_isoring: 

            for k,v in self.ring_locs.items():
                if v == loc: 
                    ks.append(k) 
            return ks

        for k,v in self.crackling_locs.items():
            if v[0] == loc: 
                ks.append(k)
        return ks


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

    def subgraph_by_radius_at_refnode(self,r,radius):
        ns = {r}

        dfsc_ = self.sp[r] 

        # fetch the `subgraph`
        for k in dfsc_.min_paths.keys():
            # fetch the path 
            s = dfsc_.min_paths[k]
            if len(s) == 0: continue
            s = np.sum(s[0].pweights)
            if s <= radius:
                ns = ns | {k}

        nla = {} 
        if type(self.ring_locs) != type(None): 
            for k,v in self.ring_locs.items():
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

        occm = {}
        for k,v in self.crackling_locs.items():
            if v[0] in ns: 
                occm[k] = v

        sgc = SNGraphContainer(sg,snx,nla,occm)
        sgc.sp = spx
        return sgc 

    def subgraph(self,ns): 
        assert type(ns) == set
        dx = {} 
        for n in ns:
            dx[n] = ns.intersection(self.d[n])
        return defaultdict(set,dx)

    def update_rc_agent_locs(self,rlocs,clocs):
        self.ring_locs = rlocs
        self.crackling_locs = clocs 
        return


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

        self.active_stat = True 
        return

    def load_path(self,G):
        assert type(G) == SNGraphContainer
        assert self.location in G.d
        ##print("RING LOCS: ",G.ring_locs)
        ##print("...")

        target_loc = self.search_for_target_node(G)

        if type(target_loc) == type(None):
            print("error: target node not found")
            return

        ##print("LOC->TARGET")
        ##print(G.d)
        ##print("##")
        dfsc = G.sp[self.location]
        nodePath = dfsc.min_paths[target_loc][0]
        self.node_path = nodePath.invert() 
        self.index = 0

    def search_for_target_node(self,G):
        assert type(G) == SNGraphContainer
        ##print("TARGET NODE: ",self.target_node)

        # case: vp = I; literal node destination 
        if self.vantage_point == "I": 
            return self.target_node

        # case: vp = C; target <IsoRing> idn. 
        if self.target_node not in G.ring_locs:
            return None

        loc = G.ring_locs[self.target_node]
        return loc 

    # TODO: test 
    def __next__(self):
        if not self.active_stat: return None 

        if self.index >= len(self.node_path):
            return None
        self.index += self.velocity
        q = min([len(self.node_path) - 1,self.index])
        self.index = q 
        q = self.node_path.p[self.index]
        self.location = q

        if self.index == len(self.node_path) - 1:
            self.index += 1
            self.active_stat = False 

        return q

    def scaled__next__(self,scale = 1.0):
        if not self.active_stat: 
            return self.location,False
            
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

TDIR_SETTINGS = {'I':{"avoid target",\
                    "radar null"},\
                'C': {"search for target",\
                    "capture target"}}

"""
<TDirector> directs an associated agent <IsoRing>
or <Crackling> based on objectives.

<TDirector> requires a feedback loop with the 
connected <SecNet> it operates on.

The <SecNet> handles the subgraph of information 
given to the <TDirector> instance for any process 
step that demands that information. 
"""
class TDirector:

    def __init__(self,loc,target_node,\
        vantage_point,vantage_idn,radius=4,velocity=1):
        self.td = TDir(loc,target_node,vantage_point,\
            radius,velocity) 
        self.vantage_idn = vantage_idn

        self.ref_nodes = [deepcopy(self.td.location)]

        self.td_log = [] 

        # Resource subgraph, of type <SNGraphContainer>
        self.resource_sg = None

        self.obj_stat = None
        if self.td.vantage_point == "I":  
            self.obj_stat = "radar null"
        else:
            self.obj_stat = "search for target"
        return

    def vp(self):
        return self.td.vantage_point 

    def load_graph(self,G:SNGraphContainer):
        assert type(G) == SNGraphContainer
        self.resource_sg = G 

    def loc(self):
        return self.td.location

    def update_tdir(self):
        return -1

    def check_obj(self):


        return -1

    """
    checks for chaser <Crackling> if 
    vantage_point = I, otherwise checks 
    for target <IsoRing>.
    """
    def check_radar(self):
        
        assert type(self.resource_sg) != type(None) 
        
        # check for any Cracklings that are 
        # targetting
        if self.vp() == "I":
            q = []
            for k,v in self.resource_sg.crackling_locs.items():
                stat = v[1] == self.vantage_idn
                if stat:
                    q.append(k)
            return q 
        
        # check for target isoring loc
        else: 
            if self.td.target_node in self.sg.ring_locs:
                return [self.sg.ring_locs[self.td.target_node]]
            return []

    #######################
    ######### loc search mode: used by <Crackling> to locate
    #########                  target <IsoRing>.

    """
    searches for a target node to travel the agent I|C. 
    """
    def loc_search__set_TDir(self,sgc:SNGraphContainer,\
        rnd_struct=random):#,check_completion=True):
        q = self.loc_search_at_ref(sgc,rnd_struct)
        if type(q) == type(None):
            return 

        l = self.td.location
        tn = self.td.target_node
        vp = self.td.vantage_point
        r = self.td.radius
        v = self.td.velocity

        tn = q 
        tdx = TDir(l,tn,vp,r,v)
        self.td_log.append(self.td)
        self.td = tdx 
        return

    """
    searches for the next location. 

    """
    def loc_search_at_ref(self,sgc:SNGraphContainer,rnd_struct,\
        exclude_refs:bool=True):
        q = self.ref_nodes[-1] 
        candidate_list = self.loc_search_set(sgc,q)

        if exclude_refs:
            candidate_list = candidate_list - set(self.ref_nodes)

        if len(candidate_list) == 0:
            return None 

        candidate_list = sorted(list(candidate_list))
        xi = rnd_struct.randrange(0,len(candidate_list))
        cl = candidate_list[xi]
        return cl 

    """
    return:
    - set, nodes of the greatest distance d
    """
    def loc_search_set(self,sgc:SNGraphContainer,ref):
        if len(sgc.sp.keys()) == 0:
            return set()
        
        rkx = sorted([(k,sum(v.pweights)) for k,v in sgc.sp.items()],\
            key=lambda x:x[1])[::-1]
        
        mx = rkx[0]
        st = set()
        st = st | {mx[0]}
        for i in range(1,len(rkx)):
            if rkx[i][1] == mx[1]: 
                st = st | {rkx[i][0]}
            else:
                break  
        return st 

# TODO: tests
'''
loc search: used for objective C.search_for_target
'''
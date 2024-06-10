from dfs_struct import *

"""
an <SNGraphContainer> (SN := SecNet) holds data of 
type dict|defaultdict. 
"""
class SNGraphContainer:

    def __init__(self,m,sec_nodeset,ring_locs,clocs,\
        entry_points,path_size=10):
        assert type(m) in {dict,defaultdict}
        self.d = m
        self.sn = sec_nodeset 
        ## node idn -> corresponding <DFSCache> instance 
        self.sp = defaultdict(DFSCache)
        # <IsoRing> idn -> location 
        self.ring_locs = ring_locs
        # <Crackling> idn -> (location,target node)
        self.crackling_locs = clocs
        # set, node idn. for entry points
        self.entry_points = entry_points
        # 

        self.path_size = path_size

    def nsec_nodeset(self):
        return set(self.d.keys()) - self.sn

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
        dfsc.store_minpaths(num_paths=self.path_size)
        self.sp[n] = dfsc
        return 

    def DFSCache_fullproc(self):
        for k in self.d.keys():
            ##print('dfs proc for {}'.format(k))
            self.DFSCache_proc(k)

    def subgraph_by_radius_at_refnode(self,r,radius):
        ns = {r}

        dfsc_ = self.sp[r] 

        # fetch the `subgraph`
        for k in dfsc_.min_paths.keys():
            # fetch the path 
            s = dfsc_.min_paths[k]
            if len(s) == 0: continue
            s = s[0].cost()
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

        entry_points = self.entry_points.intersection(ns)

        sgc = SNGraphContainer(sg,snx,nla,occm,entry_points)
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

    def __str__(self):
        s = "<TDIR>\n"
        s += "\t* loc: {}\n".format(self.location)
        s += "\t* target: {}\n".format(self.target_node)
        s += "\t* vantage: {}\n".format(self.vantage_point)
        s += "\t* radius: {}\n".format(self.radius)
        s += "\t* velocity: {}\n".format(self.velocity)
        s += "\t* nodepath:\n{}".format(str(self.node_path))
        return s 


    def load_path(self,G):
        assert type(G) == SNGraphContainer
        assert self.location in G.d
        ##print("RING LOCS: ",G.ring_locs)
        ##print("...")
        ##print("LOADING PATH")

        target_loc = self.search_for_target_node(G)

        if type(target_loc) == type(None):
            self.active_stat = False 
            print("error: target node not found")
            return

        p = self.load_path_(G,target_loc)

        if type(p) != NodePath:
            print("no path")
            return
        self.node_path = p
        self.index = 0

        ## TODO: delete? 
        ##print("LOC->TARGET")
        ##print(G.d)
        ##print("##")
        """
        dfsc = G.sp[self.location]
        nodePath = dfsc.min_paths[target_loc][0]
        self.node_path = nodePath.invert() 
        self.index = 0
        """ 
    
    def load_path_(self,G,loc):
        if loc not in G.sp:
            print("no location in DFS")
            return

        dfsc = G.sp[loc]
        if self.location not in dfsc.min_paths:
            print("no path in DFS")
            return

        nodePath = dfsc.min_paths[self.location][0]
        return nodePath


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
        # log of <TDir> instances
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
        self.td.load_path(G) 

    def loc(self):
        return self.td.location

    def switch_obj_stat(self):
        if self.obj_stat == "avoid target":
            self.obj_stat = "radar null"
        elif self.obj_stat == "radar null":
            self.obj_stat = "avoid target"
        elif self.obj_stat == "search for target":
            self.obj_stat = "capture target"
        else:
            self.obj_stat = "search for target"

    def update_tdir(self):
        return -1

    # TODO: 
    """
    return: 
    - bool, if `obj_stat` in {radar null,search for target},
      otherwise <NodePath> instance to update the <TDir>
    """
    def check_obj(self):
        assert self.vp() in {"I","C"}

        if self.obj_stat in {"radar null",\
                "search for target"}:
            q = self.check_radar()
            return len(q) > 0

        if self.obj_stat == "avoid target":
            return -1

        if self.obj_stat == "capture target": 
            return -1
        return

    """
    checks for chaser <Crackling> if 
    vantage_point = I, otherwise checks 
    for target <IsoRing>.

    return:
    - set, node locations
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
            return set(q) 
        
        # check for target isoring loc
        else: 
            if self.td.target_node in self.resource_sg.ring_locs:
                return set([self.resource_sg.ring_locs[\
                    self.td.target_node]])
            return set()

    #######################
    ######### loc search mode: used by agent A to achieve its objective.
    #########               - <Crackling> to locate target <IsoRing>.
    #########               - <IsoRing> to distance itself from pursuing <Crackling> instances. 


    ########################## section: extremum locator
    """
    searches for a target node to travel the agent I|C. 
    """
    def extloc_search__set_TDir(self,extf=max,rnd_struct=random): 
    #(self,sgc:SNGraphContainer,extf=max,rnd_struct=random):#,check_completion=True):
        print("INITIATE EXTLOC")

        assert type(self.resource_sg) == SNGraphContainer
        q = self.extloc_search_at_ref(self.resource_sg,rnd_struct,extf)
        if type(q) == type(None):
            return 

        print("node EXTLOC: ",q)

        l = self.td.location
        tn = self.td.target_node
        vp = self.td.vantage_point
        r = self.td.radius
        v = self.td.velocity

        ##tn = q 
        print("target <Sec>: {} node dest: {}".format(tn,q))

        tdx = TDir(l,tn,vp,r,v)
        print("loading EXTLOC path")
        node_path = tdx.load_path_(self.resource_sg,q) 
        if type(node_path) != NodePath: 
            print("no node path")
            return
        tdx.node_path = node_path
        tdx.index = 0

        self.td_log.append(self.td)
        self.td = tdx 
        return

    """
    searches for the next location. 

    return: 
    - node, of extreme distance d selected by `rnd_struct`.
    """
    def extloc_search_at_ref(self,sgc:SNGraphContainer,rnd_struct,\
        exclude_refs:bool=True,extf=max):
        q = self.ref_nodes[-1] 
        candidate_list = self.extloc_search_set(sgc,q,extf)

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
    - set, nodes of the most extreme distance d (depending on `extf`). 
    """
    def extloc_search_set(self,sgc:SNGraphContainer,ref,extf):
        assert extf in {min,max}

        if len(sgc.sp.keys()) == 0:
            return set()
        
        ##self.sp = defaultdict(DFSCache)
        rkx = []
        for k,v in sgc.sp.items():
            if ref not in v.min_paths: 
                continue

            q = v.min_paths[ref][0].cost()
            rkx.append((k,q))

        #rkx = sorted([(k,sum(v.pweights)) for k,v in sgc.sp.items()],\
        #    key=lambda x:x[1])
        if extf == max:
            rkx = rkx[::-1]

        mx = rkx[0]
        st = set()
        st = st | {mx[0]}
        for i in range(1,len(rkx)):
            if rkx[i][1] == mx[1]: 
                st = st | {rkx[i][0]}
            else:
                break  
        return st

    ###################### section: SEC-node locator

    """
    produces metrics for the choice of 
    target node `tn` as a destination,
    given the vantage point of the 
    <TDirector>. 
    """
    def targetnode_analysis(self,tn,rnd_struct):

        if self.vp() == "C": 
            return self.targetnode_analysis__VPC(tn,rnd_struct)
        return self.targetnode_analysis__VPI(tn,rnd_struct)

    """
    analysis of node `tn` as a destination node
    given vantage point <IsoRing>.

    return:
    - dict, <Crackling> idn -> 
        (distance to tn / mean distance to SEC nodes).
    """
    def targetnode_analysis__VPI(self,tn):

        # cracklings detected by radar
        sn = self.check_radar()

        if len(sn) == 0: 
            return None

        # get the distance of each crackling 
        # to `tn`.
        d = {}
        dfsc = self.resource_sg.sp[tn]#self.td.location]
        for sn_ in sn:
            if sn_ not in dfsc.min_paths:
                continue
            d[sn_] = dfsc.min_paths[sn_][0].cost()

        # get <Crackling> distance to SEC node
        # crackling -> mean distance to SEC node
        
        ## MAYBE: min? 
        mean_distances = self.mean_mindistance_of_nodeset(deepcopy(sn))
        d2 = {}
        for k,v in d.items():
            mds = mean_distances[k] 
            d2[k] = v / mds  
        return d2

    def mean_mindistance_of_nodeset(self,ns):
        d = {}
        for n in ns:
            dx = self.secnode_distance_map(n)
            sx = measures.zero_div(sum(dx.values()),len(dx),0.0)
            d[n] = sx
        return d 

    # TODO: needs to be demonstrated. 
    """
    """
    def targetnode_analysis__VPC(self,tn,rnd_struct):

        if tn not in self.resource_sg.sp:
            return None
        dfsc = self.resource_sg.sp[tn]

        cs = self.check_radar()

        if len(cs) == 0: 
            return None


        # find the shortest path to cs
        assert len(cs) == 1

        # spot the security of I's location
        cs = cs.pop()
        stat = cs in self.resource_sg.sn 

        # decision: if location is SEC, then 
        #           travel to it. Otherwise,
        #           select an NSEC for site of 
        #           interception.
        if stat:

            if cs not in dfsc.min_paths:
                return None
            return dfsc.min_paths[cs][0] 
        else:
            sndict = self.secnode_distance_map(tn,is_sec=False)
            if len(sndict) == 0: return None

            sndict_ = [(k,v) for k,v in sndict.items()] 
            sndict_ = sorted(sndict_,key=lambda x: x[1])

            qs = sndict_.pop(0)
            sx = set([qs[0]])
            sc = qs[1]
            while len(sndict_) > 0:
                snd = sndict_.pop(0)
                if snd[1] != sc:
                    break
                else: 
                    sx = sx | {snd[0]}

            # choose a random nsec node of min distance to
            # tn
            sx = list(sx)
            sxi = rnd_struct.randrange(0,len(sx))
            target = sx[sxi]

            if target not in dfsc.min_paths:
                return None
            return dfsc.min_paths[target][0]

    # TODO: future
    def load_predicted_travel_at_node(self,n,npath):
        return -1 

    """
    calculates distance of node `n` to each 
    of the sec nodes. 

    return:
    - dict, sec node idn -> distance
    """
    def secnode_distance_map(self,n,is_sec:bool=True):

        if n not in self.resource_sg:
            return {}
        
        sns = deepcopy(self.resource_sg.sn)
        if not is_sec:
            sns = self.resource_sg.nsec_nodeset()

        d = {} 
        dfsc = self.resource_sg.sp[n]
        for s in sns:
            if s not in dfsc.min_paths: 
                continue
            npath = dfsc.min_paths[s][0]
            d[s] = npath.cost()
        return d
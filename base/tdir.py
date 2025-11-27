from .dfs_struct import *
from morebs2 import measures

DEFAULT_TDIRECTOR_OBJ_FUNCTIONS = [np.min,\
                    np.max,np.mean,std_iscaled_mean]

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

    def display(self,verbose=1):
        assert verbose in {0,1}

        if verbose: 
            print("DICTIONARY")
            for k,v in self.d.items():
                print("* {} --> {}".format(k,v))
            print()

        print("SEC NODESET")
        print(self.sn)
        print()
        print("ISORING")
        print(self.ring_locs)
        print()
        print("CRACKLING")
        print(self.crackling_locs)
        print()
        print("ENTRY POINTS")
        print(self.entry_points)


    @staticmethod
    def unpickle_thyself(fp): 
        fx = open(fp,"rb")

        fobj = pickle.load(fx) 
        assert type(fobj) == SNGraphContainer
        
        for v in fobj.sp.values():
            v.ecf = DEFAULT_EDGE_COST_FUNCTION
        return fobj 

    def pickle_thyself(self,fp):
        fx = open(fp,"wb")

        qr = {}

        for k,v in self.sp.items():
            qr[k] = v.ecf
            v.ecf = None

        pickle.dump(self,fx)
        fx.close()

        for k,v in qr.items():
            self.sp[k].ecf = v 

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
    return:
    - set(other agents at same location as `aidn`)
    """
    def coloc_for_agent(self,aidn,is_crackling:bool):
        l,q = None,None
        if is_crackling:
            if aidn not in self.crackling_locs: assert False
            l = self.crackling_locs[aidn][0]
            q = self.ring_locs 

        else: 
            if aidn not in self.ring_locs: assert False
            l = self.ring_locs[aidn]
            q = self.crackling_locs
            
        sx = set()
        print("[BUG] ALL CO-LOCS FOR {}:{}".format(aidn,is_crackling))
        print(q)

        for k,v in q.items():
            v2 = v if is_crackling else v[0] 
            if l == v2: sx = sx | {k}

        return sx 

    #####################################
    ############# methods for calculating boundary points
    ############# of simple, undirected graphs

    def eccentricity_map(self):
        dx = {}

        for k,v in self.sp.items():
            qxs = [v2[0].cost() for v2 in v.min_paths.values()]
            dx[k] = max(qxs)
        return dx 

    """
    return: 
    - list(set::(peripheral nodeset of i'th order)) 
    """
    def peripheral_nodes_by_order(self,orderionos=1):
        assert type(orderionos) == int 
        assert orderionos >= 1 

        em = self.eccentricity_map()
        if len(em) == 0: return None
        dx = []
        sv = sorted(list(em.values()))[::-1]
        sv2 = [sv.pop(0)]

        i = 1
        stat = True 
        candidates = list(em.keys())
        while stat:
            ixs = [] 
            ns = [] 
            for (i2,c) in enumerate(candidates):
                v_ = em[c]

                if round(abs(v_ - sv2[-1]),5) == 0.0:
                    ixs.append(i2)
                    ns.append(c)

            ns2 = [c for (i2,c) in enumerate(candidates) if \
                i2 not in ixs]
            candidates = ns2
            dx.append(set(ns)) 

            sv2.append(sv.pop(0))
            i += 1

            if len(candidates) == 0:
                stat = False
                continue 

            if len(sv) == 0: 
                stat = False 
                continue

            if i > orderionos: 
                stat = False 

        return dx 

    def radius(self):
        return self.graph_distance_ext(min)

    def diameter(self):
        return self.graph_distance_ext(max) 

    def graph_distance_ext(self,extf):
        assert extf in {min,max}

        em = self.eccentricity_map()
        q = list(em.values())

        if len(q) == 0: return -1
        return extf(q) 

    #####################################

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

        self.active_stat = True # False 
        self.open_info_var = [] 
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
        
        ##print("TARGET NODE FOUND! ",target_loc)

        p = self.load_path_(G,target_loc)

        if type(p) != NodePath:
            print("no path")
            return
        self.node_path = p
        self.index = 0
    
    def load_path_(self,G,loc):
        if loc not in G.sp:
            print("no location in DFS")
            return

        dfsc = G.sp[loc]
        if self.location not in dfsc.min_paths:
            print("no path in DFS")
            return
        ##print("YESLOC, loc {} target {}".format(self.location,loc))

        nodePath = dfsc.min_paths[self.location][0]
        return nodePath
        
    def search_for_target_node(self,G):
        assert type(G) == SNGraphContainer
        ##print("LOCA: ",self.location)
        ##print("TARGET NODE: ",self.target_node)
        ##print("RING LOCS: ",G.ring_locs)

        # case: vp = I; literal node destination 
        if self.vantage_point == "I": 
            return self.target_node

        # case: vp = C; target <IsoRing> idn. 
        if self.target_node not in G.ring_locs:
            return None

        loc = G.ring_locs[self.target_node]
        return loc 

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
            print("NOT ACTIVE")
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
        vantage_point,vantage_idn,radius=4,\
        velocity=0,tdts=DEFAULT_TDIRECTOR_TIMESTAMP_SIZE):
        self.td = TDir(loc,target_node,vantage_point,\
            radius,velocity) 
        self.tdts = tdts
        self.vantage_idn = vantage_idn

        self.ref_nodes = None 
        self.reset_ref_nodes()
        # log of <TDir> instances
        self.td_log = [] 

        # Resource subgraph, of type <SNGraphContainer>
        self.resource_sg = None

        self.obj_stat = None
        if self.td.vantage_point == "I":  
            self.obj_stat = "radar null"
        else:
            self.obj_stat = "search for target"

        # performance scores
        self.ps = np.array([float('inf')])
        return

    def vp(self):
        return self.td.vantage_point 

    def load_graph(self,G:SNGraphContainer):
        assert type(G) == SNGraphContainer
        self.resource_sg = G
        self.td.load_path(G) 

    def load_new_path(self,npath:NodePath):
        assert type(npath) == NodePath 
        self.td.node_path = npath
        self.td.index = 0 
        self.td.active_stat = True 

    def loc(self):
        return self.td.location

    def is_loc_SEC(self):
        return self.is_node_SEC(self.loc())

    def is_node_SEC(self,n):
        return n in self.resource_sg.sn 

    def coloc(self):
        if type(self.resource_sg) != SNGraphContainer:
            return None

        is_crackling = True if self.vp() == "C" else False
        return self.resource_sg.coloc_for_agent(\
            self.vantage_idn,is_crackling)

    def switch_obj_stat(self):
        if self.obj_stat == "avoid target":
            self.obj_stat = "radar null"
        elif self.obj_stat == "radar null":
            self.obj_stat = "avoid target"
        elif self.obj_stat == "search for target":
            self.obj_stat = "capture target"
        else:
            self.obj_stat = "search for target"
        self.clear_data() 

    """
    clears the following variables: 
    - ps (the performance scores by reflection)
    - ref_nodes (collected during `search for target` by <Crackling>)
    """
    def clear_data(self):
        if self.obj_stat in {"avoid target","capture target"}:
            ##print("cannot delete data during avoid/capture ops.")
            return

        self.ps = np.array([])
        self.reset_ref_nodes()
        self.td.active_stat = False
        return

    def update_tdir(self):
        return -1

    # TODO: use this method to aid in agent decisions. 
    """
    return: 
    - bool, ?switch objective?
    """
    def check_obj(self):
        q = self.check_radar()
        if self.obj_stat in {"radar null",\
                "search for target"}:
            return len(q) > 0

        return len(q) == 0

    """
    checks for chaser <Crackling> if 
    vantage_point = I, otherwise checks 
    for target <IsoRing>.

    return:
    - set, node locations
    """
    def check_radar(self,is_agent_idn:bool=False):
        
        assert type(self.resource_sg) != type(None) 
        
        # check for any Cracklings that are 
        # targetting
        if self.vp() == "I":
            q = []
            for k,v in self.resource_sg.crackling_locs.items():
                stat = v[1] == self.vantage_idn
                if stat:
                    v_ = v[0] if not is_agent_idn else k
                    q.append(v_)
            return set(q) 
        
        # check for target isoring loc
        else: 
            if self.td.target_node in self.resource_sg.ring_locs:
                v_ = self.td.target_node if is_agent_idn else \
                    self.resource_sg.ring_locs[self.td.target_node]
                return set([v_])
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

        assert type(self.resource_sg) == SNGraphContainer
        q = self.extloc_search_at_ref(self.resource_sg,rnd_struct,extf)
        if type(q) == type(None):
            return 

        l = self.td.location
        tn = self.td.target_node
        vp = self.td.vantage_point
        r = self.td.radius
        v = self.td.velocity

        tdx = TDir(l,tn,vp,r,v)
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
        
        rkx = []
        for k,v in sgc.sp.items():
            if ref not in v.min_paths: 
                continue

            q = v.min_paths[ref][0].cost()
            rkx.append((k,q))

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

    def reset_ref_nodes(self):
        self.ref_nodes = [deepcopy(self.td.location)]


    ###################### section: SEC-node locator

    def default_node_analysis(self):
        ta = self.targetnode_analysis(np.min)
        td = self.destnode_sec_analysis_dict(np.mean)
        
        ##
        """
        print("TA")
        print(ta)
        print("TD")
        print(td)
        """
        ## 

        k = set(list(ta.keys()) + list(td.keys()))
        d = {}
        for k_ in k:
            d[k_] = 0 
            if k_ in ta: d[k_] = d[k_] + ta[k_]
            if k_ in td: d[k_] = d[k_] + td[k_]
        
        return d

    # TODO: test 
    """
    return: 
    - dict, node -> `targetnode_analysis_(node,objf)`
    """
    def targetnode_analysis(self,objf):
        qn = set(self.resource_sg.d.keys())
        d = {}
        for n in qn:
            score = self.targetnode_analysis_(n,objf)
            d[n] = score
        return d

    """
    return:
    - float, inf if no targets in sight, otherwise
            `objf(distance vector)`. 
    """
    def targetnode_analysis_(self,n,objf):
        assert objf in DEFAULT_TDIRECTOR_OBJ_FUNCTIONS

        sn = self.check_radar()

        if len(sn) == 0: return float('inf') 

        # TODO: refactor 
        d = {}
        dfsc = self.resource_sg.sp[n]
        for sn_ in sn:
            if sn_ not in dfsc.min_paths:
                continue
            d[sn_] = dfsc.min_paths[sn_][0].cost()

        vx = np.array(list(d.values()))

        return objf(vx) 

    """
    return: 
    - dict, node -> objf1(destnode_sec_analysis(node).keys)
    """
    def destnode_sec_analysis_dict(self,objf1): 
        assert objf1 in {np.min,np.mean}
        d = {}
        for n in self.resource_sg.d.keys():
            dx = self.destnode_sec_analysis(n)
            dxv = np.array(list(dx.values())) 
            assert len(dxv) > 0

            v = objf1(dxv)

            # case: v == inf
            if np.isinf(v): v = 0.0 

            d[n] = v 
        return d

    # TODO: test 
    """
    analysis of node `tn` as a destination node
    given vantage point <IsoRing>.

    return:
    - dict, location of opponent -> 
        (distance to tn / mean distance to SEC nodes).
    """
    def destnode_sec_analysis(self,tn):

        # cracklings detected by radar
        sn = self.check_radar()

        if len(sn) == 0:
            return None

        # get the distance of each crackling 
        # to `tn`.
        d = {}
        dfsc = self.resource_sg.sp[tn]
        for sn_ in sn:
            if sn_ not in dfsc.min_paths:
                continue
            d[sn_] = dfsc.min_paths[sn_][0].cost()

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
            sx = measures.zero_div(sum(dx.values()),len(dx),float('inf')) 
            d[n] = sx
        return d 

    """
    calculates distance of node `n` to each 
    of the sec nodes. 

    return:
    - dict, sec node idn -> distance
    """
    def secnode_distance_map(self,n,is_sec:bool=True):

        if n not in self.resource_sg.d:
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

    ###################################################
    ######## methods used in <TDir> nodepath deltas. 
    
    # TODO: un-used,test 
    """
    NOTE: class instances of <TDirector> do not remember
          `locset` as destination nodes for navigation.

    locset := set, node locations that may or may not
              be in the `resource_sg`.

    return:
    - float, in [0.,1.] such that 1. signifies the performance 
             as better than all previous, 0. signifies performance 
             as worse than all previous. 
    """
    def reflect_on_target_performance(self,locset,f=np.min):
        assert f in DEFAULT_TDIRECTOR_OBJ_FUNCTIONS
        dfsc = self.resource_sg.sp[self.loc()]

        def fx():
            distances = []
            for l in locset:
                if l not in dfsc.min_paths:
                    continue 
                q = dfsc.min_paths[0].cost()
                distances.append(q) 
            distances = np.array(distances)
            if len(distances) == 0: return None
            return f(distances)

        if self.vp() == "I":
            assert len(locset) >= 1
        else: 
            assert len(locset)  == 1

        sx = fx()
        if type(sx) == None: return None 
        self.ps = np.append(self.ps,sx)
        while len(self.ps) > self.tdts:
            self.ps.pop(0)        
        return self.reflect_on_performance()

    """
    required method of calculation by class 
    method above. 
    """
    def reflect_on_performance(self):
        q = self.ps[-1]
        r = 0 
        for i in range(len(self.ps) - 1): 
            if q > self.ps[i]:
                r += 1
        return measures.zero_div(r,len(self.ps) - 1,0.0)

    ################################################
    
    """

    Default path selection to a SEC node. 

    if no SEC node exists for path, outputs None, 
    otherwise outputs the path to the nearest 
    SEC node. 

    Uses `rnd_struct` as a selector in the
    event of ties. 
    """
    def defsec_path(self,rnd_struct):
        # case: 

        q = []
        if type(self.resource_sg) != SNGraphContainer:
            print("NADA")
            return None 

        dfsc = self.resource_sg.sp[self.loc()]
        assert type(dfsc) == DFSCache

        for qs in self.resource_sg.sn:
            if qs not in dfsc.min_paths:
                print("?missing path?")
                continue
            px = dfsc.min_paths[qs][0]
            ex = (qs,px)
            q.append(ex)

        if len(q) == 0: 
            return None
        
        s = sorted(q,key=lambda xr:xr[1].cost())
        q = [s.pop(0)]
        stat = True
        while stat:
            if len(s) == 0: 
                stat = False
                continue 

            if s[0][1].cost() == q[0][1].cost(): 
                q.append(s.pop())
            else:
                stat = not stat
        
        qi = rnd_struct.randrange(0,len(q))
        return deepcopy(q[qi][1].invert()) 

    # TODO: test 
    def default_crackling_pathdec(self,predicted_distance:int=None,\
        rnd_struct=random,objf=np.min):
        ## NOTE: open-ended question here. 
        if not (type(predicted_distance) == int and predicted_distance >= 0):
            g_rad = self.resource_sg.radius()
            if g_rad < 1:
                predicted_distance = 1
            else:  
                predicted_distance = rnd_struct.randrange(1,g_rad+1) 

            predicted_distance = 0


        print("CRCK PRED-DIST ",predicted_distance)

        cs = self.check_radar()
        if len(cs) == 0: return None
        assert len(cs) == 1

        # spot the security of I's location
        cs = cs.pop()
        stat = cs in self.resource_sg.sn 

        dfsc = self.resource_sg.sp[self.loc()]

        # decision: if location is SEC, then 
        #           travel to it. Otherwise,
        #           select an NSEC for site of 
        #           interception.
        if stat:
            print("CRCK YES SEC")
            if cs not in dfsc.min_paths:
                return None
            return dfsc.min_paths[cs][0].invert()

        # run analysis on all nodes
        xs = self.targetnode_analysis(objf)
        assert len(xs) > 0
        print("TARGETNODE ANALYSIS")
        print(xs) 
        xs = [(k,abs(v - predicted_distance)) for k,v in xs.items()]
        nx_ = random_tiebreaker(xs,rnd_struct,False)         
        nx = nx_[0] 
        print("\tNODE: ",nx) 
        print("\t\t\t&&& &&& &&&\t\t\t")
        if nx not in dfsc.min_paths: return None
        return dfsc.min_paths[nx][0].invert()

    def load_path_to_node(self,loc):
        p = self.td.load_path_(self.resource_sg,loc)
        assert type(p) != type(None)
        self.load_new_path(p)

    # TODO: test
    """
    outputs 
    """
    def open_info_pathdec(self,rnd_struct):

        if len(self.open_info_var) == 0:
            return None

        # check `open_info_var` size based on `vp()`
        v = self.vp()
        if v == "C":
            assert len(self.open_info_var) == 1

        # collect all relevant nodes
        ns = set() 
        for (info_type,idn,value) in self.open_info_var:
            if info_type == 0: continue 

            if info_type == 2: 
                ns = ns | {value}
                continue
            ns_ = self.open_info_velocity_prediction((info_type,idn,value))
            ns = ns | ns_

        # select the (node,path) based on `ns` and `vp()`
        dfsc = self.resource_sg.sp[self.loc()]

        # case: Crackling
        if v == "C":
            # case: certain
            if len(ns) == 1:
                ns_ = ns.pop()
                if ns_ not in dfsc.min_paths:
                    return None
                return dfsc.min_paths[ns_][0].invert()
            
            # case: not certain
            #       use rnd_struct
            ns = sorted(list(ns))
            i = rnd_struct.randrange(0,len(ns))
            return dfsc.min_paths[ns[i]].invert() 

        # case: IsoRing
        #       calculate the nodes to avoid
        avoid = ns - self.resource_sg.sn
        c = self.td.default_node_analysis()
        c_ = [(k,v) for k,v in c.items() if k not in avoid] 

            # subcase: none, revert to default
        if len(c_) == 0:
            c_ = [(k,v) for k,v in c.items()]
    
        l = random_tiebreaker(c,rnd_struct,max)[0]
        return dfsc.min_paths[l].invert() 

    # TODO: test 
    """
    arguments:
    - one_open_info_sample := (info_type,idn,value)

    return:
    - set, probable nodes that the complementary agent/s
           will be on at the end of a timespan.
    """
    def open_info_velocity_prediction(self,one_open_info_sample):
        assert one_open_info_sample[0] == 1
        q = self.resource_sg.ring_locs if \
            self.vp() == "C" else \
            self.resource_sg.crackling_locs
        idn = one_open_info_sample[1] 
        assert idn in q

        f = lambda x: x if self.vp() == "C" else x[0]
        x0 = f(q[idn])

        # use `resource_sg` to determine all possible 
        # nodes for next
        dfsc = self.resource_sg[x0]
        nsx = dfsc.nodeset_of_distance_d(\
            one_open_info_sample[2])
        return nsx
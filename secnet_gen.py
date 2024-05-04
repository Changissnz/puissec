from dfs_struct import * 
from sec_seq import * 

"""
serves as a descriptor for a methodology that 
<SecNetFrameGen> uses to generate the `frame` 
of a <SecNet>, meaning its nodes and edges. 
"""
class SecNetGenScheme:

    def __init__(self,description,rnd_seed,additional_args=None):
        self.description = description
        self.rnd_seed = rnd_seed
        self.additional_args = additional_args
        self.check_args()
        return

    def check_args(self):
        assert self.description in {"spine frame", "pairing frame","pseudo random"}
        assert type(self.rnd_seed) == int
        if self.description == "pseudo random":
            assert type(self.additional_args) == float
            assert self.additional_args >= 0.0 and self.additional_args <= 1.0

"""
generates frame according to instructions given by 
<SecNetGenScheme>. 
"""
class SecNetFrameGen:

    def __init__(self,sec_nodevec,nsec_nodevec,sngs):
        assert type(sngs) == SecNetGenScheme
        assert len(sec_nodevec) == len(set(sec_nodevec))
        assert len(nsec_nodevec) == len(set(nsec_nodevec))

        self.sec_nodevec = sec_nodevec
        self.nsec_nodevec = nsec_nodevec
        self.sngs = sngs
        self.init_rand()

        self.node_components = []
        self.d = defaultdict(list)
        self.node2sec_distance_map = defaultdict(defaultdict)
        self.init_search()
        return

    @staticmethod
    def generate_node_idns_by_sec_partition(num_sec,num_nsec):
        return -1

    def init_rand(self):
        random.seed(self.sngs.rnd_seed) 

    def init_search(self):
        self.node_components.extend([\
            set([n]) for n in self.sec_nodevec])
        
        for x in self.nsec_nodevec:
            self.node2sec_distance_map[x] = defaultdict(int)
        for x in self.sec_nodevec:
            self.node2sec_distance_map[x] = defaultdict(int)
        
        if self.sngs.description == "spine frame":
            self.init_spine_frame() 
            return

        self.node_components.extend([\
            set([n]) for n in self.nsec_nodevec])

        self.max_connections()

    def is_secnode(self,n):
        if n in self.sec_nodevec:
            return True
        if n in self.nsec_nodevec:
            return False
        return -1

    def node_to_component_index(self,node):
        for (i,x) in enumerate(self.node_components):
            if node in self.node_components[i]:
                return i
        return -1

    def component_index_to_d(self,component_index):
        nodes = self.node_components[component_index]
        d = defaultdict(set)
        for n in nodes:
            q = nodes.intersection(self.d[n])
            d[n] = q
        return d

    def component_nsec_node_distance_map(self,component_index,score_type=1):
        assert score_type in {1,2}
        q = self.node_components[component_index]
        dx = defaultdict(int)
        fx = self.nsec_node_distance_score if score_type == 1 \
                else self.nsec_node_distance_score_type2

        for q_ in q:
            dx[q_] = fx(q_)
        return dx

    ## NOTE: needs to be re-written 
    """
    distances to secnodes - degree of node
    """
    def nsec_node_distance_score(self,nsec_node):
        q = set(self.node2sec_distance_map[nsec_node].keys()) \
            - set(self.nsec_nodevec)
        l1 = 0
        for k in q:
            if self.is_secnode(k):
                l1 += self.node2sec_distance_map[nsec_node][k]
        l2 = len(self.d[nsec_node])
        return l1 - l2 

    """
    min distance to secnodes - degree of node 
    """
    def nsec_node_distance_score_type2(self,nsec_node):
        q = set(self.node2sec_distance_map[nsec_node].keys()) \
            - set(self.nsec_nodevec)
        l1 = []
        for k in q:
            if self.is_secnode(k):
                l1.append(self.node2sec_distance_map[nsec_node][k])
        if len(l1) == 0:
            l1 = -1 * float("inf") if nsec_node in self.nsec_nodevec else 0 
        else: 
            l1 = sum(l1) / len(l1)
        l2 = len(self.d[nsec_node])
        return l1 - l2 

    def connectivity(self):
        d = 0
        for v in self.d.values():
            d += len(v)
        d = int(d / 2)
        return d / self.max_conn 

    def max_connections(self):
        self.max_conn = 0 
        l = len(self.sec_nodevec) + len(self.nsec_nodevec)
        for i in range(l):
            self.max_conn += i 

    ############## next connection by mode

    def construct_frame(self):
        while self.next_conn():
            continue
        return

    def next_conn(self):
        if self.sngs.description == "spine frame":
            return self.next_conn_sf() 

        if self.sngs.description == "pairing frame":
            return self.next_conn_pf()

        if self.sngs.description == "pseudo random":
            return self.next_conn_pr()
        return False 

    def make_connection(self,node1,node2):
        self.d[node1].append(node2)
        self.d[node2].append(node1)

        self.node2sec_distance_map[node1][node2] = 1
        self.node2sec_distance_map[node2][node1] = 1 
        return

    ############### next connection: pairing frame 
    ############################################################

    def next_conn_pf(self):
        pc = self.pair_of_components()

        if type(pc) == type(None):
            return False

        pf1 = [(n,self.nsec_node_distance_score_type2(n)) for n in self.node_components[pc[0]]]
        random.shuffle(pf1)

        pf2 = [(n,self.nsec_node_distance_score_type2(n)) for n in self.node_components[pc[1]]]
        random.shuffle(pf2)

        n1 = sorted(pf1,key=lambda x: x[1])[-1][0]
        n2 = sorted(pf2,key=lambda x: x[1])[-1][0]
        self.df_update_distance__type_bicomponent(pc[0],pc[1],\
            n1,n2)
        return True

    """
    fetches indices of two components C1 and C2 such that 
    C1 has a sec node and C2 has an nsec node 
    """
    def pair_of_components(self):
        if len(self.node_components) == 1:
            return None

        # search for the sec node
        index1 = -1
        for i in range(len(self.node_components)): 
            if len(self.node_components[i]) != 1: 
                continue 
            x = next(iter(self.node_components[i]))
            if x in self.sec_nodevec:
                index1 = i
                break

        # search for the component with nsec node
        qs = [(i,sum(set(self.component_nsec_node_distance_map(i,score_type=2).values()))) \
            for (i,c) in enumerate(self.node_components)]
        random.shuffle(qs)
        qs = sorted(qs,key=lambda x: x[1])
        if index1 == -1:
            index1 = qs[0][0]
            index2 = qs[1][0]
            return index1,index2

        index2 = qs[0][0]
        return index1,index2

    ######### next connection: spine frame method
    ############################################################

    def init_spine_frame(self):
        self.node_components.append(set(deepcopy(self.nsec_nodevec)))
        l  = len(self.nsec_nodevec)
        for i in range(l - 1):
            n0 = self.nsec_nodevec[i]
            n1 = self.nsec_nodevec[i+1]
            self.d[n0].append(n1)
            self.d[n1].append(n0)

            for j in range(i+1,l):
                n2 = self.nsec_nodevec[j] 
                self.node2sec_distance_map[n0][n2] = j - i
                self.node2sec_distance_map[n2][n0] = j - i
        self.spine_counter = 0 
        return

    def next_conn_sf(self):
        sn = self.find_next_secnode_component()
        if sn == -1: return False
        sn_node = next(iter(self.node_components[sn]))
        component_index,min_node = self.next_conn_sf_()
        l = len(self.node_components) - 1
        self.df_update_distance__type_bicomponent(\
            sn,component_index,sn_node,min_node)
        return True 

    """
    makes next connection for spine frame
    """
    def next_conn_sf_(self):

        component_index = len(self.node_components) - 1
        dx = self.component_nsec_node_distance_map(component_index,score_type=1)
        dx = [(k,v) for (k,v) in dx.items()]
        random.shuffle(dx)
        
        f = min if self.spine_counter == 0 else max
        self.spine_counter = (self.spine_counter + 1) % 2
        min_node = f(dx,key=lambda x: x[1])
        return component_index,min_node[0] 

    """
    finds the first component that has 1 node
    and the node is a secnode
    """
    def find_next_secnode_component(self):
        component_index = -1 
        for i in range(len(self.node_components)):
            if len(self.node_components[i]) != 1: continue
            n = next(iter(self.node_components[i]))
            if self.is_secnode(n):
                component_index = i
                break 
        return component_index

    ######### next connection: pseudo-random frame method
    ############################################################

    def next_conn_pr(self):
        conn = self.connectivity() 
        ##print("next conn, pr: ", conn)
        nodepair = self.pr_nodepair() 

        stat1 = conn >= self.sngs.additional_args
        stat2 = type(nodepair) == type(None) 

        if stat1 or stat2:
            ## ?? 
            ##print("number of components: ", len(self.node_components))
            for i in range(len(self.node_components)):
                ##print("updating component {}".format(i))
                self.df_update_distance__type_component(i)
            return False 

        # case: connectivity level satisfied,finalize by updating
        #       distances for all 
        """
        if conn >= self.sngs.additional_args:
            return False 

        nodepair = self.pr_nodepair() 
        print("node-pair: ",nodepair)
        if type(nodepair) == type(None): 
            return False
        """

        ci1 = self.node_to_component_index(nodepair[0])
        ci2 = self.node_to_component_index(nodepair[1])

        ## ??
        # merge components
        if ci1 != ci2:
            ci = sorted([ci1,ci2])[::-1]
            q = self.node_components.pop(ci[0])
            q = q | self.node_components.pop(ci[1])
            self.node_components.append(q) 

        self.make_connection(nodepair[0],nodepair[1])
        return True 
        
        ###########################
        
        ci1 = self.node_to_component_index(nodepair[0])
        ci2 = self.node_to_component_index(nodepair[1])

        # case 1: uni-component
        if ci1 == ci2:
            self.make_connection(nodepair[0],nodepair[1])
            self.df_update_distance__type_component(ci1)
        # case 2: bi-component 
        else: 
            self.df_update_distance__type_bicomponent(ci1,\
                ci2,nodepair[0],nodepair[1])
        return True 

    def pr_nodepair(self):
        q = self.nsec_nodevec + self.sec_nodevec
        qx = [q_ for q_ in q if not self.maximally_connected_node(q_)]

        if len(qx) == 0:
            return None

        # randomly choose the first qualifying node 
        i = random.randrange(0,len(qx))
        n = qx.pop(i)

        # randomly choose the second qualifying node
        qx2 = list(set(qx) - set(self.d[n]))
        j = random.randrange(0,len(qx2))
        n2 = qx2.pop(j)
        return n,n2 

    def maximally_connected_node(self,node):
        l = len(self.d[node])
        return l == len(self.nsec_nodevec) + len(self.sec_nodevec) - 1  


    ######### distance update functions  
    ############################################################

    """
    updates distance for a component C after an edge has 
    been added to it, in which the nodeset of C remains 
    identical before and after the edge addition. 
    """
    def df_update_distance__type_component(self,component_index,path_cost_func= CUMULATIVE_PATH_COST_FUNC):
        d = self.component_index_to_d(component_index)
        prev_nodes = []
        c = self.node_components[component_index]
        for n in c: 
            ##print("\t\tupdating node ",n)
            dfsc = DFSCache(n,deepcopy(d),edge_cost_function=CUMULATIVE_EDGE_COST_FUNCTION)
            ##print("\t\texec dfs")
            dfsc.exec_DFS()
            ##print("\t\tstoring min paths")
            dfsc.store_minpaths(ns=set(c) - set(prev_nodes),cost_func = path_cost_func)
            self.update_paths_by_dfsc(dfsc,path_cost_func)
            prev_nodes.append(n) 
        return

    """
    updates the minimum distances of the map 
    `node2sec_distance_map` using the <DFSCache> 
    argument. 
    """
    def update_paths_by_dfsc(self,dfsc,path_cost_func= CUMULATIVE_PATH_COST_FUNC):
        head = dfsc.start_node
        ##print("head node: ",head)
        for (k,v) in dfsc.min_paths.items():
            """
            print("UPDATING")
            for v_ in v:
                print("\t cost: ", v_.cost(path_cost_func))
            """
            d = v[0].cost(path_cost_func)
            self.node2sec_distance_map[k][head] = d
            self.node2sec_distance_map[head][k] = d
        return

    """
    updates the distances between two components that 
    become connected by the edge `c1-c2`. 
    """
    def df_update_distance__type_bicomponent(self,c1_index,c2_index,\
        c1_node,c2_node):

        comp1 = self.node_components[c1_index]
        comp2 = self.node_components[c2_index]
        assert c1_node in comp1 and c2_node in comp2

        # add the edge `c1-c2`
        self.make_connection(c1_node,c2_node)

        for x in comp1:
            # get distance of x to c1_node
            d1 = self.node2sec_distance_map[x][c1_node]

            for y in comp2:
                # get distance of y to c2_node
                d2 = self.node2sec_distance_map[y][c2_node]
                d3 = d1 + d2 + 1
                self.node2sec_distance_map[x][y] = d3
                self.node2sec_distance_map[y][x] = d3

        q = sorted([c1_index,c2_index])[::-1]
        cx = comp1 | comp2 
        self.node_components.pop(q[0])
        self.node_components.pop(q[1]) 
        self.node_components.append(cx) 
        return

"""
generates the connections for Pr. dependencies 
of each <Sec> instance in `sec_seq`. See the 
class method `make_conn` for more information.

For each <Sec> S, its (co?)-dependency map has elements
of the form 
```
key: 
    "index of the local optimum of S",
    "index of other <Sec> S2"."index of local optimum of S2"
value:
    f in [0.,1.]
```

For example, for two secrets <S1,S2> with 5 and 6 local 
optima,respectively: 
```
4,1.3
```
is the key for the fifth local optimum of S1 with the 
 fourth local optimum of 2nd secret S2. 

The co-dependency map generated is reflexive. 
"""
class SecNetDepGen:

    def __init__(self,sec_seq,rnd_struct,\
        min_components:int,max_nconn_ratio:float,
        dlength_range,depconn_ratio=None,\
        conn_candidate_size=float('inf')):
        assert type(sec_seq) == list
        for s in sec_seq: assert type(s) == Sec 
        assert type(min_components) == int and min_components > 0
        assert len(sec_seq) >= min_components
        assert max_nconn_ratio >= 0.0 and max_nconn_ratio <= 1.0

        assert len(dlength_range) == 2
        assert dlength_range[0] < dlength_range[1]
        assert type(dlength_range[0]) == type(dlength_range[1])
        assert type(dlength_range[0]) == int
        assert dlength_range[1] <= len(sec_seq) 
        assert type(depconn_ratio) in {type(None),float}

        self.sq = sec_seq 
        self.clear_sec_pr_maps()

        self.rnd_struct = rnd_struct 
        self.min_components = min_components
        self.max_nconn_ratio = max_nconn_ratio
        # TODO: barely used, delete? 
        self.dlength_range = dlength_range 
        self.depconn_ratio = depconn_ratio
        # used to narrow the selection space 
        # for quicker processing of depconn. 
        self.conn_candidate_size = conn_candidate_size

        # each element is a set representing a 
        # co-dependent component
        self.cd_comp_sets = []
        # index of <Sec> -> indices of other <Sec> 
        #               dependent on it.
        self.dep_map = defaultdict(list)

        for i in range(len(sec_seq)):
            self.cd_comp_sets.append(set([i])) 
            self.dep_map[i] = [] 
        
        # every element is of the form
        #   key= node idn.
        #   value= dict
        #       key=connection b/t node optimum and other
        #       value=strength of the connection
        self.dep = defaultdict(defaultdict)
        self.codep = defaultdict(defaultdict)
        return

    def clear_sec_pr_maps(self):
        for s in self.sq: 
            s.dm.clear()
            s.cdm.clear()
        return 

    def assign_conn(self,size_limt=float('inf'),\
        l=[1,2,3]):
        stat = True
        i = 0
        while stat: 
            stat = self.make_conn(deepcopy(l))
            stat = stat and i < size_limt
            i += 1
        self.write_conn_to_Sec()
        return

    def fstat(self):

        # d-length satisfaction
        dstat = max([len(v) for v in self.dep_map.values()])
        dstat = dstat >= self.dlength_range[0] and dstat <= self.dlength_range[1] 

        # number of components 
        cstat = self.available_C2C_conn()

        # max_nconn_ratio
        r = max([self.codep_conn_of_node(i) for i in range(len(self.sq))])
        rstat = r <= self.max_nconn_ratio
        return dstat,cstat,rstat 

    ###################### make conn. method. 
    def make_conn(self,options = [1,2,3]):
        assert set(options).issubset({1,2,3})
        if len(options) == 0:
            return False

        # choose a random element in options
        i = self.rnd_struct.randrange(0,len(options))
        o = options.pop(i)

        stat = True
        # dependency conn.
        if o == 1:
            stat = self.make_dep_conn()
        # co-dep. conn. (b/t two components)
        elif o == 2: 
            stat = self.make_codep_C2C_conn()
        # co-dep. conn. (in one component) 
        else:  
            stat = self.make_codep_CInC_conn()
        ##print("stat for {}:{}".format(o,stat))
        if type(stat) != type(None):
            return True
        return self.make_conn(options)
        
    def make_dep_conn(self):
        # get list of candidates
        l = self.available_for_dependency(self.conn_candidate_size)
        ##print("# available for dep.: ",len(l))
        if len(l) == 0: 
            ##print("NADA")
            return None

        # choose a pair of non-connected nodes
        q = self.rnd_struct.randrange(0,len(l))
        l = list(l)
        x = l[q]

        self.dep_conn_process(x) 
        return True

    def dep_conn_process(self,pn):
        c = None
        if type(self.depconn_ratio) == type(None):
            c = 1
        else: 
            assert self.depconn_ratio >= 0.0 and self.depconn_ratio <= 1.0

            x = len(self.sq[pn[0]].opm) * len(self.sq[pn[1]].opm)
            c = int(math.ceil(self.depconn_ratio * x)) 
        
        ##print("making {} dep.".format(c))
        for i in range(c):
            ##print("-- {}".format(i))
            self.one_dep_conn(pn)

    """
    adds an optimum-to-optimum dependency connection 
    so that pn[1] depends on pn[0].
    """
    def one_dep_conn(self,pn):
        # get the available optima indices for each
        # of the nodes w.r.t. each other
        dmx = self.dep[pn[1]]
        sz0 = len(self.sq[pn[1]].opm)
        sz1 = len(self.sq[pn[0]].opm)
        l1_,l2_ = available_for_dep(dmx,pn[0],sz0,sz1)

        # choose a local optima from each of the elements
        i1 = l1_[self.rnd_struct.randrange(0,len(l1_))]
        i2 = l2_[self.rnd_struct.randrange(0,len(l2_))]
        prv = round(self.rnd_struct.uniform(0.,1.),5)
        self.add_dependency(pn[0],pn[1],i1,i2,prv)

    def make_codep_C2C_conn(self,attempts=10):
        if attempts <= 0: 
            return None 

        stat = self.available_C2C_conn()
        if not stat: return None

        cind = [i for i in range(len(self.cd_comp_sets))]
        # choose a (node,component)
        ci1 = self.rnd_struct.randrange(0,len(cind)) 
        ci1 = cind.pop(ci1) 
        x = list(self.cd_comp_sets[ci1])
        ci1_ = self.rnd_struct.randrange(0,len(x))
        n1 = x[ci1_]

        # get available nodes for n 
        avail = self.available_for_codep_with_node(n1)
        if len(avail) == 0:
            return self.make_codep_C2C_conn(attempts-1)
        q = self.rnd_struct.randrange(len(avail))
        n2 = list(avail)[q]

        # make a co-dep b/t n1 and n2 
        prv = round(self.rnd_struct.uniform(0.,1.),5)
        self.add_codep(n1,n2,prv)
        return True 

    def make_codep_CInC_conn(self):
        q = self.available_CInC_conn()
        x = len([q_ for q_ in q if len(q_) > 1])
        if x == 0: return None

        n1,n2 = None,None
        for (i,q_) in enumerate(q):
            if len(q_) < 2: continue

            q2 = list(q_) 

            ix = self.rnd_struct.randrange(0,len(q2))
            n1 = q2.pop(ix) 
    
            ix = self.rnd_struct.randrange(0,len(q2))
            n2 = q2.pop(ix) 
            break

        if type(n1) == type(None) or type(n2) == type(None):
            return None

        prv = round(self.rnd_struct.uniform(0.,1.),5)
        self.add_codep(n1,n2,prv)
        return True 

    def add_codep(self,n1:int,n2:int,prv:float): 
        # choose two local optima, one 
        # from each of n1,n2
        ps1 = self.sq[n1].optima_points().shape
        ps2 = self.sq[n2].optima_points().shape

        i1 = self.rnd_struct.randrange(0,ps1[0]) 
        i2 = self.rnd_struct.randrange(0,ps2[0]) 

        s1 = str(i1) + "," + str(n2) + "." + str(i2)
        s2 = str(i2) + "," + str(n1) + "." + str(i1)
        self.codep[n1][s1] = prv
        self.codep[n2][s2] = prv 

        c1 = self.component_of_node(n1)
        c2 = self.component_of_node(n2)
        if c1 != c2: 
            q = sorted([c1,c2])[::-1]
            x1 = self.cd_comp_sets.pop(q[0])
            x2 = self.cd_comp_sets.pop(q[1])
            self.cd_comp_sets.append(x1 | x2) 

    def available_optimapair_for_nodepair(self,n1,n2,map_type='c'):
        assert map_type in {'c','d'}


        return -1 

    ################ for co-dependency

    def available_C2C_conn(self):
        l = len(self.cd_comp_sets)
        return l > self.min_components

    def available_CInC_conn(self):
        # first, get each node's codep conn.
        # measure 
        ms = []
        for i in range(len(self.sq)):
            ms.append(self.codep_conn_of_node(i))

        # get the open Sec instances
        ms_ = [i for (i,m_) in enumerate(ms) if m_ < \
            self.max_nconn_ratio]

        # group the open Sec instances based 
        # on component
        qs = [set() for i in range(len(self.cd_comp_sets))]
        for si in ms_: 
            j = self.component_of_node(si)
            assert j != -1
            qs[j] = qs[j] | {si} 
        return qs 

    def codep_conn_of_node(self,n):
        q = self.codep[n]

        j = self.component_of_node(n) 
        if j == -1: return 0.0

        q2 = self.cd_comp_sets[j] - {n}
        sm = sum([len(self.sq[q_].opm) for q_ in q2])
        sm = sm * len(self.sq[n].opm)

        if sm == 0.0: return 0.
        return round(len(q)/ sm,5)

    def component_of_node(self,n):
        j = -1
        for (i,c) in enumerate(self.cd_comp_sets):
            if n in c:
                j = i
                break 
        return j 
        
    #################################################

    """
    adds a connection b/t `ref` and `dependent`,
    and updates all dependents of `dependent` w/
    `ref`. 
    """
    def add_dependency(self,ref,dependent,i1,i2,prv):
        ##print("adding dependency")
        ##print("{}.{}-->{}.{}".format(ref,i1,dependent,i2))

        assert ref not in self.dep_map[dependent]
        qx = deepcopy(self.dep_map[ref]) 
        ##print(qx) 
        q = [dependent] 
        q = q + self.dep_map[dependent]
        q1 = deepcopy(q) 
        self.dep_map[ref].extend(q)
        self.dep_map[ref] = list(set(self.dep_map[ref]))
        for q_ in qx:
            self.dep_map[q_].extend(q1)
            self.dep_map[q_] = list(set(self.dep_map[q_]))

        if dependent not in self.dep: 
            self.dep[dependent] = defaultdict(float)
        
        k = str(i2) + "," + str(ref) + "." + str(i1) 
        self.dep[dependent][k] = prv
        return

    """
    return:
    - set, each is a tuple (master node,dependent node)
    """
    def available_for_dependency(self,sz_limit=float('inf')):
        l = set()

        for i in range(len(self.sq)):
            sx = self.available_for_dependency_on_node(i)
            qx = set([(i,sx_) for sx_ in sx])
            for qx_ in qx: 
                l = l | {qx_} 
                if len(l) >= sz_limit: 
                    break 

        return l

    """
    for a node n2 to be able to depend on 
    `n`, node `n` cannot be dependent on n2
    and also cannot have any codependencies 
    with n2. 
    """
    def available_for_dependency_on_node(self,n,verbose=False):
        q = set([i for i in range(len(self.sq))])
        ##print("L: ",len(q))
        already = self.dep_map[n]

        # get those that n depends on
        qs = []
        for (k,v) in self.dep_map.items():
            if k == n: continue
            if n in v: qs.append(k)
        qs = set(qs) 

        if verbose:
            print("N={} depends on \n\t{}".format(n,qs))
            print("A-set: {}".format(already))

        ##print("L: ",len(qs))

        # no dependencies or duplicates. 
        q = q - {n}
        q = q - set(already)
        q = q - qs
        ##print("L: ",len(q))

        # no codependencies
        j = self.component_of_node(n)
        q = q - self.cd_comp_sets[j] 
        ##print("L: ",len(q))

        # iterate through candidates and 
        # see which ones cannot due to 
        # contradiction
        qsx = []
        for q_ in q: 
            qr = self.dep_map[q_]
            sx = set(qr).intersection(qs) 
            if len(sx) > 0: 
                if verbose: 
                    print("intersection w/ {}: {}".format(qs,qr))
                continue
            qsx.append(q_)
        qsx = set(qsx)
        
        ##print("L: ",len(qsx))
        ##print("available for dep on node {}:{}".format(n,qsx))
        return qsx

    ## ?? 
    def available_for_codep_with_node(self,n):

        q = set([i for i in range(len(self.sq))])
        q = q - set(self.dep_map[n])

        for (k,v) in self.dep_map.items():
            if n == k: continue 
            if n in v: 
                q = q - {k}
        return q 

    def write_conn_to_Sec(self):
        for (i,s) in enumerate(self.sq):
            s.dm = self.dep[i]
            s.cdm = self.codep[i]
        return

    #################################

    """
    calculates the codependent nodes of 
    `n` if `is_dep`=False, otherwise calculates
    the nodes that n depends on.
    """
    def connected_nodes(self,n,is_dep:bool=True): 
        dx = self.dep if is_dep else self.codep 
        q = dx[n]
        s = set()
        for k,v in q.items():
            k_ = parse_dconn(k)
            s = s | {k_[1]} 
        return s 



def Sec_list_sample1(): 
    random.seed(12)
    np.random.seed(12)

    singleton_range = [0.,1.] 
    dimension = 5
    num_optima = 2
    countermeasure = (0.6,0.5) 
    secs = []

    for i in range(5): 
        sec = Sec.generate_bare_instance(singleton_range,dimension,num_optima,\
        countermeasure,rnd_struct=np.random)
        secs.append(sec)
    return secs 

def Sec_list_sample2(num_secs=12): 
    random.seed(14)
    np.random.seed(19)

    singleton_range = [0.,1.] 
    num_optima = 12 
    dimension = 5
    countermeasure = (0.7,0.5) 
    secs = []

    for i in range(num_secs): 
        sec = Sec.generate_bare_instance(singleton_range,dimension,num_optima,\
        countermeasure,rnd_struct=np.random)
        secs.append(sec)
        ##print("sec {}".format(i))
    return secs 

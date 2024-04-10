from dfs_struct import * 

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

"""
<Sec> reference map used to compare <Cracker> 
solution candidates to actual solution or 
assign each of those solution candidates a Pr. 
value for their being the solution. 
"""
from sec_mapr import * 

######################## filter-functions for 
######################## sequences of possible-decision maps. 
######################## Used in the frequency-counting phase 
######################## of <SRefMap>.

def filterf_idnmap(prmap_seq): 
    return prmap_seq

def filterf_ext_(prmap_seq,extf):
    assert extf in {min,max}
    if len(prmap_seq) == 0: return [] 

    # retrieve the pertinent Pr. idn 
    # for each of the maps
    prv = [] 
    for (i,p) in enumerate(prmap_seq):
        x = (i,p[i])
        prv.append(x) 

    # get the extremum
    extr = extf(prv,key=lambda q:q[1])
    return [prmap_seq[extr[0]]] 

def generate_filterf_ext_instance(extf): 

    def f(prmap_seq): 
        return filterf_ext_(prmap_seq,extf) 
    return f

"""
the default binop. 
"""
def absdiff(x1,x2):
    return abs(x1 - x2) 

def filterf_ext_with_binary_op_accum_(prmap_seq,extf,binop): 
    l = len(prmap_seq)
    if l == 0: return 0 

    def binop_on_index(i): 
        x = prmap_seq[i]
        ks = set(x.keys()) - {i} 
        a = 0.0 
        for k in ks: 
            a += binop(x[i],x[k]) 
        return a 
    
    v = [(j,binop_on_index(j)) for j in range(l)]
    extr = extf(v,key=lambda x: x[1])
    return [prmap_seq[extr[0]]]

def generate_filterf_ext_with_binary_op_accum(extf,binop):

    def f(prmap_seq): 
        return filterf_ext_with_binary_op_accum_(prmap_seq,extf,binop)
    return f 

class SRefMap: 

    PRISM_VERTEX_LABELS = {'greedy-lone',\
            'greedy-d','greedy-c','greedy-dc',\
            'actual'}

    def __init__(self,opmn,dms,cdms,ndmaptype='cd'):
        assert ndmaptype in {'c','d','cd'}
        self.opmn = opmn 
        self.dms = dms
        self.cdms = cdms 
        self.ndmt = ndmaptype

        # vertex label -> vector value
        self.v = {}
        # node idn ->
        # identifier of optima in <Sec> `s` -> 
        # (min. possible-decision map,max. possible-decision map)
        self.preproc_map = None

        self.dcnt = defaultdict(Counter)
        self.preprocess() 
        return 

    def load_prism_vertices(self,d):
        for (k,v) in d.items():
            assert k in self.PRISM_VERTEX_LABELS
            self.v[k] = v
        return

    ############# pre-process function and calculating
    ############# possible-decision maps 

    def preprocess(self):
        self.preproc_on_seclist()
        return

    def reprocess(self,ndmt):
        assert ndmt in {'c','d','cd'}
        self.ndmt = ndmt
        self.preprocess()  

    """
    preprocess function on all <Sec> identifiers; 
    creates a map with min&max possible-decision 
    maps for every (<Sec> idn.,optima idn.) pair. 
    """
    def preproc_on_seclist(self):
        self.preproc_map = defaultdict(defaultdict)
        for k in self.opmn.keys():
            q = self.sec_to_pd_ext_map(k)
            self.preproc_map[k] = q 
        return 

    # TODO: test. 
    """
    return:
    - dict, identifier of optima in <Sec> `s` -> 
            (min. possible-decision map,max. possible-decision map)
    """
    def sec_to_pd_ext_map(self,s):
        assert s in self.opmn
        opt2pdchain = {}

        for k in self.opmn[s].keys():
            opt2pdchain[k] = self.greedy_pd_chain_ext(s,k,self.ndmt)
        return opt2pdchain

    ##################### greedy solutions using 
    ##################### extremum functions 

    """
    return: 
    - dict, dependent|codependent node of `node` -> 
            set of local opt. indices by ext. func. 
    """
    def greedy_pd_chain_ext(self,node,dec,map_types):
        assert map_types in {'c','d','cd'}

        if map_types == 'c': 
            return self.greedy_pd_chain_ext_(node,dec,False)
        
        if map_types == 'd': 
            return self.greedy_pd_chain_ext_(node,dec,True)

        dm1 = self.greedy_pd_chain_ext_(node,dec,True)
        dm2 = self.greedy_pd_chain_ext_(node,dec,False)
        dm1[0].update(dm2[0])
        dm1[1].update(dm2[1])
        return dm1 
        
    '''
    calculates the possible-decision chain 
    extremum (min,max) for `node` with `dec`. 
    '''
    def greedy_pd_chain_ext_(self,node,dec,is_dm=True):
        q = None 
        if is_dm:
            q = deepcopy(self.dms[node])
        else: 
            q = deepcopy(self.cdms[node])
        minny = extdec_dmap_set(q,dec,min)
        maxxy = extdec_dmap_set(q,dec,max)
        return minny,maxxy 

    #################### probability calculations for 
    #################### decision-chains. 

    # get the Pr. range of a node
    def pr_range_of_node_dec(self,node_idn,dec_idn,
        opmi,dm,decision_chain):
        prmap = dep_weighted_Pr_for_node_dec(node_idn,\
            dec_idn,opmi,dm,decision_chain)
        return 

    # TODO: test.
    '''
    frequency-counter process that implements
    a methodology for determining the best 
    (node,dec) map. 

    return: 
    - dict, sec idn -> index of selected opt.
    '''
    def fc_proc__best_nodedec_map(self,indices=[0,1]): 
        ## TODO: re-write 
        # identifier of optima in <Sec> `s` -> 
        self.dcnt = defaultdict(Counter)
        ks = list(self.opmn.keys())

        for ks_ in ks:
            self.extfc_proc_on_node(ks_,indices)

        d = {}
        for k,v in self.dcnt.items():
            v_ = sorted([(k2,v2) for k2,v2 in v.items()])
            if len(v_) == 0:
                d[k] = -1 
                continue 
            d[k] = max(v_,key=lambda x:x[1])[0]
        return d 

    # TODO: test.
    '''
    return:
    - 
    '''
    def extfc_proc_on_node(self,n,indices): 
        assert type(self.preproc_map) != type(None) 
        assert set(indices).issubset({0,1}) 

        q = self.preproc_map[n]
        for v in q.values(): 
            for i in indices: 
                self.dcnt = update_SRefMap_counter(self.dcnt,v[i])
        return 

    '''
    calculates the optima map for node `n` based on 
    its decision `dec` using function `F`, in which 
    `fi` is an index in {0,1} specifying which 
    possible-decision map to use (MIN|MAX). 
    '''
    def prmap_for_nodedec(self,n,dec,fi,pr_type):
        assert fi in {0,1}
        F = self.prfunc_by_type(pr_type)
        return F(n,dec,fi)

    # TODO: test.
    """
    outputs the Pr-function to use for `greedy-?` 
    functions. 
    """
    def prfunc_by_type(self,pr_type="greedy-lone"):
        f = None
        if pr_type in ['greedy-d','greedy-c']:
            def f(n,dec,fi):
                q = self.preproc_map[n]
                pdx = q[dec][fi]
                pmi = PDMapIter(pdx)
                ix = next(pmi)
                ix_ = [(k,v) for (k,v) in ix.items()]
                dm = self.dms[n] if pr_type == 'greedy-d' else self.cdms[n]
                return dep_weighted_Pr_for_node_dec(n,dec,self.opmn[n],\
                    dm,ix_)

        elif pr_type == 'greedy-lone': 
            def f(n,dec,fi):
                dm = self.opmn[n]
                return dm 

        elif pr_type == 'greedy-dc': 
            def f(n,dec,fi):
                f1 = self.prfunc_by_type('greedy-d')
                f2 = self.prfunc_by_type('greedy-c')
                m1 = f1(n,dec,fi)
                m2 = f2(n,dec,fi)

                c1 = Counter(m1)
                c2 = Counter(m2)
                c1 = c1 + c2 
                for k in c1.keys():
                    c1[k] = c1[k] / 2.0 
                return defaultdict(float,c1)
        return f 


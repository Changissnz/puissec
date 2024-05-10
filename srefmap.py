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


######################### methods for use with: 
######################### - decision and probability maps

# symmetric operations
mapdiff_discrete = lambda x1,x2: round(abs(x1 - x2),10) == 0 if \
                                type(x1) != type(None) and \
                                type(x2) != type(None) else False 

mapdiff_cont = lambda x1,x2: round(abs(x1-x2),5) if \
                    type(x1) != type(None) and \
                    type(x2) != type(None) else np.nan 

def mapdiff_(d1,d2,f):
    assert f in {mapdiff_cont,mapdiff_discrete}

    ks = set(d1.keys())
    ks = ks | set(d2.keys()) 

    q = defaultdict(bool if f \
        == mapdiff_discrete else \
        float) 
    for k in ks:
        v = d1[k] if k in d1 else None
        v2 = d2[k] if k in d2 else None
        q[k] = f(v,v2)
    return q 

################################################################

"""
<SRefMap> is a data structure that serves as a
solution-reference map. It outputs values based
on 
"""
class SRefMap: 

    PRISM_VERTEX_LABELS = {'greedy-lone',\
            'greedy-d','greedy-c','greedy-dc'} 

    def __init__(self,opmn,dms,cdms,ndmaptype='cd'):
        assert ndmaptype in {'c','d','cd'}
        self.opmn = opmn 
        self.dms = dms
        self.cdms = cdms 
        self.ndmt = ndmaptype

        # node idn ->
        # identifier of optima in <Sec> `s` -> 
        # (min. possible-decision map,max. possible-decision map)
        self.preproc_map = None

        # map structures that constitute prisms
        # for probability values and their associated
        # decisions made. 
        """
        key := str identifier for prism
        value := dict, sec idn. -> optima idn. -> PRMap(optima idn. -> Pr in [0,1])
        """
        self.prism_typePr = None
        """
        key := str identifier for prism
        value := dict,sec idn. -> optima idn.
        """
        self.prism_typeDec = None
        
        """
        comparison map for prism of type Pr. 

        Each key is a comparison-idn str.: 
            sec idn,opt. #1 idn,opt. #2 idn,F
        Each value is the map-difference of #1 and #2
            according to a function F, one of 
                {<mapdiff_cont>,<mapdiff_discrete>} 
        """
        self.prism_typePr_CMP = defaultdict(None)
        """
        comparison map for prism of type Dec. 

        Each key is a comparison-idn str.:
            x_0&x_1&F s.t. for x_i a 3-tuple,
            0 in {c,d,cd}
            1 in greedy-*
            2 in {0,1}.

        and F is identical to the one for 
        `prism_typePr_CMP`. 
        """
        self.prism_typeDec_CMP = None

        self.actual_dec = None
        self.actual_Pr = None

        self.dcnt = defaultdict(Counter)
        self.preprocess() 
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
    def fc_proc__best_nodedec_map(self,f=max,indices=[0,1]): 
        self.fc_proc__best_nodedec_map_(indices)

        d = {}
        for k,v in self.dcnt.items():
            v_ = sorted([(k2,v2) for k2,v2 in v.items()])
            if len(v_) == 0:
                d[k] = -1 
                continue 
            d[k] = f(v_,key=lambda x:x[1])[0]
        return d 

    def fc_proc__best_nodedec_map_(self,indices):
        ## TODO: re-write 
        # identifier of optima in <Sec> `s` -> 
        self.dcnt = defaultdict(Counter)
        ks = list(self.opmn.keys())

        for ks_ in ks:
            self.extfc_proc_on_node(ks_,indices)

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

    def first_decision(self):
        d = defaultdict(int)
        for k in self.opmn.keys():
            d[k] = {0,1} 
        return d 

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
                dx1 = deepcopy(self.opmn[n])

                q = self.preproc_map[n]
                deflt = self.first_decision() 
                pdx = q[dec][fi]
                if len(pdx) == 0: 
                    return dx1

                deflt.update(pdx) 
                pdx = deflt

                pmi = PDMapIter(pdx)
                ix = next(pmi)
                if type(ix) == type(None): 
                    return dx1

                ix_ = [(k,v) for (k,v) in ix.items()]
                dm = self.dms[n] if pr_type == 'greedy-d' else self.cdms[n]
                dx2 = dep_weighted_Pr_for_node_dec(\
                    n,dec,self.opmn[n],dm,ix_)
                dx1.update(dx2)
                return dx1 

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

    #################################################

    def build_prism(self):
        self.build_prism_('pr')
        self.build_prism_('dec')

    def build_prism_(self,prism_type):
        assert prism_type in {'pr','dec'}
        ptx = None 
        if prism_type == 'pr': 
            self.prism_typePr = defaultdict(defaultdict)
            ptx = self.prism_typePr
        else: 
            self.prism_typeDec = defaultdict(defaultdict)
            ptx = self.prism_typeDec

        map_types = ['c','d','cd']
        pr_types = ['greedy-lone','greedy-d','greedy-c',\
            'greedy-dc']
        F = [min,max] 
        pdmi = [0,1]
        pdmi2 = [(0,),(1,),(0,1)]

        X = None
        PF = None 
        if prism_type == 'pr':
            X = defaultdict(set,\
                {0:deepcopy(map_types),1:deepcopy(pr_types),\
                2:deepcopy(pdmi)})
            PF = self.collect_prism_points__PrMap
        else: 
            X = defaultdict(set,\
                {0:map_types,1:F,2:pdmi2})
            PF  = self.collect_prism_points__DecMap

        pdmi = PDMapIter(X)
        stat = True
        while stat:
            x = next(pdmi)
            stat = not (type(x) == type(None))
            if not stat: continue 
            prism = PF(x[0],x[1],x[2])
            key = []
            for i in range(3):
                key.append(str(x[i]))
            key = ",".join(key) 
            ptx[key] = prism 

    """

    return:
    - dict, sec idn. -> optima idn. -> PRMap(optima idn. -> Pr in [0,1])
    """
    def collect_prism_points__PrMap(self,map_type,pr_type,pdmi):
        assert pr_type in self.PRISM_VERTEX_LABELS 
        assert pdmi in {0,1} 

        def f(n):
            q = list(self.opmn[n].keys())
            d = defaultdict(defaultdict)
            for q_ in q: 
                prv = self.prmap_for_nodedec(n,q_,pdmi,pr_type)
                d[q_] = prv
            return d

        if map_type != self.ndmt: 
            self.reprocess(map_type)

        sk = set(self.opmn.keys())
        dx = defaultdict(None)
        for sk_ in sk: 
            dx[sk_] = f(sk_)
        return dx 

    """
    return:
    - dict, sec idn. -> optima idn.
    """
    def collect_prism_points__DecMap(self,ndmaptype,f,indices=[0,1]): 
        assert ndmaptype in {'c','d','cd'}
        assert f in {min,max}

        if ndmaptype != self.ndmt: 
            self.reprocess(ndmaptype)

        return self.fc_proc__best_nodedec_map(f,indices=indices)

    # TODO: test
    """
    """
    def binarycmp_prism_points__typePr(self,prism_key,\
        sec_idn,opt1_index,opt2_index,F):

        if prism_key not in self.prism_typePr: 
            return None

        d = self.prism_typePr[prism_key]
        d2 = d[sec_idn]
        dx1 = d2[opt1_index]
        dx2 = d2[opt2_index]

        m = mapdiff_(dx1,dx2,F)

        ox = sorted([opt1_index,opt2_index])
        key = prism_key + "&" + str(sec_idn) + "," +\
                ox[0] + "," + ox[1] + "&" + str(F)  
        self.prism_typePr_CMP[key] = m
        return

    # TODO: test
    """
    compares 2 dictionaries d1,d2 in Pr-map; comparative 
    measures used for introspection of 1 particular 
    """
    def binarycmp_prism_points__typeDec(self,prism_key1,prism_key2,F):
        assert F in {mapdiff_cont,mapdiff_discrete}

        d1 = self.prism_typeDec[prism_key1]
        d2 = self.prism_typeDec[prism_key2] 

        if type(d1) == type(None) or type(d2) == type(None): 
            return None

        m = mapdiff_(d1,d2,F) 

        key = prism_key1 + "&" + prism_key2 + "&" + str(F)
        self.prism_typeDec_CMP[key] = m
        return m 

    ##################################################################

    '''
    n2opt_map := dict, node idn. -> 
                        optima idn. -> 
                        Pr. value. 
    nd_map := dict, sec idn. -> optima index
    ''' 
    def load_actual(self,n2opt_map,nd_map): 
        self.actual_Pr = n2opt_map
        self.nd_map = nd_map 
        return


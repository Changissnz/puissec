"""
<Sec> reference map used to compare <Cracker> 
solution candidates to actual solution or 
assign each of those solution candidates a Pr. 
value for their being the solution. 
"""
from sec_mapr import * 

class SRefMap: 

    PRISM_VERTEX_LABELS = {'greedy-lone',\
            'greedy-d','greedy-c','greedy-dc',\
            'actual','conn-derived'}

    def __init__(self,opmn,dms,cdms):
        self.opmn = opmn 
        self.dms = dms
        self.cdms = cdms 

        self.ssm = secseq_map
        # vertex label -> vector value
        self.v = {}
        # node idn -> index of opt. -> 
        # identifier of optima in <Sec> `s` -> 
        # (min. possible-decision map,max. possible-decision map)
        self.preproc_map = None
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

        for k in self.opmn[s].items():
            opt2pdchain[k] = self.greedy_lone_pd_chain_ext(s,k,no_intersecting_keys=True)
        return opt2pdchain

    ##################### greedy solutions using 
    ##################### extremum functions 

    """
    return: 
    - dict, dependent|codependent node of `node` -> 
            set of local opt. indices by ext. func. 
    """
    def greedy_lone_pd_chain_ext(self,node,dec,no_intersecting_keys=True):

        dm1 = self.greedy_lone_pd_chain_ext_(node,dec,True)
        dm2 = self.greedy_lone_pd_chain_ext_(node,dec,False)

        if no_intersecting_keys:
            k1 = set(dm1[0].keys())
            k2 = set(dm1[1].keys())
            assert len(k1.intersection(k2)) == 0
        
        dm1[0].update(dm2[0])
        dm1[1].update(dm2[1])
        return dm1 
        
    '''
    calculates the possible-decision chain 
    extremum (min,max) for `node` with `dec`. 
    '''
    def greedy_lone_pd_chain_ext_(self,node,dec,is_dm=True):
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

    '''
    frequency-counter process. 
    '''
    def fc_proc(self): 
        return -1

    def fc_proc_on_node(self,i): 
        return -1 

    # get the Pr. range of a node
    def pr_range_of_node_dec(self,node_idn,dec_idn,
        opmi,dm,decision_chain):
        prmap = dep_weighted_Pr_for_node_dec(node_idn,\
            dec_idn,opmi,dm,decision_chain)
        return 


    """
    d := dict, sec idn -> local optima index
    """
    def pr_of_nodedec_(self,d,pr_type):
        return -1


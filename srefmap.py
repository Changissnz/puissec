"""
<Sec> reference map used to compare <Cracker> 
solution candidates to actual solution or 
assign each of those solution candidates a Pr. 
value for their being the solution. 
"""
from sec_mapr import * 

class SRefMap: 

    def __init__(self,opmn,dms,cdms):
        self.opmn = opmn 
        self.dms = dms
        self.cdms = cdms 

        self.ssm = secseq_map
        # vertex label -> vector value
        self.v = {}
        return 

    def load_prism_vertices(self,d):
        for (k,v) in d.items():
            assert k in {'greedy-lone',\
            'greedy-d','greedy-c','greedy-dc',\
            'actual','conn-derived'}
            self.v[k] = v
        return

    def preprocess(self):
        return -1 

    def sec_to_pd_ext_map(self,s):
# def greedy_lone_pd_chain_ext(node,dec,is_dm=True)
        assert s in self.opmn
        sec2pdchain = {}

        for k in self.opmn[s].items():
            greedy_lone_pd_chain_ext(s,k)

        return -1

    def pd_ext_for_sec(self,s): 
        assert s in self.opmn
        return -1 

    # get the Pr. range of a node
    def pr_range_of_node_dec(self,node_idn,dec_idn,
        opmi,dm,decision_chain):
        prmap = dep_weighted_Pr_for_node_dec(node_idn,\
            dec_idn,opmi,dm,decision_chain)
        return 

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

    """
    d := dict, sec idn -> local optima index
    """
    def pr_of_selection(self,d):
        return -1


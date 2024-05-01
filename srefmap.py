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

    # get the Pr. range of a node
    def pr_range_of_node_dec(self,node_idn,dec_idn,
        opmi,dm,decision_chain):
        prmap = dep_weighted_Pr_for_node_dec(node_idn,\
            dec_idn,opmi,dm,decision_chain)
        return 

    """
    d := dict, sec idn -> local optima index
    """
    def pr_of_selection(self,d):
        return -1


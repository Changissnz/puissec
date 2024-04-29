"""
<Sec> reference map used to compare <Cracker> 
solution candidates to actual solution or 
assign each of those solution candidates a Pr. 
value for their being the solution. 
"""
from sec_mapr import * 

class SrefMap: 

    def __init__(self,secseq_map):
        self.ssm = secseq_map
        return 

    def preprocess(self):
        return -1

    """
    d := dict, sec idn -> local optima index
    """
    def pr_of_selection(self,d):
        return -1



    


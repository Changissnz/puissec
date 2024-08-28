from morebs2 import matrix_methods 
import numpy as np 

def any_intersecting_bounds(bounds_seq):
    for i in range(0,len(bounds_seq) - 1):
        for j in range(i+1,len(bounds_seq)):
            b1 = bounds_seq[i] 
            b2 = bounds_seq[i+1]
            q = matrix_methods.intersection_of_bounds(b1,b2)
            if type(q) != type(None):
                return True
    return False 

"""
The Hypothesis structure.

uses vector logs for 3rd-party scores
S with the `suspected_subbounds`, to 
produce an output
P(s in S) --> is s relevant?
"""
class HypStruct:

    def __init__(self,seq_idn:int,targeted_index:int,suspected_subbounds,
        sb_pr=None,hs_vec=None):
        ##assert type(seq_idn) == int and seq_idn >= 0         
        assert type(targeted_index) == int and targeted_index >= 0 
        assert len(suspected_subbounds) > 0 
        assert not any_intersecting_bounds(suspected_subbounds)
        
        self.seq_idn = seq_idn
        self.target_index = targeted_index
        self.suspected_subbounds = suspected_subbounds

        # the probability/weight of each sub-bound
        assert type(sb_pr) == type(None) or \
            matrix_methods.is_vector(sb_pr)

        if type(sb_pr) == type(None):
            sb_pr = np.ones((len(suspected_subbounds),),dtype=float)
            sb_pr = sb_pr * 1.0 / len(suspected_subbounds)
        if type(hs_vec) == type(None):
            hs_vec = np.ones((len(suspected_subbounds),),dtype=int)
            hs_vec = hs_vec * 5

        assert len(sb_pr) == len(suspected_subbounds)
        self.sb_pr = sb_pr
        self.hs_vec= hs_vec 
        return

    def new_HypStruct_by_next_subbound(self):
        if len(self.suspected_subbounds) == 0:
            return None

        sb2 = [self.suspected_subbounds.pop(0)]
        pr2 = np.array([self.sb_pr[0]])
        hsz = np.array([self.hs_vec[0]],dtype=int)
        self.sb_pr = self.sb_pr[1:] 
        self.hs_vec = self.hs_vec[1:]

        hs2 = HypStruct(self.seq_idn,self.target_index,\
            sb2,pr2,hsz)  
        return hs2

    def secdim(self):
        if not len(self.suspected_subbounds): 
            return None 

        q = self.suspected_subbounds[0].shape
        return q[0] 

    def mark_subbound(self,i,v=0.0):
        assert i >= 0 and i < len(self.suspected_subbounds)
        assert type(v) == float 
        self.sb_pr[i] = v
        return

    def most_probable_subbound(self):
        qx = self.most_probable_subbound_i()
        return self.suspected_subbounds[qx[0]]

    def most_probable_subbound_i(self):
        q = list(enumerate(self.sb_pr))
        qi = max(q,key=lambda x:x[1])
        return qi[0]

    def add_subbound(self,sb,sbpr,zero_other_pr=False):
        assert matrix_methods.is_proper_bounds_vector(sb)

        if zero_other_pr:
            self.sb_pr = np.zeros((len(self.suspected_subbounds),),dtype=float)
        
        self.suspected_subbounds.append(sb)
        self.sb_pr = np.append(self.sb_pr,sbpr)

        s = np.sum(self.sb_pr) 
        sx = np.array([measures.zero_div(s_,s,0.0) for s_ in self.sb_pr])
        self.sb_pr = sx
        return 

    def __str__(self):
        q = "sec: " + str(self.seq_idn) + " opt.: " + str(self.target_index)
        q += "\n" + "sub-bounds"
        for s in self.suspected_subbounds:
            q += "\n\t" + str(s) 
        q += "\n"
        q += "\n" + "sub-bound pr."
        q += str(self.sb_pr)
        q += "\n"
        q += "\n" + "hop vec"
        q += str(self.hs_vec)
        q += "\n"
        return q 




"""
the closest int|float multiple of v1
to v2
"""
def closest_single_multiple__v2v(v1,v2):
    bs,_ = CVec__scan_kmult_search(v2,v1,depth=1)
    bsx = max([bs,1.0])
    return int(round(bsx,0.))

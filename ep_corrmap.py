from pr_catgen import * 

################################################### 

"""
Map struct that uses a sequence of 
Python3 <collections.Counter> objects 
to derive probability values of type
* lone probability
* dependent probability
* codependent probability

according to the probability maps
- exact-correlation
- partial-correlation. 

"""
class DerivatorPrMap:

    def __init__(self):

        self.cnt = []
        return

    """
    counts the frequency of elements 
    of two 2-d index vectors alongside axis 0
    """
    def process_index_pair(self,p1,p2):
        c = counter_for_index_2d_op(\
            p1,p2,axis=0)
        
        #self.cnt += c 
        self.cnt.append(c)
        return c

    def fin_count(self,corr_type:str,\
        pred_opt2pr_map:defaultdict):
        assert corr_type in {"e","p"}
        cf = partial_correlation_pr if corr_type == "p" \
            else exact_correlation_pr

        #sz = float('inf') if type(sz) == type(None) \
        #    else sz 
        i = 0
        pr_vec = []
        # lone Pr 
        while len(self.cnt) > 0:# and i < sz:
            cnt = self.cnt.pop(0)
            pr = self.output_Pr(cnt,cf,pred_opt2pr_map)
            pr_vec.append(pr)
            i += 1 

        # dependent and codependent Pr
        return pr_vec

    def output_Pr(self,cnt,cf,pred_opt2pr_map:defaultdict):
        return cf(cnt,pred_opt2pr_map)
        
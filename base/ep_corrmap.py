from .pr_catgen import * 

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
        self.cnt.append(c)
        return c

    """
    - arguments: 
    corr_type := e|p 
    d2_rsz := int, dimension of dist. D2 
    pred_opt2pr_map := 2-d index of D1 -> probability 

    - return: 
    [0] vector, lone Pr, length is number of samples (|`cnt`|)
    [1] defaultdict(float), key is 
    """
    def fin_count(self,corr_type:str,\
        d2_rsz:int,pred_opt2pr_map:defaultdict):
        assert corr_type in {"e","p"}
        cf = partial_correlation_pr if corr_type == "p" \
            else exact_correlation_pr
        i = 0

        # lone Pr 
        pr_vec = []
        exact_corr = {}
        pred_corrmap = defaultdict(float)
        while len(self.cnt) > 0:
            # collect frequency for lone Pr 
            cnt = self.cnt.pop(0)
            pr = self.output_Pr(cnt,cf,pred_opt2pr_map)
            pr_vec.append(pr)

            # collect values for e.c. pr. 
            exact_corr[i] = exact_correlation(cnt)

            # update correlation map for predecessors
            pred_corrmap[i] = map_freq2pr(cnt)
            i += 1 

        """
        # dependency map
        prq = None
        if corr_type == "p":
            prq = partial_correlation_dep_Pr(\
                d2_rsz,exact_corr,\
                pred_corrmap,pred_opt2pr_map)
        else: 
            prq = exact_correlation_dep_Pr(\
                d2_rsz,exact_corr,pred_opt2pr_map)
        """

        return pr_vec,exact_corr#,prq 

    def output_Pr(self,cnt,cf,pred_opt2pr_map:defaultdict):
        return cf(cnt,pred_opt2pr_map)
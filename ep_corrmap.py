from defaults import *

"""
Contains basic functions for calculating 
exact-correlation and partial-correlation
probabilities. 
"""


"""
opt2freq_map := optima index-> frequency of occurence as derivator
                                element.
pred_opt2pr_map := optima index -> Pr that it is the answer. 
"""
def exact_correlation_pr(opt2freq_map,pred_opt2pr_map):

    d = [(k,v) for (k,v) in opt2freq_map.items()]
    d = sorted(d,key=lambda x: x[0])
    d = sorted(d,key=lambda x: x[1])[::-1]
    ans = d[0][0] 

    if ans not in pred_opt2pr_map:
        print("[!] WARNING: null value for exact-correlation function")
        return 0.0
    return pred_opt2pr_map[ans]
    
"""
formula for partial correlation probability is 
"""
def partial_correlation_pr(opt2freq_map,pred_opt2pr_map):
    ##
    """
    print("PARTIAL")
    print("[0]")
    print(opt2freq_map)
    print("[1]")
    print(pred_opt2pr_map)
    """
    ##

    l = list(opt2freq_map.values())
    if len(l) == 0: return 0.0

    s = sum(l)
    d = {(k,pred_opt2pr_map[k]/s) for k in l} 
    d_ = defaultdict(float)

    for (k,v) in d.items():
        v_ = pred_opt2pr_map[k]
        d_[k] = v_ * v
    return sum(list(d_.values())) / len(d_)

def counter_for_index_2d_op(v1,v2,axis=0):
    assert axis in {0,1}

    cnt = Counter()
    for (i,j) in zip(v1,v2):
        cnt[i[axis]] += 1
        cnt[j[axis]] += 1
    return cnt

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
        pred_opt2pr_map:defaultdict,sz=None):
        assert corr_type in {"e","p"}
        cf = partial_correlation_pr if corr_type == "p" \
            else exact_correlation_pr

        sz = float('inf') if type(sz) == type(None) \
            else sz 
        i = 0
        pr_vec = []
        while len(self.cnt) > 0 and i < sz:
            cnt = self.cnt.pop(0)
            pr = self.output_Pr(cnt,cf,pred_opt2pr_map)
            pr_vec.append(pr)
            i += 1 
        return pr_vec

    def output_Pr(self,cnt,cf,pred_opt2pr_map:defaultdict):
        return cf(cnt,pred_opt2pr_map)
        
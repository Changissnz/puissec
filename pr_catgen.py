
from defaults import *
from morebs2 import measures,matrix_methods,aprng_gauge

"""
Contains basic functions for calculating 
exact-correlation and partial-correlation
probabilities. 
"""

def default_dotkey_func():
    f = lambda s1,s2: matrix_methods.vector_to_string(s1 + s2,int)
    return f

def map_freq2pr(m):
    assert type(m) in {defaultdict,dict,Counter}
    s = sum(list(m.values()))
    m_ = defaultdict(float,\
        [(k,measures.zero_div(v,s,0.0)) \
        for (k,v) in m.items()])
    return m_ 


def map_dot_product(d1,d2,x2x2_func=np.multiply,\
    dotkey_func=default_dotkey_func()):
    assert type(d1) in {defaultdict,dict,Counter}
    assert type(d2) in {defaultdict,dict,Counter}

    dx = defaultdict(float)
    for (k,v) in d1.items():
        for (k2,v2) in d2.items():
            k3 = dotkey_func(k,k2) 
            dx[k3] = x2x2_func(v,v2)
    return dx 

################################ lone Pr calculations

def counter_for_index_2d_op(v1,v2,axis=0):
    assert axis in {0,1}

    cnt = Counter()
    for (i,j) in zip(v1,v2):
        cnt[i[axis]] += 1
        cnt[j[axis]] += 1
    return cnt

    ################# exact corr. 
"""
opt2freq_map := optima index-> frequency of occurence as derivator
                                element.
pred_opt2pr_map := optima index -> Pr that it is the answer. 
"""
def exact_correlation_pr(opt2freq_map,pred_opt2pr_map):
    """
    d = [(k,v) for (k,v) in opt2freq_map.items()]
    d = sorted(d,key=lambda x: x[0])
    d = sorted(d,key=lambda x: x[1])[::-1]
    ans = d[0][0] 
    """
    ans = exact_correlation(opt2freq_map)
    if ans not in pred_opt2pr_map:
        print("[!] WARNING: null value for exact-correlation function")
        return 0.0
    return pred_opt2pr_map[ans]

def exact_correlation(freq_map):
    d = [(k,v) for (k,v) in freq_map.items()]
    d = sorted(d,key=lambda x: x[0])
    d = sorted(d,key=lambda x: x[1])[::-1]
    return d[0][0]


    ################# partial corr. 
    
"""
formula for partial correlation probability p is 

l := sum(pr_i * f_j)
p := l / |pr_i|.

The output map <opt2freq_map>

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

    pm = partial_correlation_map(opt2freq_map,pred_opt2pr_map)
    return partial_correlation_pr_v2(pm)
    

def partial_correlation_pr_v2(pm):
    if len(pm) == 0.0: return 0.0 
    return sum(list(pm.values())) / len(pm)

def partial_correlation_map(opt2freq_map,pred_opt2pr_map):
    ##
    l = list(opt2freq_map.values())
    if len(l) == 0: return defaultdict(float)

    s = sum(l)
    d = {(k,pred_opt2pr_map[k]/s) for k in l} 
    d_ = defaultdict(float)

    for (k,v) in d.items():
        v_ = pred_opt2pr_map[k]
        d_[k] = v_ * v
    return d_ 


################################ dependency Pr calculations

"""
calculates the map of

    pair of D2 row indices --> bond score in [0.,1.]

pred_corrmap := index of D2 -> index of D1
pred_opt2pr_map := predecessor map; pair-set index -> bond measure. 
"""
#def exact_correlation_dep_Pr(opt2freq_map,pred_corrmap,pred_opt2pr_map):
def exact_correlation_dep_Pr(d2_rsz,pred_corrmap,pred_opt2pr_map):
    # make the index gen. 
    bis = aprng_gauge.BatchIncrStruct(d2_rsz,is_perm=False,\
        is_reflective=False,subset_size=2)
    
    pr_index = None
    stat = True
    pr_map = defaultdict(float)
    while stat: 
        pr_index = next(bis)
        stat = type(pr_index) != type(None)
        if not stat: continue 

        # new key 
        vstr = matrix_methods.vector_to_string(pr_index,int)

        # old key 
        ##print("PR INDEX: ",pr_index)
        ##print("PRED CORRMAP\n{}\n".format(pred_corrmap))

        pm1 = pred_corrmap[pr_index[0]]
        pm2 = pred_corrmap[pr_index[1]]
        v2 = sorted([pm1,pm2])
        vstr2 = matrix_methods.vector_to_string(v2,int)
        
        pr_val = 0.0 
        if vstr2 in pred_opt2pr_map:
            pr_val =  pred_opt2pr_map[vstr2] 

        pr_map[vstr] = pr_val

    return pr_map

"""
pred_exact_corrmap := index of D2 -> index of D1
pred_corrmap := index of D2 -> index of D1 -> bond strength 
pred_opt2pr_map := 2-d index of D1 -> probability 
"""
def partial_correlation_dep_Pr__process_func(index,pred_exact_corrmap,\
    pred_corrmap,pred_opt2pr_map): 
    assert len(index) == 2

    i0 = pred_exact_corrmap[index[0]]
    i1 = pred_exact_corrmap[index[1]]

    x0 = pred_corrmap[i0]
    x1 = pred_corrmap[i1] 
    map_dot = map_dot_product(x0,x1)


    # 
    md2 = {}
    i = 0 
    for (k,v) in map_dot.items():
        md2[i] = v * pred_opt2pr_map[k] 
        i += 1

    if len(md2) == 0: return 0.0
    q = sum(list(md2.values())) / len(md2)
    return round(q,5)

"""
d2_rsz := number of rows belonging to dist. D2 
pred_exact_corrmap := index of D2 -> index of D1
pred_corrmap := index of D2 -> index of D1 -> bond strength 
pred_opt2pr_map := 2-d index of D1 -> probability 

return:
- 
"""
def partial_correlation_dep_Pr(d2_rsz,pred_exact_corrmap,\
    pred_corrmap,pred_opt2pr_map):

    # make the index gen. 
    bis = aprng_gauge.BatchIncrStruct(d2_rsz,is_perm=True,\
        is_reflective=True,subset_size=2)
    
    ##partial_correlation_pr_v2(pm)
    stat = True
    pr_map = defaultdict(float)
    l = 0
    while stat:

        pr_index = next(bis)
        stat = type(pr_index) != type(None)
        if not stat: continue 

        pfc = partial_correlation_dep_Pr__process_func(pr_index,pred_exact_corrmap,\
            pred_corrmap,pred_opt2pr_map)

        key = matrix_methods.vector_to_string(pr_index,int)
        pr_map[key] = pfc 
    return pr_map 
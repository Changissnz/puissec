
#from defaults import *
from morebs2 import measures,numerical_generator
from imod import * 

"""
Contains basic functions for calculating 
exact-correlation and partial-correlation
probabilities. 
"""

def generate_2d_match_pr(dim):
    assert len(dim) == 2
    assert min(dim) > 0

    bis = aprng_gauge.BatchIncrStruct(max(dim),\
        True,True,2)
    d_ = np.array([[0,dim[0]],\
            [0,dim[1]]])
    i2dm = Index2DMod(d_,bis)

    stat = True
    pr_map = defaultdict(float)
    while stat: 
        ind = next(i2dm)
        stat = not type(ind) == type(None)
        if not stat: continue

        sind = matrix_methods.vector_to_string(ind,int)
        pr_map[sind] = random.random() 
    return pr_map

def generate_pr_dist_D1_to_D2(d1_sz,d2_sz):
    prdist = defaultdict(None)
    for i in range(d2_sz):
        prdist[i] = defaultdict(float)
        for j in range(d1_sz):
            prdist[i][j] = random.random()

        s = sum(list(prdist[i].values()))
        for j in range(d1_sz):
            prdist[i][j] = measures.zero_div(prdist[i][j],s,0.0)
            prdist[i][j] = round(prdist[i][j],5)
    return prdist

def replace_map_key(m,f):
    m_ = defaultdict(None)
    for (k,v) in m.items():
        m_[f(k)] = v 
    return m_ 

def split_map(m,ratio=0.5):
    m_ = defaultdict(None)
    m2_ = defaultdict(None)

    for (k,v) in m.items():
        if random.random() >= ratio:
            m2_[k] = v
        else: 
            m_[k] = v 
    return m_,m2_ 

############################################
############# map dot products 

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
            print("K: ", k,"\tK2: ", k2)
            k3 = dotkey_func([k],[k2]) 
            dx[k3] = x2x2_func(v,v2)
    return dx 

def map_dot_product_(d1,d2):

    dx = defaultdict(float)
    for (k,v) in d1.items():
        dx[k] = v * d2[k] 
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

<opt2freq_map>
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

"""
opt2freq_map := optima index-> frequency of occurence as derivator
                                element.
pred_opt2pr_map := optima index -> Pr that it is the answer. 
"""
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

- return: 
defaultdict, index pair of D2 -> bond strength 
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

        pr_map[vstr] = round(pr_val,5)

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

    ##
    """
    print("PRED CORRMAP")
    print(pred_corrmap)
    print()
    print("PRED OPT2PR MAP")
    print(pred_opt2pr_map)
    print()
    print("----")
    """
    ##

    x0 = pred_corrmap[i0]
    x1 = pred_corrmap[i1] 

    ##
    """
    print("mapping elements")
    print(x0)
    print(x1)
    print()
    """
    ##

    map_dot = map_dot_product_(x0,x1)

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
- defaultdict, index pair of D2 -> bond strength 
"""
def partial_correlation_dep_Pr(d2_rsz,pred_exact_corrmap,\
    pred_corrmap,pred_opt2pr_map):

    # make the index gen. 
    bis = aprng_gauge.BatchIncrStruct(d2_rsz,is_perm=True,\
        is_reflective=True,subset_size=2)    
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

##################################################
########### bounds generator

"""
- arguments: 
superbound := bounds vector, the bound that all 
    output bounds have to fall within.
spacing_ratio_range := the permitted range of 
    ratios (w.r.t. superbound) allowed for the distance 
    between any two output bounds.
outlier_pr_ratio := the permitted range of variation,
    given the remaining allocatable space, 
    allowed for each bounds.
num_bounds := int,the number of bounds to be generated.

*Examples*
(1)
    spacing_ratio_range = (0.,0.) 
    - no spacing; all allocated space is contiguous (connected)
(2) 
    outlier_pr_ratio = (0.,0.1)

    The default vector length for each bound is 
        L_d = (superbound[1] - superbound[0]) / num_bounds 

    The range of values that a Python random generator 
    draws from is [L_d + 0.,L_d + 0.1 * 
                    (||superbound|| * L_d)]

(3) 
    outlier_pr_ratio = (1.0,1.0) 

    The range of values that a Python random generator
    draws from is 

        [B[1], superbound[1]]; B the previous bound.
"""
def generate_bounds_vector_sequence(superbound,\
    spacing_ratio_range,outlier_pr_ratio,num_bounds):
    assert matrix_methods.is_proper_bounds_vector(superbound)
    assert num_bounds > 0 and type(num_bounds) == int

    outlier_pr_calibrated = [1/num_bounds,1/num_bounds]
    outlier_pr_calibrated[0] = max([0.,\
        outlier_pr_calibrated[0] - outlier_pr_ratio])
    outlier_pr_calibrated[1] = min([1.,\
        outlier_pr_calibrated[1] + outlier_pr_ratio])

    sb_diff = superbound[:,1] - superbound[:,0] 
    bseq = []
    ref_point = deepcopy(superbound[:,0])
    next_point = None
    while num_bounds > 0: 

        # update the next point 
        next_point = generate_point_for_bvecseq(superbound,\
        sb_diff,ref_point,outlier_pr_calibrated)

        bd = deepcopy(np.array([ref_point,next_point]).T)

        # update the ref point 
        ref_point = generate_point_for_bvecseq(superbound,\
        sb_diff,next_point,spacing_ratio_range)

        bseq.append(bd)
        num_bounds -= 1
    return bseq 

def generate_point_for_bvecseq(superbound,\
    superbound_diff,ref_point,ratio_range):
    assert len(ratio_range) == 2 
    assert ratio_range[0] <= ratio_range[1]
    assert ratio_range[0] >= 0. and ratio_range[1] <= 1.0 

    # get the difference b/t the ref_point 
    # and the superbound max
    r_diff = superbound[:,1] - ref_point
    assert np.min(r_diff) >= 0.0 


    # calculate a ratio vecor
    rvec = numerical_generator.generate_uniform_sequence_in_bounds\
        (superbound.shape[0], np.array([ratio_range]))
    
    adder = r_diff * rvec
    return np.round(ref_point + adder,5)
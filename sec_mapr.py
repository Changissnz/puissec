from ep_corrmap import * 
from bloominhurst import *


default_rnd_boolean_index_splitter = lambda x: True if random.randrange(0,2) else False 

## TODO: test this section 
##################### methods for reading 
##################### (co?)-dependency maps
##################### belonging to <Sec> 

"""
return:
- number of components, (minimum length,maximum length)
"""
def metrics_on_depmap(dm,is_dep):

    cs = connected_subsets_of_depmap(dm)
    if is_dep:
        cs_ = []
        for c in cs: 
            cs_.append(sort_depmap_component_by_dependency(dm,c)) 
        lcs = [len(c) for c in cs_]  
    else: 
        lcs = [len(c) for c in cs] 
    return len(cs), [min(lcs),max(lcs)] 


"""
return: 
- sequence, set of connected <Sec> instances 
"""
def connected_subsets_of_depmap(dm): 
    
    finished = [] 
    s = set(dm.keys())
    s_ = []

    def index_in_soln(n):
        for (j,s2) in enumerate(s_):
            if s2 == n: return j 
        return -1 

    for r in s: 
        px = parse_dconn(r) 

        i1 = index_in_soln(px[0])
        i2 = index_in_soln(px[1]) 

        if i1 == -1 and i2 == -1:
            s_.append(set([px[0],px[1]]))
        elif i1 != -1 and i2 != -1:
            q = sorted([i1,i2]) 
            t1 = s_.pop(q[1])
            t0 = s_.pop(q[0]) 
            s_.append(t0 | t1) 
        elif i1 == -1 and i2 != -1:
            s_[i2] = s_[i2] | {px[0]} 
        elif i1 != -1 and i2 == -1:
            s_[i1] = s_[i1] | {px[1]} 
        else:
            assert False
    return s_

"""
return: 
- sequence, each element is a set, right-most 
  element is most dependent
"""
def sort_depmap_component_by_dependency(dm,component:set):
    
    cs = deepcopy(component)
    s_ = []
    def index_in_soln(n):
        for (j,s2) in enumerate(s_):
            if s2 == n: return j 
        return -1 

    for k in dm.keys():
        px = parse_dconn(k)

        i1 = index_in_soln(px[0])
        i2 = index_in_soln(px[1]) 

        if i1 == -1 and i2 == -1:
            q1 = {px[0]}
            q2 = {px[1]}
            s_.extend([q2,q1])

        if i1 != -1 and i2 != -1:
            assert i2 < i1

        if i1 == -1 and i2 != -1:
            s_.append({px[0]}) 
        
        if i1 != -1 and i2 == -1: 
            s_.insert(0,{px[1]}) 
    return s_ 

"""
return: 
- number of connections, 
- map, sec index -> set of its optima indices 
"""
def metrics_on_node_in_depmap(dm,n):
    assert type(dm) == defaultdict

    count = 0
    # sec index -> set of its optima indices 
    other_sec = defaultdict(set)
    for (k,v) in dm.items():
        px = parse_dconn(k)
        if px[0] != n: continue
        count += 1
        other_sec[px[1]] = other_sec[px[1]] | {px[2]} 
    return count,other_sec 

####################### generator functions
####################### for <Sec> vars 

def generate_default_szmult_mod(m,drange=[2,11]):
    assert len(drange) == 2 and drange[0] <= drange[1]

    def f(x):
        q = int(round(x * m) % (drange[1] - drange[0]))
        return q + drange[0] 
    return f 

def default_AltBaseFunc_for_IsoRing():
    q = generate_efunc_type_q(1,1)#,1,1)
    q2 = generate_efunc_type_q(0.5,1.0)#,0.5,1.0)
    mfs = deepcopy(DEFAULT_PAIRWISE_VEC_FUNCS)
    mfs = mfs + [q,q2] 
    return AltBaseFunc(mfs,random)

from ep_corrmap import * 
from bloominhurst import *


default_rnd_boolean_index_splitter = lambda x: True if random.randrange(0,2) else False 

## TODO: test this section 
##################### methods for reading 
##################### (co?)-dependency maps
##################### belonging to <Sec> 

"""
return: 
- number of connections, 
- set, of its optima indices 
"""
def metrics_on_node_in_depmap(dm,n):
    assert type(dm) == defaultdict

    count = 0
    # sec index -> set of its optima indices 
    other_sec = set()
    for (k,v) in dm.items():
        px = parse_dconn(k)
        print("parsed results: ",px)
        if px[1] != n: continue
        count += 1
        other_sec = other_sec | {px[2]}
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

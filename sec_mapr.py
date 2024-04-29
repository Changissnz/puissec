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

# TODO: test this
"""
sec2dm := dict, sec. idn -> dependency idn -> Pr
"""
def connected_subsets_of_codepmap(sec2dm):
    s = set(sec2dm.keys())
    s_ = [{i} for i in s] 

    def index_in_soln(n):
        for (j,s2) in enumerate(s_):
            if n in s2: return j 
        return -1 

    dun = set() 
    for s2 in s:
        # get all codep. Secs
        q = sec2dm[s2]
        ql = list(q.keys())
        qls = set([parse_dconn(ql_)[1] for ql_ in ql])

        # skip the ones already iterated;
        qls = qls - dun 
        for qls_ in qls:
            i1 = index_in_soln(s2)
            i2 = index_in_soln(qls_)
            if i1 == i2: continue

            x = s_[i1] | s_[i2]
            s_[i1] = x
            s_.pop(i2) 

        dun |= {s2} 
    return s_ 

# TODO: test this
"""
sec2dm := dict, sec. idn -> dependency idn -> Pr
sec_id := 
"""
def depchain_for_Sec(sec2dm,sec_id):

    q = [sec_id]
    chain = [{sec_id}]

    def in_chain(idn):
        for x in chain: 
            if idn in x: return True
        return False 

    while len(q) > 0:
        si = q.pop(0)

        d = sec2dm[si]
        ql = list(d.keys())
        qls = set([parse_dconn(ql_)[1] for ql_ in ql])

        # check for circular dependencies
        for q_ in qls: 
            if in_chain(q_):
                return None 
        if len(qls) > 0: 
            chain.append(qls)
        q.extend(list(qls)) 

    return chain 

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

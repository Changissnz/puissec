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
def metrics_on_node_in_depmap(dm,n,full_output=False):
    assert type(dm) == defaultdict

    count = 0
    # optima indices
    other_sec = set()
    this_sec = set()
    for (k,v) in dm.items():
        px = parse_dconn(k)
        if px[1] != n: continue
        count += 1
        other_sec = other_sec | {px[2]}

        if full_output:
            this_sec = this_sec | {px[0]}

    if not full_output:
        return count,other_sec
    return count,other_sec,this_sec 

"""
dm := dict, dep. map for dependent node.
nondep_node := int, identifier for node depended on.
os1 := int, number of optima for dependent 
os2 := int, number of optima for non-dependent

return:
- available optima for dependent node
- available optima for non-dependent node 
"""
def available_for_dep(dm,nondep_node,os1,os2):
    # number of connections, 
    # local optima for nondep_node, 
    # local optima for dep node 
    count,other_sec,this_sec = metrics_on_node_in_depmap(dm,nondep_node,True)
    
    q1 = sorted(list(set(range(os1)) - this_sec))
    q2 = sorted(list(set(range(os2)) - other_sec))
    return q1,q2

"""
dm := dict, dep. map for a <Sec> instance. 

return: 
- list, elements are 
    (other <Sec> idn,index of other <Sec> optimum, bond strength).
"""
def filter_optimum_conn_in_depmap(dm,o): 
    q = []
    for (k,v) in dm.items():
        k_ = parse_dconn(k) 
        if k_[0] != o: continue
        s = (k_[1],k_[2],v)
        q.append(s) 
    return q 

"""
return: 
- number of conn.,cumulative bond strength
"""
def metrics_on_optimum_in_depmap(dm,o):
    fo = filter_optimum_conn_in_depmap(dm,o)
    n = len(fo) 
    bs = sum(list(fo.values()))
    return n,bs

# TODO: test. 
"""
possible-decision chain by extremum 
function `extf` for optimum `o` in 
(co?)-dep. map `dm`.

return:
- dict, sec idn -> Set(decision with extreme Pr. value)
"""
def extdec_dmap_set(dm,o,extf):
    assert extf in {min,max}

    oconn = filter_optimum_conn_in_depmap(dm,o)

    # get the max
    def sort_next():
        if len(oconn) == 0: 
            return
        x = oconn.pop(0)
        max0 = x[2] 
        max_set = set([x[1]])
        ref_sec = x[0] 
        i = 0 
        while i < len(oconn):
            x2 = oconn[0]

            # case: not relevant
            if x2[0] != ref_sec:
                i += 1
                continue

            x2 = oconn.pop(0)
            q = [max0,x2[2]]
            
            # case: equals
            if abs(q[0] - q[1]) < 10 ** -4:
                max_set = max_set | {x2[1]} 
                continue

            # case: extremum
            mx = extf(q)
            if x2[2] == mx:
                max_set = set([x2[1]])
                max0 = x2[2]
        return ref_sec,max_set 

    d = defaultdict(set)
    while len(oconn) > 0:
        q = sort_next()
        if type(q) == type(None): 
            continue
        d[q[0]] = q[1] 

    return d

"""
pdc := dict, sec idn. -> set(<opt.>)
dc := dict, sec idn. -> (set(<opt.>)|opt idn.)
fullkey_req := bool, set to True, method asserts
            both `pdc.
            
return:
- is dc a subset of pdc? 
"""
def is_in_pd_chain(pdc,dc,fullkey_req=False):

    if fullkey_req:
        stat = set(dc.keys()) == set(pdc.keys()) 
        if not stat: 
            return False 

    for k,v in dc.items():
        assert type(v) in {int,np.int32,set}
        v_ = None
        if type(v) != set:
            v_ = {v}
        else:
            v_ = v 

        if not v_.issubset(pdc[k]):
            return False
    return True 

"""
calculates the connected subsets of <Sec>
instances based on non-null co-dependency 
values using `sec2dm`.

sec2dm := dict, sec. idn -> dependency idn -> Pr

return: 
- list, each element a Set(codependent) 
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

"""
sec2dm := dict, sec. idn -> dependency idn -> Pr
sec_id := int, identifier for <Sec>

return:
- list, each element a Set(dependency)
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

####################### range-of-interpretation
####################### formulations

# TODO: test this 
'''
(co?)-dependent weighted Pr for node

n := int, identifier for node 
ndec := int, identifier for decision on `n`. 
opm := dict, int of optimum -> float, Pr value. 
dm := dict, (co?)-dependency map pertaining to `n`. 
decision_chain := list, (sec index, index of selected optimum)

'''
def dep_weighted_Pr_for_node_dec(n,ndec,opm,dm,decision_chain): 
    # gather the relevant probability values
    parsed_dc = []
    prvs = []
    for d in decision_chain:
        key = str(ndec) + ","\
            + d[0] + "." + d[1]
        prvs.append(dm[key])
    assert ndec in opm

    # add weighted additions to Pr(ndec) 
    v = opm[ndec]
    v2 = v

    for x in prvs: 
        v2 = v2 + (v * x) 

    # normalize the value 
    qk = set(opm.keys()) - {ndec}
    s = 0.0
    for qk_ in qk:
        s += opm[qk_]
    s += v2

    opm_ = defaultdict(float) 
    for (k,v) in opm.items():
        opm_[k] = measures.zero_div(v,s,0.0)
    return opm_

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

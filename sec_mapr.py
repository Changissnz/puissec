from ep_corrmap import * 
from bloominhurst import *
from morebs2 import search_space_iterator

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
        qls_ = [] 
        for q_ in qls: 
            if in_chain(q_):
                continue
            else:
                qls_.append(q_) 
        qls = set(qls_)

        if len(qls) > 0: 
            chain.append(qls)
        q.extend(list(qls)) 

    return chain 

########################### Order-of-Cracking permutation
########################### PR. map noise-adder 

# sec idn -> sec dim. -> opt. idn -> Pr. 
def noise_on_opt_pr_map_SEEDTYPE_PythonANDNumpy(opm,rnd_struct): 

    def apply_noise_on_idnANDdim(k,k2):
        # retrieve dict
        vx = opm[k][k2]
        qr = [(v1,v2) for (v1,v2) in vx.items()]
        q2 = [v[1] for v in qr]
        if len(q2) == 0:
            return 

        # add noise
        m = np.zeros((len(q2),2),dtype="float")
        m[:,1] = 1.0 
        no = Noiseano(m,rnd_struct=rnd_struct)
        q2 = no.noisha(np.array(q2))


        qr = [(v[0],q2[i]) for (i,v) in enumerate(qr)]

        # normalize
        s = sum([qr_[0] for qr_ in qr])

        if s == 0.0:
            s = 1.0
        qr = [(v[0],v[1]/s) for v in qr] 

        # add back to dict
        opm[k][k2] = dict(qr)

    for k,v in opm.items():
        for k2 in v.keys():
            apply_noise_on_idnANDdim(k,k2)
    return opm 

"""
***
specialized for use with output from 
`depchain_for_Sec`,`connected_subsets_of_codepmap`
***

- list, each element a Set(codependent) 
"""
def permute_setseq(setseq,rnd_struct,num_swaps:int):
    assert type(setseq) == list
    for s in setseq: assert type(s) == set and len(s) > 0 
    assert num_swaps >= 0 and type(num_swaps) == int
    actual_swaps = 0

    """
    exclude := None|(i'th index in `setseq`)

    return:
    - (i'th set of `setseq`, element in `setseq`[i])
    """
    def choose_with_exclude(exclude=None): 
        # choose one i'th index. 
        l = set([i for i in range(len(setseq))]) 
        if type(exclude) != type(None):
            l = l - set([exclude])

        if len(l) == 0: 
            return None,None
        i = rnd_struct.randrange(0,len(l)) 
        l = list(l)
        i = l[i]

        # choose a j'th index in the i'th set (ordered)
        jx = list(setseq[i])
        j = rnd_struct.randrange(0,len(jx))
        return i, jx[j]

    def swap_one():
        i1,j1 = choose_with_exclude()
        i2,j2 = choose_with_exclude(i1)

        if type(i1) == type(None):
            return False

        if type(i2) == type(None):
            return False

        setseq[i1] = setseq[i1] - {j1}
        setseq[i1] = setseq[i1] | {j2} 
        setseq[i2] = setseq[i2] - {j2}
        setseq[i2] = setseq[i2] | {j1} 
        return True
    
    while num_swaps > 0:
        stat = swap_one()
        num_swaps -= 1 
        if stat: actual_swaps += 1 
    return setseq,actual_swaps



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
            + str(d[0]) + "." + str(d[1])
        if key in dm:
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
    opm_[ndec] = measures.zero_div(v2,s,0.0)
    return opm_


"""
used in frequency-counting for <SRefMap>; updates 
`dcnt` w/ `pd`. 

arguments:
- dcnt := defaultdict, sec. idn. -> Counter(opt.index -> frequency)
- pd := dict, possible-decision map

return: 
- defaultdict, sec. idn. -> Counter(opt.index -> frequency)
"""
def update_SRefMap_counter(dcnt:defaultdict,pd): 
    assert type(dcnt) == defaultdict
    assert type(pd) in {defaultdict,dict} 
    for (k,v) in pd.items():
        x = dcnt[k]
        for v_ in v:
            x[v_] += 1 
    return dcnt

"""
iterator for possible-decision maps, a 
dictionary structure w/ 

    key := sec idn.
    value := set of opt. indices 
"""
class PDMapIter: 

    def __init__(self,pd):
        assert type(pd) in {dict,defaultdict}
        for v in pd.values(): assert type(v) in {list,set}
        
        self.pd = defaultdict(list)
        for k,v in pd.items(): 
            self.pd[k] = list(v) 

        self.set_iter()
        self.sz = 0 
        return 

    def set_iter(self):

        ls = sorted(list(self.pd.keys())) 
        u = [] 
        for ls_ in ls:
            lx = len(self.pd[ls_])
            u.append((ls_,lx)) 
        u = np.array(u)
        v = deepcopy(u[:,1])
        bounds = np.zeros((len(u),2))
        bounds[:,1] = v 
        start_point = np.zeros((len(v),))
        column_order = [i for i in range(len(v))][::-1]
        ssihop = deepcopy(v)
        cycleOn = False
        self.ssi_ref = u 
        self.ssi = search_space_iterator.SearchSpaceIterator(bounds,\
            start_point,column_order,SSIHop=ssihop,cycleOn=cycleOn)

    def reached_end(self):
        return self.ssi.reached_end()

    def __next__(self): 

        if self.reached_end():
            return None 

        # fetch the index 
        i = next(self.ssi) 

        # fetch the element 
        pdx = defaultdict(None)

        for (j,s) in enumerate(self.ssi_ref):
            pdx[s[0]] = self.pd[int(s[0])][int(i[j])]
        self.sz += 1 
        # 
        return pdx 

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

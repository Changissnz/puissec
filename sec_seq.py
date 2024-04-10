from defaults import *
from collections import defaultdict 
from bloominhurst import *

def generate_default_szmult_mod(m,drange=[2,11]):
    assert len(drange) == 2 and drange[0] <= drange[1]

    def f(x):
        q = int(round(x * m) % (drange[1] - drange[0]))
        return q + drange[0] 
    return f 

"""
used to generate points by Python standard
`random` library. 
"""
def generate_pointgen_stdrand(bounds,num_points,rnd_struct):
    assert matrix_methods.is_proper_bounds_vector(bounds)
    
    ps = []
    for i in range(num_points):
        rx = rnd_struct.random(bounds.shape[0])
        p = matrix_methods.point_on_bounds_by_ratio_vector(bounds,rx)
        ps.append(p)
    return np.array(ps)

def default_AltBaseFunc_for_IsoRing():
    q = generate_efunc_type_q(1,1,1,1)
    q2 = generate_efunc_type_q(0.5,1.0,0.5,1.0)
    mfs = deepcopy(DEFAULT_PAIRWISE_VEC_FUNCS)
    mfs = mfs + [q,q2] 
    return AltBaseFunc(mfs,random)

############################################


class OptimaBloomFuncSecRep:

    def __init__(self,obf,sz_map,bloom_sz_limit:int=1000):
        assert type(obf) == OptimaBloomFunc
        assert type(bloom_sz_limit) == int and bloom_sz_limit > 0

        self.obf = obf
        self.sz_map = sz_map
        self.bloom_sz_limit = bloom_sz_limit

        self.sz = None
        self.bpoints = {} 
        self.bpoints[self.sz] = deepcopy(\
            self.obf.oseeds)

        self.set_sz()
        self.bpoints = {}

    def set_sz(self):
        # get dim of OptimaBloomFunc
        self.sz = self.obf.oseeds.shape[1]
        self.sz = self.sz_map(self.sz)
        print("SZ IS: {}".format(self.sz))
        return self.sz

    def iso_appear(self,i):
        assert i in self.bpoints
        return deepcopy(self.bpoints[i])

    def __next__(self):
        stat = self.stat_at_i(self.sz) 

        # case: reset the sz
        if not stat: 
            self.set_sz()
            return None

        self.obf.d = self.sz


        d = next(self.obf)
        if type(d) == type(None):
            return None

        if self.sz not in self.bpoints:
            self.bpoints[self.sz] = np.array([d])
        else:
            v = self.bpoints[self.sz]
            self.bpoints = np.vstack((v,d))

        return d

    def stat_at_i(self,i):
        sample_sz = self.sample_size_at_i(i)
        return sample_sz < self.bloom_sz_limit

    def sample_size_at_i(self,i):
        if i not in self.bpoints: return 0
        s = self.bpoints[i].shape 
        return s[0] * s[1] 

class Sec:

    def __init__(self,sequence,singleton_range,\
        optima_pr_map,dep_map,codep_map):

        assert matrix_methods.is_vector(sequence)
        assert np.array(singleton_range).ndim == 1 and \
            singleton_range[0] <= singleton_range[1] 

        assert type(optima_pr_map) == type(dep_map) and \
            type(dep_map) == type(codep_map)
        assert type(optima_pr_map) == defaultdict

        self.seq = sequence
        self.singleton_range = singleton_range
        self.opm = optima_pr_map
        self.dm = dep_map
        self.cdm = codep_map

        self.declare_obf()

    def declare_obf(self):

        abf = default_AltBaseFunc_for_IsoRing()
        bloom_func = abf.load_AltBaseFunc_function(abf)

        # generate the size-delta function for
        # the <IsoRingBloomer>
        multiplier = 3.13
        szm = generate_default_szmult_mod(multiplier)

        # initialize the <OptimaBloomFunc>
        optima_points = [matrix_methods.string_to_vector(v,float) for v in \
                    self.opm.keys()]
        optima_points = np.array(optima_points)
        bounds = np.array([deepcopy(self.singleton_range) for _ in \
            range(len(self.seq))])

        obf = OptimaBloomFunc(optima_points,bounds,\
            None,bloom_func,optima_points.shape[0],8)
        self.obfsr = OptimaBloomFuncSecRep(obf,szm,DEFAULT_OBF_SECREP_BLOOM_SZ_LIMIT)
        return

    """
    return: 
    - bloom value, [index for element 1, index for element 2]
    """
    def __next__(self):
        q = next(self.obfsr)

        if type(q) == type(None):
            return None,None

        q2 = np.array([self.obfsr.obf.prev_i1,self.obfsr.obf.prev_i2])
        return q, q2

    @staticmethod
    def generate_Sec(seq,rnd_struct):
        return -1

class SecSeq:

    def __init__(self,sequence,reference_map=None): 

        for s in sequence: assert type(s) == Sec
        assert len(reference_map) == len(sequence)
        assert type(reference_map) == defaultdict

        lm = list(reference_map.values())
        minny,maxxy = min(lm),max(lm)
        assert minny >= 0 and maxxy < len(sequence)
        self.sequence = sequence
        self.reference_map = reference_map
        return
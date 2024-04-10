from defaults import *
from collections import defaultdict 
from bloominhurst import *

def generate_default_szmult_mod(m,drange):
    assert len(drange) == 2 and drange[0] <= drange[1]

    def f(x):
        q = (x * m) % (drange[1] - drange[0])
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

class Sec:

    def __init__(self,sequence,singleton_range=DEFAULT_SINGLETON_RANGE,\
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

    def declare_obf():

        abf = default_AltBaseFunc_for_IsoRing()
        bloom_func = abf.load_AltBaseFunc_function(abf)

        # generate the size-delta function for
        # the <IsoRingBloomer>
        drange = deepcopy(self.singleton_range)
        multiplier = 3.13
        szm = generate_default_szmult_mod(multiplier,drange)

        # initialize the <OptimaBloomFunc>
        optima_points = [matrix_methods.vector_to_string(v,float) for v in \
                    self.opm.values()]
        optima_points = np.array(optima_points)
        bounds = np.array([deepcopy(self.singleton_range) for _ in \
            range(len(self.sequence))])

        obf = OptimaBloomFunc(optima_points,bounds,\
            None,bloom_func,optima_points.shape[0],8)
        self.obf = obf
        return

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
"""
for use with single elements of a secret
"""
from sec_seq import * 
from defaults import * 
from obj_func import * 
from bloominhurst import *

def generate_default_szmult_mod(m,drange):
    assert len(drange) == 2 and drange[0] <= drange[1]

    def f(x):
        q = (x * m) % (drange[1] - drange[0])
        return q + drange[0] 
    return f 

class IsoRingBloomer:

    def __init__(self,obf,sz_map):
        assert type(obf) == OptimaBloomFunc
        self.obf = obf
        self.sz_map = sz_map
        self.sz = None
        self.set_sz()
        self.bpoints = []
        return -1 

    def stat(self):
        return self.obj.finished_stat

    def set_sz(self):
        # get dim of OptimaBloomFunc
        self.sz = self.obf.shape[1]
        self.sz = self.sz_map(self.sz)
        return self.sz

    def bloom_one(self):
        if self.stat(): return 

        self.obf.d = self.sz
        d = next(self.obf)
        self.bpoints.append(deepcopy(d))
        return d 

    def reset_obf(self):
        self.bpoints = np.array(self.bpoints)

        return -1 

class IsoRingMem:

    def __init__(self,bounds,optima_points,actual_optima):
        self.bounds = bounds
        self.optima_points = optima_points
        self.actual_optima = actual_optima

##############################################

def default_AltBaseFunc_for_IsoRing():
    q = generate_efunc_type_q(1,1,1,1)
    q2 = generate_efunc_type_q(0.5,1.0,0.5,1.0)
    mfs = deepcopy(DEFAULT_PAIRWISE_VEC_FUNCS)
    mfs = mfs + [q,q2] 
    return AltBaseFunc(mfs,random)
    return

"""
See arguments for <SecSeq> 
"""
class IsoRing:

    def __init__(self,bounds,optima_points,\
        actual_optima,objf):
        assert matrix_methods.is_proper_bounds_vector(bounds)
        for p in optima_points:
            assert matrix_methods.point_in_bounds(bounds,p)
        assert matrix_methods.point_in_bounds(bounds,\
            actual_optima)
        assert type(objf) == ObjFunc

        self.bounds = bounds
        self.optima_points = optima_points
        self.actual_optima = actual_optima
        self.irm = IsoRingMem(bounds,optima_points,\
            actual_optima)

        self.objf = objf
        self.irb = None
        self.set_irb()
        return

    """
    outputs a score for each of the optima points
    """
    def register_attempt(self,p):
        q = np.zeros((self.optima_points.shape[0],))
        rx = [self.objf.output(p_,p) for p_ in self.optima_points]
        return np.array(rx)
        
    def set_irb(self):
        # set the default <AltBaseFunc> as the
        # bloom function. 
        abf = default_AltBaseFunc_for_IsoRing()
        bloom_func = abf.load_AltBaseFunc_function(abf)

        # generate the size-delta function for
        # the <IsoRingBloomer>
        drange = deepcopy(DEFAULT_SINGLETON_RANGE)
        multiplier = 3.13
        szm = generate_default_szmult_mod(multiplier,drange)

        # initialize the <OptimaBloomFunc>
        obf = OptimaBloomFunc(deepcopy(self.optima_points),\
            deepcopy(self.bounds),selector_func=None,\
            bloom_func,d=self.bounds.shape[0],\
            split_sz=8,\
            splitsz_delta=DEFAULT_OPTIMA_BLOOM_SZ_DELTA_FUNC)

        # initialize the IsoRingBloomer
        self.irb = IsoRingBloomer(obf,szm)
        return

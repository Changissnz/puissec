from ep_corrmap import * 
#from collections import defaultdict 
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

    def __init__(self,obf,sz_map,immutable_dim:int,\
        corr_type="e"):#,bloom_sz_limit:int=1000):
        assert type(obf) == OptimaBloomFunc
        assert type(immutable_dim) == int and immutable_dim >= 0
        #assert type(bloom_sz_limit) == int and bloom_sz_limit > 0

        self.obf = obf
        self.sz_map = sz_map

        self.tdim = [immutable_dim]
        self.corr_type = corr_type

        self.sz = self.obf.oseeds.shape[1]
        self.bpoints = {} 
        self.bpoints[self.sz] = deepcopy(\
            self.obf.oseeds)
        
        # dimension d -> 
        # [0] reference dimension of optima 
        #     points used to generate d-points
        # [1] vector of greatest similarity. 
        self.bpoint_dim_ref = {}

        # terminate at dim `sz` 
        self.fstat = False
        self.tstat = False
        self.set_sz()

        self.set_dpm() 

    def set_dpm(self):
        self.dpm = DerivatorPrMap() 

    def set_sz(self,attempts=10):
        # get dim of OptimaBloomFunc
        self.sz = self.obf.oseeds.shape[1]

        self.sz = self.sz_map(self.sz)
        if self.sz in self.tdim:
            if attempts > 0:
                return self.set_sz(attempts - 1)
            return None 

        ##print("SZ IS: {}".format(self.sz))
        return self.sz

    def iso_appear(self,i):
        assert i in self.bpoints
        return deepcopy(self.bpoints[i])

    def __next__(self):
        self.obf.d = self.sz
        if self.sz in self.tdim:
            return None

        d = next(self.obf)
        if type(d) == type(None):
            ##print("[!!] TSTAT")
            self.tstat = True
            self.tdim.append(self.sz)
            return None

        ##print("D: ", d)
        if self.sz not in self.bpoints:
            ##print("D: ",d)
            self.bpoints[self.sz] = np.array([d])
        else:
            v = self.bpoints[self.sz]
            self.bpoints[self.sz] = np.vstack((v,d))

        self.assign_next_value_to_pred(d)
        ##
        """
        print("BPOINTS")
        print(self.bpoints)
        print("------------")
        """
        ##
        return d

    def assign_next_value_to_pred(self,next_val):
        px1 = self.obf.prev_i1
        px2 = self.obf.prev_i2
        self.dpm.process_index_pair(px1,px2)
        return

    def finalize_dcount(self,pred_opt2pr_map:defaultdict):
        # get the dimension of predecessor
        ##
        """
        print("CHECK HERE")
        print(pred_opt2pr_map)
        """
        pr_vec = self.dpm.fin_count(self.corr_type,pred_opt2pr_map)
        self.set_dpm()
        return pr_vec  

    def sample_size_at_i(self,i):
        if i not in self.bpoints: return 0
        s = self.bpoints[i].shape 
        return s[0] * s[1] 

    def check_fstat(self,dimset):
        self.fstat = set(self.tdim) == dimset 

class Sec:

    def __init__(self,sequence,singleton_range,\
        optima_pr_map,dep_map,codep_map,obfsr=None):

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

        self.obfsr = obfsr
        if type(self.obfsr) == type(None):
            self.declare_obf()
        
        self.declared_dir = [] 

    def optima_points(self):
        ks = sorted(list(self.opm.keys()))
        optima_points = [matrix_methods.string_to_vector(v,float) for v in \
                    ks]
        optima_points = np.array(optima_points)
        return optima_points

    def optima_points_to_index_pr_map(self):

        ks = sorted(list(self.opm.keys()))
        ks = [(i,self.opm[k]) for (i,k) in enumerate(ks)]

        d = defaultdict(float,ks) 
        return d

    def declare_obf(self):

        abf = default_AltBaseFunc_for_IsoRing()
        bloom_func = abf.load_AltBaseFunc_function(abf)

        # generate the size-delta function for
        # the <IsoRingBloomer>
        multiplier = 3.13
        szm = generate_default_szmult_mod(multiplier)

        # initialize the <OptimaBloomFunc>
        optima_points = self.optima_points()
        bounds = np.array([deepcopy(self.singleton_range) for _ in \
            range(len(self.seq))])

        obf = OptimaBloomFunc(optima_points,bounds,\
            None,bloom_func,optima_points.shape[0],8)
        self.obfsr = OptimaBloomFuncSecRep(obf,szm,\
            optima_points.shape[1])
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

    def process_one_bloomiso(self):
        stat = True 
        while stat: 
            # next value, indices of derivators
            b1,b2 = self.__next__()
            stat = not (type(b1) == type(None))
            if not stat:
                ##self.declare_another_Sec(self.obfsr.sz)
                # case: obfsr needs to be reset
                self.obfsr.set_sz()
            return -1 

        return -1 

    def lone_pr_vec_for_bloom(self):
        x = self.optima_points_to_index_pr_map()
        return self.obfsr.finalize_dcount(x)

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
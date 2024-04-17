from ep_corrmap import * 
from bloominhurst import *

default_rnd_boolean_index_splitter = lambda x: True if random.randrange(0,2) else False 

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
        # finished blooming
        self.fstat = False
        # terminate at dim `sz` 
        self.tstat = False
        self.set_sz()

        self.set_dpm() 

    def set_dpm(self):
        self.dpm = DerivatorPrMap() 

    def set_sz(self,attempts=10):
        # get dim of OptimaBloomFunc
        self.sz = self.obf.oseeds.shape[1]
        self.tdim.append(self.sz) 

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

        #if self.sz in self.tdim:
        #    return None

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
        ##print("SZ IS: ",self.sz)
        pr_vec,dep_map = self.dpm.fin_count(self.corr_type,\
            self.sz,pred_opt2pr_map)
        self.set_dpm()
        return pr_vec,dep_map 

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
        assert matrix_methods.vector_to_string(sequence,float) in optima_pr_map

        self.seq = sequence
        self.singleton_range = singleton_range
        self.opm = optima_pr_map
        self.dm = dep_map
        self.cdm = codep_map

        self.obfsr = obfsr
        if type(self.obfsr) == type(None):
            self.declare_obf()
        
        self.declared_dir = []
        self.next_sz = None 

        self.dimso = None 

    def __str__(self):
        s = "** sequence\n"
        s += matrix_methods.vector_to_string(self.seq,float)
        s += "\n" + "** optima pr." + "\n"
        s += str(self.opm)
        s += "\n" + "** dep. map" + "\n"
        s += str(self.dm)
        s += "\n" + "** co-dep. map" + "\n"
        s += str(self.cdm)

        return s + "\n"

    def optima_points(self):
        ks = sorted(list(self.opm.keys()))
        optima_points = [matrix_methods.string_to_vector(v,float) for v in \
                    ks]
        optima_points = np.array(optima_points)
        return optima_points

    """
    converts opm map 
        key := stringized point
        value := Pr. value
    to a map w/ 
        key := index of ordering for stringized point
        value := Pr. value 
    """
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
        
        self.dimso = len(q)
        return q, q2

    def process_one_bloomiso(self):
        stat = True 
        while stat: 
            ##print("one bloom pt.")
            # next value, indices of derivators
            b1,b2 = self.__next__()
            ##
            """
            print("SEC NEXT")
            print(b1)
            print()
            print(b2)
            print("-----")
            """
            ##
            stat = not (type(b1) == type(None))
            if not stat:
                continue
        return

    # TODO: bug here. 
    def lone_pr_vec_for_bloom(self):
        
        x = self.optima_points_to_index_pr_map()
        return self.obfsr.finalize_dcount(x)
        
        #####

    def generate_next_Sec(self):
        assert not self.obfsr.fstat
        assert self.obfsr.tstat 
        
        q = self.lone_pr_vec_for_bloom()

        # get the optima points for the `sz`
        # dimension
        bps = deepcopy(self.obfsr.bpoints[self.obfsr.sz])
        assert bps.shape[0] == len(q[0])

        optima_pr_map = defaultdict(float)
        for i,q_ in enumerate(q[0]):
            vs = matrix_methods.vector_to_string(bps[i],float)
            optima_pr_map[vs] = np.round(q_,5)

        # get the dep. maps
        #   split the dependency map w/ co-dep map
        dep_map,codep_map = {},{}
        for (k,v) in q[1].items():
            stat = default_rnd_boolean_index_splitter(v)
            if stat: 
                dep_map[k] = v
            else: 
                codep_map[k] = v
        dep_map,codep_map = defaultdict(float,dep_map),\
            defaultdict(float,codep_map)

        # get the new <seq> vector 
            # Pr. of old ans. 
        vsm = matrix_methods.vector_to_string(self.seq,float)
        pr = self.opm[vsm]

        qopm = [(k,v) for k,v in optima_pr_map.items()] 
        qopm = sorted(qopm,key=lambda x: x[0])
        nu_pt = min(qopm,key=lambda x: abs(x[1] - pr)) 
        nu_pt = nu_pt[0]
        nu_pt = matrix_methods.string_to_vector(nu_pt,float)
        #nu_pt = np.round(nu_pt,5)

        s = Sec(nu_pt,deepcopy(self.singleton_range),\
            optima_pr_map,dep_map,codep_map,obfsr=self.obfsr)
        sz0 = s.obfsr.sz
        ##
        """
        print("-- RESETTING OSEEDS")
        print("\t-- ORIGINAL")
        print(self.obfsr.obf.oseeds)
        print("\t-- NEW")
        print(bps)
        print()
        """
        ##
        sr_mat = np.array([deepcopy(self.singleton_range) for _ in range(sz0)])
        self.obfsr.obf.reset_oseeds(bps,sr_mat,True)
        self.obfsr = None
        s.obfsr.set_sz()
        s.obfsr.tstat = False 
        return s,sz0,s.obfsr.sz
        

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
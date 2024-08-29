# functions to modulate index vectors
from morebs2 import relevance_functions,aprng_gauge,\
    poly_interpolation,matrix_methods
from custom_rnd import *

def roundian(x,r):
    rd = r[1] - r[0]
    q0 = int(int(round(x)) % rd + r[0])
    return q0 

############# functions to generate starting points for 
############# <DefaultLPSMod> 

def generate_2d_points__prng(x0_range,x1_range,prng_func,num_points:int):
    assert len(x0_range) == 2 and x0_range[0] < x0_range[1]
    assert len(x1_range) == 2 and x1_range[0] <= x1_range[1]

    points = []
    for i in range(num_points):

        # get the ratio for each of the 
        # two dimensions
        r1 = prng_func()
        r2 = prng_func()
        assert min([r1,r2]) >= 0.0 and max([r1,r2]) <= 1.0
        
        # get the 2-d point
        p1 = (x0_range[1] - x0_range[0]) * r1 + x0_range[0]
        p2 = (x1_range[1] - x1_range[0]) * r2 + x1_range[0]
        points.append((p1,p2))
    
    points.append([x0_range[0],x1_range[0]])
    points.append([x0_range[1],x1_range[1]])
    points = np.array(points)
    return np.sort(points,axis = 0)

########################################################### preliminary functions 
"""
outputs all points in range specified by default 
hop value of <morebs2>. 
"""
def generate_correspondent_func_type__lagrange_poly(ndim=2,\
    frange=deepcopy(DEFAULT_SINGLETON_RANGE),\
    rnd_seed_int:int=None,num_points = 8):
    assert len(frange) == 2 and frange[0] <= frange[1]

    if type(rnd_seed_int) == int:
        random.seed(rnd_seed_int)

    pointseq = generate_2d_points__prng(deepcopy(frange),deepcopy(frange),\
        random.random,num_points)
    ##print("pointseq")
    ##print(pointseq)
    ##print()

    lps = poly_interpolation.LagrangePolySolver(pointseq, prefetch = True)

    is_reverse = bool(random.randrange(0,2))
    pseq = list(lps.form_point_sequence(reverse = is_reverse))

    def out_func():
        if len(pseq) == 0:
            return None
        #rindex = random.randrange(0,int(len(pseq) / 2))
        rindex = 0 
        r = pseq.pop(rindex)
        return np.array([r[0] % (frange[1] - frange[0]) + frange[0],\
            r[1] % (frange[1] - frange[0]) + frange[0]])

    return out_func 

"""
generates a modulator that outputs

npstd_func(x,corr_func()); in 2d-range [dim0,dim1]. 

given input x. 

corr_func := corresponding function that outputs the 
            `y` value for input into numpy std. function
"""
def generate_2dmod__np_std(dim0,dim1,npstd_func,corr_func=None):
    assert len(dim0) == 2 
    assert len(dim0) == len(dim1)

    if type(corr_func) == type(None):
        corr_func = generate_correspondent_func_type__lagrange_poly()
    
    def f(x):
        c = corr_func()
        if type(c) != type(None):
            q = npstd_func(x,c)
            ##print(">>> Q: ", q)
            q0 = DEFAULT_REPLACE_INF(q[0])
            q1 = DEFAULT_REPLACE_INF(q[1])
        else:
            q0 = 0.0
            q1 = 0.0 
        ##print(">>> q0: ", q1, " q1: ", q1)
        q0 = roundian(q0,dim0)
        q1 = roundian(q1,dim1)
        return q0,q1
    return f

class DefaultLPSMod:

    def __init__(self,lps,drange,splitsz,splitsz_delta):
        self.lps = lps
        assert matrix_methods.is_proper_bounds_vector(drange)
        assert drange.shape == (2,2)
        assert type(splitsz) == int and splitsz > 0
        self.drange = drange
        self.splitsz = splitsz
        self.splitsz_delta = splitsz_delta
        self.splitf = 0.0

        if type(self.lps) == type(None):
            self.generate_lps()

        # set the point reference
        self.pref = None
        self.set_split()

    def generate_lps(self,num_points=8):
        pointseq = generate_2d_points__prng(deepcopy(self.drange[0]),\
            deepcopy(self.drange[1]),random.random,num_points)
        lps = poly_interpolation.LagrangePolySolver(pointseq, prefetch = True)
        self.lps = lps
        return deepcopy(lps)

    """
    only on axis 0. 
    """
    def set_split(self):
        """
        self.splitf = np.array([self.drange[0,1] - \
            self.drange[0,0],\
            self.drange[1,1] - self.drange[1,0]])
        self.splitf = self.splitf / self.splitsz
        """
        ##### 
        rx = self.lps.bounds_for_x("float")
        self.splitf = np.round(abs(rx[1] - rx[0]) / self.splitsz,5)
        self.pref = self.lps.bounds_for_x("float")[0]

    def reset_split(self):
        self.splitsz = self.splitsz_delta(self.splitsz)
        self.set_split()

    def __next__(self):
        # check if `pref` exceeds x-bounds
        rx = self.lps.bounds_for_x("float")
        if self.pref > rx[1]:
            self.reset_split()

        # get the value
        y = self.lps.output_by_lagrange_basis(self.pref)
        px = (self.pref,y)

        self.pref = self.pref + self.splitf
        return px

    @staticmethod
    def load_DefaultLPSMod_function(lpsm,npstd_func):
        assert type(lpsm) == DefaultLPSMod

        def de_corr_func():
            return next(lpsm)

        cf = de_corr_func
        mod2d_fuunc = generate_2dmod__np_std(\
            deepcopy(lpsm.drange[0]),\
            deepcopy(lpsm.drange[1]),\
            npstd_func,\
            corr_func=cf)
        return mod2d_fuunc 

########################## default function for mapping 
########################## 2-d point to 2-d space. 

class Index2DMod:

    def __init__(self,dim2d,bis,f=None):
        assert matrix_methods.is_proper_bounds_vector(dim2d)
        assert np.min(dim2d) >= 0 
        assert type(bis) == morebs2.aprng_gauge.BatchIncrStruct
        self.dim2d = dim2d 
        self.bis = bis

        # case: generate a 2d-modulator function
        if type(f) == type(None):
            ##print("generating 2dm func.")
            q = deepcopy(DEFAULT_PAIRWISE_VEC_FUNCS)
            i = random.randrange(len(q))
            npstd_func = q[i]
            mod2d = generate_2dmod__np_std(self.dim2d[0],self.dim2d[1],npstd_func,None)
            self.f = mod2d 
        else: 
            self.f = f
        
        self.fin_stat = False
        return

    def __next__(self):
        if self.fin_stat: 
            return None

        q = next(self.bis)
        ##print("q: ",q)

        if type(q) == type(None): 
            # case: terminated. 
            self.fin_stat = True
            return q

        q2 = self.f(q)
        q2_ = [roundian(q2[0],self.dim2d[0]),\
            roundian(q2[1],self.dim2d[1])]
        return q2_ 


class IndexVecPermuter:

    def __init__(self,vecsize_seq,head=0):
        assert head in {0,-1}

        self.vss = vecsize_seq
        self.head = head 
        self.p = [0 for _ in range(len(self.vss))]
        self.initialized = False
        self.finished = False

    def __next__(self):
        if self.finished:
            return None

        """
        if set(self.p) == {0} and self.initialized:
            self.finished = True
            return None
        """
        
        px = deepcopy(self.p)


        index = 0 if self.head == 0 else len(self.p) - 1
        index_delta = 1 if self.head == 0 else -1
        px = self.next__(px,index,index_delta)
        
        px2 = deepcopy(self.p)
        self.p = px
        return px2

    def next__(self,px,index,index_delta):
        """
        if index == -1: 
            index = len(px) - 1
        elif index >= len(self.p):
            index = 0
        """
        if index == -1 or index >= len(self.p): 
            self.finished = True
            return None

        px[index] += 1

        if px[index] >= self.vss[index]:
            px[index] = 0
            index = index + index_delta
            return self.next__(px,index,index_delta)
        return px






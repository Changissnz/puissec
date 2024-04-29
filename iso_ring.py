from obj_func import * 
from secnet_gen import * 
from cvec import * 

class BoundedObjFunc:

    """
    takes a sequence of proper bounds `bounds_seq` w/
    each bound corresponding to a different i'th element
    in `corr_obj`. 

    Main method is `output(reference point, current point)`.
    If `current point` is not in any bounds of `bounds_seq`,
    output function is `default_objf`.
    """
    def __init__(self,bounds_seq,corr_objf,default_objf):

        for c in corr_objf: assert type(c) == ObjFunc
        for b in bounds_seq: assert matrix_methods.\
            is_proper_bounds_vector(b) 
        assert len(bounds_seq) == len(corr_objf)
        assert type(default_objf) == ObjFunc

        self.bounds_seq = bounds_seq
        self.corr_objf = corr_objf
        self.dobjf = default_objf

    def __str__(self): 
        s = ""

        for (i,b) in enumerate(self.bounds_seq):
            s += "bounds\n"
            s += str(b)
            s += "\n"
            s += str(self.corr_objf[i])
            s += "\n"
        s += "--------------"
        return s  

    @staticmethod
    def generate_BoundedObjFunc(superbound,\
        spacing_ratio_range,outlier_pr_ratio,\
        num_bounds,rnd_seed:int):

        assert type(rnd_seed) == int 
        assert type(num_bounds) == int and num_bounds > 0 

        r = {0:"e.d.",1:"s.m.",2:"r.n."}
        random.seed(rnd_seed)

        bounds_seq = generate_bounds_vector_sequence(\
            superbound,spacing_ratio_range,outlier_pr_ratio,num_bounds)

        corr_objf = []
        for i in range(num_bounds):
            # choose a random obj_type
            obj_type = random.randrange(0,3)
            obj_type = r[obj_type]

            s2 = random.randrange(100)
            objf = ObjFunc(obj_type,s2)
            corr_objf.append(objf)
        
        default_objf = deepcopy(corr_objf[0])

        bof = BoundedObjFunc(bounds_seq,corr_objf,\
            default_objf)
        return bof 

    def output(self,ref,x):
        i = self.bounds_index_of_point(x)
        if i == -1:
            cf = self.dobjf
        else: 
            cf = self.corr_objf[i]
        return cf.output(ref,x)

    """
    return: 
    - index of bounds that `x` falls in OR -1 
    """
    def bounds_index_of_point(self,x):
        for (i,b) in enumerate(self.bounds_seq):
            if matrix_methods.point_in_bounds(b,x):
                return i
        return -1 

class IsoRing:

    def __init__(self,sec:Sec,ofunc:BoundedObjFunc,bounds,\
        cvecl=None):
        assert type(sec) == Sec
        assert type(ofunc) == BoundedObjFunc
        assert matrix_methods.is_proper_bounds_vector(bounds)
        if type(cvecl) != type(None):
            for c in cvecl: assert type(c) == CVec
            assert len(cvecl) == len(sec.opm)

        self.sec = sec
        self.sec_cache = [self.sec] 

        self.ofunc = ofunc 
        self.bounds = bounds 
        # bloom stat 
        self.bstat = True
        # index of repr in `sec_cache`
        self.repi = 0

        # cracked stat 
        self.cstat = False 

    def explode_contents(self,optima_size_limit=1000):
        s = len(self.sec_cache[-1].opm)
        ##print("starting length for {}: {}".format\
        ##    (len(self.sec_cache[-1].seq),s))
        while s < optima_size_limit:
            sc = self.sec_cache.pop(-1)
            sc.process_one_bloomiso()
            s2 = sc.generate_next_Sec()
            ts2 = s2[2]
            ##print("transition: {}->{}".format(s2[1],s2[2]))
            s2 = s2[0]
            s = len(s2.opm)
            self.sec_cache.append(sc)
            self.sec_cache.append(s2)
            if type(ts2) == type(None):
                break
        return

    def register_attempt(self,p):
        ##print("-- RA")
        q = self.register_attempt_(p)
        stat = matrix_methods.equal_iterables(p,self.sec.seq)
        outstat = []
        ops = self.rep().optima_points()
        for q in ops:
            stat2 = matrix_methods.equal_iterables(p,q) 
            outstat.append(stat2)
        ##print("got outstat")
        ##print(outstat)
        outstat = np.any(np.array(outstat) == True) 
        self.cstat = stat
        return q,outstat 


    def register_attempt_(self,p):
        assert matrix_methods.is_vector(p) 
        r = self.rep()
        ops = r.optima_points()
        if len(p) != ops.shape[1]:
            print("attempt in wrong dim {}, want {}".format(len(p),\
                ops.shape[1]))
            return None 
        
        # get scores for each of the optima 
        q = list(map(lambda x: self.ofunc.output(x,p),ops))
        ##print("-- reg: ",q)
        return np.array(q)

    def set_isorep(self,i):
        assert len(self.sec_cache) > i and i > -1
        self.repi = i 
        return

    def rep(self):
        return deepcopy(self.sec_cache[self.repi])

    def pr_of_optima_index(self,oi):

        return -1 

def Sec_list_sample1(): 
    random.seed(12)
    np.random.seed(12)

    singleton_range = [0.,1.] 
    dimension = 5
    num_optima = 2
    countermeasure = (0.6,0.5) 
    secs = []

    for i in range(5): 
        sec = Sec.generate_bare_instance(singleton_range,dimension,num_optima,\
        countermeasure,rnd_struct=np.random)
        secs.append(sec)
    return secs 

def Sec_list_sample2(): 
    random.seed(14)
    np.random.seed(19)

    singleton_range = [0.,1.] 
    dimension = 5
    num_optima = 12
    countermeasure = (0.7,0.5) 
    secs = []

    for i in range(num_optima): 
        sec = Sec.generate_bare_instance(singleton_range,dimension,num_optima,\
        countermeasure,rnd_struct=np.random)
        secs.append(sec)
    return secs 


### an example of an <IsoRing> instantiation
def IsoRing_sample_1():
    secs = Sec_list_sample1()
    sndg = SecNetDepGen(secs,random,2,0.75,[1,4])
    sndg.assign_conn()
    sq = sndg.sq

    superbound = np.ones((5,2)) * np.array([0.,1.])
    spacing_ratio_range = [0.,0.2]
    outlier_pr_ratio = 0.4#1.0
    num_bounds = 8
    sb = deepcopy(superbound)

    obf = BoundedObjFunc.generate_BoundedObjFunc(\
            superbound,spacing_ratio_range,\
            outlier_pr_ratio,num_bounds,3) 

    return IsoRing(sq[0],obf,sb)
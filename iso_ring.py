from obj_func import * 
from sec_seq import * 
from cvec import * 

class BoundedObjFunc:

    def __init__(self,bounds_seq,corr_objf,default_objf):

        for c in corr_objf: assert type(c) == ObjFunc
        for b in bounds_seq: assert matrix_methods.\
            is_proper_bounds_vector(b) 
        assert len(bounds_seq) == len(corr_objf)
        assert type(default_objf) == ObjFunc

        self.bounds_seq = bounds_seq
        self.corr_objf = corr_objf
        self.dobjf = default_objf

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
        # sequence of <CVec> instances, each 
        # i'th <CVec> corresponds to the i'th
        # 
        self.cvecl = cvecl
        self.declare_cvecl()

        # cracked stat 
        self.cstat = False 

    def declare_cvecl(self):
        if type(self.cvecl) != type(None):
            return

        self.cvecl = [] 
        for i in range(0,len(self.sec.opm)):
            cvseq = default_cvec_iselector_seq()
            cv = CVec(cvis=cvseq)
            self.cvecl.append(cv)

    def operating_sec(self): 
        if self.bstat:
            return self.sec_cache[-1]
        return self.sec_cache[self.repi] 

    def run_Sec_till_fstat(self,guesser_func):
        sec = self.operating_sec()
        stat = sec.obfsr.fstat

        while stat:
            # get the length of output 
            p = guesser_func(sec.obfsr.sz)
            q,s = self.register_attempt(p)
            print("guess is")
            print(q)
            print(s)
            print()

            sec = self.operating_sec()
            stat = sec.obfsr.fstat 

        print("done")

    def register_attempt(self,p):
        sq = self.sec_cache[-1]
        stat1 = sq.obfsr.fstat
        stat2 = sq.obfsr.tstat 

        # case: generate an element for next repr. 
        if not stat2: 
            x1,x2 = next(sq) 

            if type(x1) == type(None):
                print("recursing on register")
                return self.register_attempt(p)
        #       declare another <Sec> instance
        else:
            sec2,previous_sz,current_sz = sq.generate_next_Sec()
            self.sec_cache.append(sec2) 

        # calculate feedback for attempt
        q = self.register_attempt_(p)
        assert len(q) == len(self.cvecl)

        # log the feedback into each of the CVecs 
        self.cvecl[0].append(q[0],p)
        for i in range(1,len(self.cvecl)): 
            self.cvecl[i].append(q[i],None)

        # get the actual stat
        stat = matrix_methods.equal_iterables(p,self.sec.seq)
        self.cstat = stat 
        return q,stat

    """
    return: 
    - vector of distance calculations, finished status 
    """ 
    def register_attempt_(self,p):
        assert matrix_methods.is_vector(p) 
        ops = self.sec.optima_points() 
        if len(p) != ops.shape[1]:
            print("attempt in wrong dim {}, want {}".format(len(p),\
                ops.shape[1]))
            return None 
        
        # get scores for each of the optima 
        q = list(map(lambda x: self.ofunc(x,p),ops))
        return np.array(q)

    def cvec_feedback(self):
        x = []
        for c in self.cvecl:
            x.append(c.cmp())
        return np.array(x,dtype=bool)

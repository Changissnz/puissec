from obj_func import * 
from sec_seq import * 

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

    def __init__(self,sec:Sec,ofunc:BoundedObjFunc,bounds):
        assert type(sec) == Sec
        assert type(ofunc) == BoundedObjFunc
        assert matrix_methods.is_proper_bounds_vector(bounds)
        self.sec = sec
        self.sec_cache = [self.sec] 

        self.ofunc = ofunc 
        self.bounds = bounds 

    def register_attempt(self,p):
        sq = self.sec_cache[-1]
        stat1 = sq.obfsr.fstat
        stat2 = sq.obfsr.tstat 

        # case: generate an element for next repr. 
        if not stat2: 
            x1,x2 = next(sq) 
        #       declare another <Sec> instance
        else:
            sec2 = sq.generate_next_Sec()
            self.sec_cache.append(sec2) 

        # calculate feedback for attempt
        q = self.register_attempt_(p)

        # get the actual stat
        stat = matrix_methods.equal_iterables(p,self.sec.seq)
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


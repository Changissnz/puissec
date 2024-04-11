from imod import *

DEFAULT_OPTIMA_BLOOM_SZ_DELTA_FUNC = lambda x: int(x * 2)

class OptimaBloomFunc:

    """
    Uses a sequence of vector::float S1 to output another 
    sequence of vector::float S2 so that dim(S1) != dim(S2)
    possible. 

    Procedure first takes a random sampling of elements 
    from `optima_seedlings`. Then it produces the remaining 
    elements from an accumulation of single pair inputs,
    of size equal to `target_vsize`. 

    optima seedlings := sequence, of vectors in real numbers R.
    bloom_func := AltBaseFunc, operates on the integers
    """
    def __init__(self,optima_seedlings,drange,selector_func,bloom_func,\
        d:int,split_sz:int,splitsz_delta=DEFAULT_OPTIMA_BLOOM_SZ_DELTA_FUNC):
        assert type(d) == int and d > 0
        assert type(split_sz) == int and split_sz > 0
        assert type(selector_func) in {Index2DMod,type(None)}

        self.oseeds = None
        self.drange = None
        self.selector_func = selector_func
        # function that takes two float values as input
        self.bloom_func = bloom_func
        # wanted dimension of vector output
        self.d = d
        self.split_sz = split_sz
        self.splitsz_delta = splitsz_delta
        self.reset_oseeds(optima_seedlings,drange,False)

        if type(selector_func) == type(None):
            self.declare_sfunc()

        self.finished_stat = False
        self.prev_i1,self.prev_i2 = None,None

    # declares a selector function
    def declare_sfunc(self):
        # declare <BatchIncrStruct>
        m_ = max(self.oseeds.shape)
        ##print("batch incrementor on: [{},2]".format(m_))
        bis = aprng_gauge.BatchIncrStruct(m_,True,True,2)

        # declare <DefaultLPSMod>
        drange = np.array([(0,self.oseeds.shape[0]),\
            (0,self.oseeds.shape[1])])
            
        lpsm = DefaultLPSMod(None,deepcopy(drange),\
            self.split_sz,self.splitsz_delta)
        fx = DefaultLPSMod.load_DefaultLPSMod_function(lpsm,\
            np.add)
        i2dm = Index2DMod(drange,bis,f=fx)
        self.selector_func = i2dm
        return i2dm

    def reset_oseeds(self,optima_seedlings,drange,update_sfunc:bool):
        assert type(optima_seedlings) == np.ndarray
        assert optima_seedlings.ndim == 2 and \
            min(optima_seedlings.shape) > 0 
        assert type(drange) == np.ndarray
        assert drange.shape[0] == optima_seedlings.shape[1]
        assert matrix_methods.is_proper_bounds_vector(drange)

        self.oseeds = optima_seedlings
        self.drange = drange

        if update_sfunc: 
            self.declare_sfunc() 

    def shape(self):
        return self.oseeds.shape

    """
    return:
    - a vector 
    """
    def __next__(self):

        self.prev_i1,self.prev_i2 = [],[]

        def f():
            if self.finished_stat:
                print("ERROR: no more output!")
                return 

            # select two indices
            i1 = next(self.selector_func)
            i2 = next(self.selector_func)
            self.prev_i1.append(i1)
            self.prev_i2.append(i2)

            ##print("next indices: {} and {}".format(i1,i2))

            stat1 = type(i1) == type(None)
            stat2 = type(i2) == type(None)

            if stat1 or stat2:
                self.finished_stat = True 
                return

            # retrieve the values of the 
            # two indices            
            v1 = deepcopy(self.oseeds[i1[0],i1[1]])
            v2 = deepcopy(self.oseeds[i2[0],i2[1]])
            q = self.bloom_func(v1,v2)

            ##
            """
            print("OSEEDS")
            print(self.oseeds)
            print("OSEED INDEX: ", i1, i2)
            print("OSEED ON:\n{}\n{}\n".format(v1,v2))
            print("performing bloom on: {}, {} --> {}".format(v1,v2,q))
            """
            ##
            
            return q 

        def g():
            v = []
            for i in range(self.d):
                q = f()
                if type(q) == type(None):
                    break

                i_ = i % self.drange.shape[0]
                q = q % (self.drange[i_,1] - self.drange[i_,0]) +\
                    self.drange[i_,0]
                v.append(q)

            if len(v) != self.d:
                return None 
                v = [0. for _ in range(self.d)]
            ##print("bloom value: ",v)
            return np.array(v) 

        return g() 
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
        d:int,split_sz:int,splitsz_delta=DEFAULT_OPTIMA_BLOOM_SZ_DELTA_FUNC,\
        ):
        assert type(d) == int and d > 0
        assert type(split_sz) == int and split_sz > 0
        assert type(selector_func) in {Index2DMod,type(None)}

        self.oseeds = None
        self.drange = None
        self.reset_oseeds(optima_seedlings,drange)
        self.split_sz = split_sz
        self.splitsz_delta = splitsz_delta

        # function that takes two float values as input
        self.bloom_func = bloom_func
        # wanted dimension of vector output
        self.d = d
        self.selector_func = selector_func
        if type(selector_func) == type(None):
            self.declare_sfunc()

        self.finished_stat = False

    # declares a selector function
    def declare_sfunc(self):
        # declare <BatchIncrStruct>
        m_ = max(self.oseeds.shape)
        bis = aprng_gauge.BatchIncrStruct(m_,True,True,2)

        # declare <DefaultLPSMod>
        lpsm = DefaultLPSMod(None,deepcopy(self.drange),\
            self.split_sz,self.splitsz_delta)
        fx = DefaultLPSMod.load_DefaultLPSMod_function(lpsm,\
            self.bloom_func)
        i2dm = Index2DMod(self.drange,bis,f=fx)
        self.selector_func = i2dm
        return i2dm

    def reset_oseeds(self,optima_seedlings,drange):
        assert type(optima_seedlings) == np.ndarray
        assert optima_seedlings.ndim == 2 and \
            min(optima_seedlings.shape) > 0 
        assert type(drange) == np.ndarray
        assert drange.shape[0] == optima_seedlings.shape[1]
        assert matrix_methods.is_proper_bounds_vector(drange)

        self.oseeds = optima_seedlings
        self.drange = drange

    def __next__(self):

        def f():
            if self.finished_stat:
                print("ERROR: no more output!")
                return 

            # select two indices
            i1 = next(self.selector_func)
            i2 = next(self.selector_func) 

            stat1 = type(i1) == type(None)
            stat2 = type(i2) == type(None)

            if stat1 or stat2:
                self.finished_stat = True 
                return

            # retrieve the values of the 
            # two indices
            v1 = deepcopy(self.oseeds[i1[0],i1[1]])
            v2 = deepcopy(self.oseeds[i2[0],i2[1]])
            return self.bloom_func(v1,v2)
        return f() 
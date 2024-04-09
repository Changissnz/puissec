"""
for use with single elements of a secret
"""
from sec_seq import * 
from defaults import * 
from obj_func import * 
from bloominhurst import *

class IsoRingBloomer:

    def __init__(self,obf):
        assert False
        return -1 

class IsoRing:

    def __init__(self,contents,rnd_seed = None,\
                 local_optima=None, obj_func=None,\
                    is_entry_point=False):
        assert morebs2.is_proper_bounds_vector(contents)
        self.contents = contents
        self.rnd_seed = rnd_seed
        if type(self.rnd_seed) == type(int):
            random.seed(self.rnd_seed)
            np.random.seed(self.rnd_seed) 
        self.local_optima = local_optima
        self.ofunc = obj_func
        self.is_entry_point = is_entry_point
        self.instantiate()

    def instantiate(self):
        print("instantiating")
        if type(self.local_optima) == type(None):
            print("generating L.O.")
            self.generate_pr_local_optima()
        
        if type(self.ofunc) == type(None):
            return "TODO"
        return 

    def generate_pr_local_optima(self):
        print("-->")
        nr = random.randrange(DEFAULT_ISO_RING_LOCAL_OPTIMA_SIZE_RANGE[0],\
            DEFAULT_ISO_RING_LOCAL_OPTIMA_SIZE_RANGE[1] + 1)
        self.local_optima = []
        print("NR: ",nr)
        for i in range(nr):
            p = self.generate_pr_local_optimum()
            self.local_optima.append(p)
        self.local_optima = np.array(self.local_optima)
        
        return

    """
    """
    def generate_pr_local_optimum(self):
        l = len(self.contents)

        diff = self.contents[::,1] - self.contents[::,0]
        ##
        x = np.random.choice(deepcopy(DEFAULT_SINGLETON_RANGE),\
            (l,)) 
        additive = np.round(diff * x,5)
        return np.round(self.contents[::,0] + additive,5)

    def register_hypothesis(self,p):
        i = self.closest_optimum(p)
        r = self.local_optima[i,:]
        # TODO: check this
        return self.ofunc(r,p)

    def closest_optimum(self,p):
        assert len(p) == len(self.contents) 

        # get the closest local optimum to p
        x = np.array([morebs2.matrix_methods.euclidean_point_distance(\
            lo,p) for lo in self.local_optima])
        i = -1
        if len(x) > 0:
            i = np.argmin(x)
        return i

    """
    mark if the previous cracking attempt 
    is detected by the instance's owner. 
    """
    def is_detected(self):
        return -1
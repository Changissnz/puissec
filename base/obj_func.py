from morebs2 import matrix_methods,modular_labeller
from .defaults import *

class ObjFunc:
    
    """
    a class with main method `output` that operates on 
    2 vectors, a reference and another point.

    The `obj_type` is one of {"e.d.","s.m.","r.n."}.
    The rnd_seed is either an integer or None. 
    """
    def __init__(self,obj_type,rnd_seed=None):
        # e.d. := Euclidean distance
        # s.m. := scalar modular 
        # r.n. := random noise
        assert obj_type in {"e.d.","s.m.","r.n."}
        #assert morebs2.matrix_methods.is_proper_bounds_vector(target_range)

        self.obj_type = obj_type
        #self.target_range = target_range
        self.rnd_seed = rnd_seed
        if type(self.rnd_seed) == int: 
            random.seed(self.rnd_seed)
            np.random.seed(self.rnd_seed) 
        self.instantiate()
        return

    def __str__(self):
        s = "OFUNC\n\t- {}\n\t- {}".format(self.obj_type,self.rnd_seed)
        return s 
        
    def instantiate(self):
        if self.obj_type == "e.d.":
            self.ofunc = matrix_methods.euclidean_point_distance
        elif self.obj_type == "s.m.":
            self.ldp = self.generate_random_LDP()
            self.instantiate_sm_func()
        else:
            self.instantiate_rn_func() 
        return

    def instantiate_sm_func(self):
        self.ofunc = lambda ref,point: matrix_methods.euclidean_point_distance(ref,point) * self.ldp.output_label()

    def instantiate_rn_func(self):
        self.ofunc = lambda ref,point: matrix_methods.euclidean_point_distance(ref,\
                        point + \
                        ((np.zeros(len(ref),) + DEFAULT_RN_RANDOM_NOISE_RANGE[0]) + \
                        (np.random.rand(len(ref),) * (DEFAULT_RN_RANDOM_NOISE_RANGE[1] - \
                            DEFAULT_RN_RANDOM_NOISE_RANGE[0]))))
        return

    def generate_random_LDP(self):
        sl = random.randrange(DEFAULT_SM_LABEL_RANGE[0],\
            DEFAULT_SM_LABEL_RANGE[1] + 1)

        tl = random.choice(DEFAULT_SM_TICKLENGTH_RANGE)
        ldp = modular_labeller.LabelDeltaPattern(tl,\
            deepcopy(DEFAULT_SM_DEVIATION_PR),\
            deepcopy(DEFAULT_SM_DEVIATION_RANGE),\
            deepcopy(DEFAULT_SM_LABEL_RANGE),\
            sl)
        return ldp 

    def output(self,reference_point,point):
        return round(self.ofunc(reference_point,point),5)
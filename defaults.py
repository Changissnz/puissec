import numpy as np 
from copy import deepcopy 
import random 
from collections import defaultdict,Counter 
import math 
import morebs2 
import pickle 


############## default variables

DEFAULT_SINGLETON_RANGE = [-0.0,1.0]

DEFAULT_SM_TICKLENGTH_RANGE = [1,10]
DEFAULT_SM_DEVIATION_PR = 0.15
DEFAULT_SM_DEVIATION_RANGE = [-2,2]
DEFAULT_SM_LABEL_RANGE = [-10,10]

DEFAULT_RN_RANDOM_NOISE_RANGE = [-0.2,0.2]

DEFAULT_ISO_RING_LOCAL_OPTIMA_SIZE_RANGE = [2,12]
DEFAULT_EDGE_COST_FUNCTION = lambda u,v,c: 1
CUMULATIVE_EDGE_COST_FUNCTION = lambda u,v,c: 1 + c  
CUMULATIVE_PATH_COST_FUNC = lambda x: len(x) 

DEFAULT_PAIRWISE_VEC_FUNCS = [\
        np.subtract,\
        np.add,\
        np.multiply,\
        np.divide]

DEFAULT_REPLACE_INF = lambda x: round(x,5) if (not np.isnan(x) \
    and not np.isinf(abs(x))) else 0.0 

DEFAULT_TERMINATE_DETECT_FUNC = lambda l,l2: not (l == l2).all() if \
    type(l2) != type(None) else True  

DEFAULT_OBF_SECREP_BLOOM_SZ_LIMIT = 1000
DEFAULT_OBF_OPTIMA_SIZE_RANGE = [1,50] 

idn_fx = lambda x: x 
sqrt_fx = lambda x: math.sqrt(x)

# DEFAULT_BLACKBOX_FUNCTIONS 
blackbox_df1 = lambda x: (x + random.random()) % 1.0


############## miscellaneous functions

def pickle_open_with_typecheck(fp,t): 
    fobj = open(fp,"rb")
    obj = pickle.load(fobj)
    assert type(obj) == t
    return obj 

def strsplit_float(f):
    assert type(f) in {np.float64,float}
    dx = str(f)
    q = dx.split(".") 
    return q[0],q[1] 

# TODO: relocate
def all_multiples_decimal(i,rounding_depth=5):
    qi = round(i,rounding_depth)
    qi = str(qi)
    qi = qi.split(".")
    l = len(qi[1])
    x = qi[0] + qi[1]
    qi = int(x) 
    multiples = morebs2.numerical_extras.all_multiples(qi)
    multiples = sorted(multiples)
    for i in range(len(multiples)):
        multiples[i] = multiples[i] * 10 ** -l
    return multiples

#all_mu
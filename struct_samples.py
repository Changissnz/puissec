from sec_seq import * 
from secnet_gen import *

def Sec_sample_1():
    sequence = np.array([0.4,1.8,2.0])
    singleton_range = [-1.,2.]

    o1 = matrix_methods.vector_to_string([0.2,1.5,0.8],float)
    o2 = matrix_methods.vector_to_string([-0.5,1.1,0.4],float)
    o3 = matrix_methods.vector_to_string([1.,1.,0.5],float)
    o4 = matrix_methods.vector_to_string(sequence,float)
    optima_pr_map = defaultdict(float)
    optima_pr_map[o1] = 0.2
    optima_pr_map[o2] = 0.2
    optima_pr_map[o3] = 0.2
    optima_pr_map[o4] = 0.4

    dep_map = defaultdict(float)
    codep_map = defaultdict(float)

    sec = Sec(sequence,singleton_range,optima_pr_map,\
        dep_map,codep_map)
    return sec

# NOTE: invalid dep./codep. Pr. values. 
def Sec_sample_2():
    random.seed(3) 
    sequence = np.array([0.4,1.8,2])
    singleton_range = [-1.,2.]
    sec = Sec_sample_1()

    # add more points
    std_bound = np.array([deepcopy(singleton_range) \
        for _ in range(3)])
    ps = generate_pointgen_stdrand(std_bound,20,np.random)
    for p in ps:
        pstr = matrix_methods.vector_to_string(p,float) 
        sec.opm[pstr] = random.random()
    optima_pr_map = map_freq2pr(sec.opm)
    sec.opm = optima_pr_map

    # construct the dep. and co-dep. maps 
    bis = aprng_gauge.BatchIncrStruct(len(optima_pr_map),True,True,2)
    dep_map = defaultdict(float)
    codep_map = defaultdict(float)
    
    stat = True
    while stat: 
        q = next(bis) 
        stat = not type(q) == type(None) 
        if not stat: continue 
        vs = str(q[0]) + "," + "2." + str(q[1])

        is_in_dep = 1 if random.random() >= 0.5 else 0 

        if is_in_dep:
            dep_map[vs] = random.random()
        else:
            codep_map[vs] = random.random() 

    sec = Sec(sequence,singleton_range,optima_pr_map,\
        dep_map,codep_map)

    return sec 

def SecNet_graph_sample_C3():

    D = defaultdict(set)
    D[0] = {1,4}
    D[1] = {0,2,3,4}
    D[2] = {1,3}
    D[3] = {1,2,4}
    D[4] = {0,1,3,5,6}
    D[5] = {4,6}
    D[6] = {4,5}
    D[7] = {8,10}
    D[8] = {7}
    D[9] = {8}
    D[10] = {7}
    D[11] = {13}
    D[12] = {13}
    D[13] = {11,12,14,15}
    D[14] = {13,16}
    D[15] = {13}
    D[16] = {14,17}
    D[17] = {16}
    return D

def SecNet_graph_sample_CSmall():
    snv = []
    nsv = [0,1,2,3,4,5,6]
    sngs = SecNetGenScheme("spine frame",228)
    snfg = SecNetFrameGen(snv,nsv,sngs)

    snfg.construct_frame()
    D = defaultdict(set)
    for k,v in snfg.d.items():
        D[k] = set(v) 
    return D 

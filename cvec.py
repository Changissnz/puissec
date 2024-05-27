
from defaults import * 
from morebs2 import matrix_methods,rssi,measures

# these two functions are used to assign boolean
# values to float/int elements in 2-d matrices. 
metric_div = lambda x: np.less(x,1.0)
metrix_leq = lambda x: x == True 

# these functions are used to assign 1 boolean
# value for a boolean matrix

def metric_2dboolmat_to_bool(m,vx=0.65,row_wise=False):
    assert type(m) == np.ndarray
    assert str(m.dtype) in {'bool'}
    assert vx >= 0. and vx <= 1.0

    if m.size == 0: return False

    if not row_wise:
        value = np.round(m.sum() / m.size,5)
    else:
        value = 0
        for i in range(m.shape[0]): 
            q = (m[i].sum() / m.shape[1]) >= vx
            if q: value += 1
        value = value / m.shape[0] 

    return value >= vx 

########################################################

def top_percentile(v,percentile,is_reverse=False):
    assert matrix_methods.is_vector(v)
    assert type(percentile) == float and percentile >= 0.0 \
        and percentile <= 1.0

    l = int(round(percentile * len(v)))
    return sorted(v,reverse=is_reverse)[:l]

"""
converts a numerical matrix of shape 
m x n to a boolean matrix of the same 
shape using function `f` applied on 
each (i,j)'th element in m for the 
corresponding index of the output. 
"""
def std_float_matrix_to_bool_matrix(m,f):
    assert type(m) == np.ndarray
    assert str(m.dtype) in {'int32','float64'}

    x = np.zeros(m.shape,dtype=bool)

    for i in range(m.shape[0]): 
        for j in range(m.shape[1]):
            x[i,j] = f(m[i,j])
            assert type(x[i,j]) == np.bool_, "got {}".format(type(x[i,j]))
    return x 

# NOTE: unused.
def bool_matrix_to_bool(bm,f=metric_2dboolmat_to_bool):

    q = f(bm)
    assert type(q) == bool 
    return q 

class CVecISelector:

    """
    Data structure used specifically to output two 
    sequences of indices for a <CVec> to compare 
    between. 

    ref_range := [minimum int|float,maximum int|float]
    prev_range := [minimum int|float,maximum int|float]
    """
    def __init__(self,ref_range,prev_range):
        assert len(ref_range) == 2 
        assert len(prev_range) == len(ref_range)
        assert type(ref_range[0]) == type(ref_range[1])
        assert type(ref_range[0]) in {int,float}
        assert type(prev_range[0]) == type(prev_range[1])
        assert type(prev_range[0]) in {int,float}

        self.ref_range = ref_range
        self.prev_range = prev_range

    def output(self,cvec_length):
        prev,ref = None,None

        svp = self.subvec_by_type(cvec_length,"p")
        svr = self.subvec_by_type(cvec_length,"r")
        return svp,svr

    def subvec_by_type(self,cvec_length,stype="r"):
        assert stype in {"r","p"}
        
        def index_check(c):
            stat1 = True
            if c < 0:
                stat1 = abs(c) <= cvec_length
            else: 
                stat1 = c < cvec_length
            return stat1 
            
        
        prange = self.prev_range if stype == "p" else self.ref_range
        prev = None
        if type(prange[0]) == int:
            assert prange[0] <= prange[1]
            #assert max(prange) <= 0

            prev = [i for i in range(prange[0],\
                prange[1]) if index_check(i)]

            if len(prev) != prange[1] - prange[0]:
                return None
        else:
            assert prange[0] <= prange[1]
            assert prange[0] >= 0. and prange[1] <= 1.

            start_index = int(round(cvec_length * prange[0]))
            end_index = int(round(cvec_length * prange[1]))
            prev = [i for i in range(start_index,end_index)]
        return prev

def default_cvec_iselector_seq():
    cseq = [] 

    # b/t the first and second most recent
    ref_range = [-1,0]
    prev_range = [-2,-1] 
    cvis = CVecISelector(ref_range,prev_range)
    cseq.append(cvis)

    ref_range = [-3,0] 
    prev_range = [-6,-3]
    cvis = CVecISelector(ref_range,prev_range)
    cseq.append(cvis)

    ref_range = [-1,0] 
    prev_range = [0,1]
    cvis = CVecISelector(ref_range,prev_range)
    cseq.append(cvis)

    ref_range = [-2,0] 
    prev_range = [0,2]
    cvis = CVecISelector(ref_range,prev_range)
    cseq.append(cvis)

    ref_range = [0.8,1.]
    prev_range = [0.,0.2] 
    cvis = CVecISelector(ref_range,prev_range)
    cseq.append(cvis)

    ref_range = [0.8,1.]
    prev_range = [0.6,0.8] 
    cvis = CVecISelector(ref_range,prev_range)
    cseq.append(cvis)

    ref_range = [0.9,1.]
    prev_range = [0.75,0.9] 
    cvis = CVecISelector(ref_range,prev_range)
    cseq.append(cvis)
    return cseq

def basic_cvec_iselector_seq(): 
    ref_range = [-1,0]
    prev_range = [-2,-1] 
    cvis = CVecISelector(ref_range,prev_range)
    
    cseq = [] 
    cseq.append(cvis)
    return cseq 

class CVec:

    """
    a tool that records and analyzes performance 
    of a 3rd-party function F. 

    <CVec> records pairs of values 
    (input from F,performance output measure on the input). 

    At any index in the mutable sequence `v`; <CVec>
    can perform an analysis on it by outputting a 
    boolean vector with length equal to `cvis`. 
    Additionally, it can derive a lone boolean from
    the boolean vector.

    """
    def __init__(self,input_samples = [],start_vec=np.array([]),\
        cvis=[]):
        assert matrix_methods.is_vector(start_vec)
        assert len(start_vec) == len(input_samples)
        assert type(cvis) == list
        for c in cvis: assert type(c) == CVecISelector
        self.input_samples = [] 
        self.v = start_vec
        self.cvis = cvis 
        return

    def __len__(self):
        return len(self.input_samples) 

    def append(self,v,sample):
        assert type(v) in {int,float,np.float64,np.int32},"got {}".format(type(v))
        if len(self.input_samples) > 0:
            assert type(sample) == type(self.input_samples[-1]) 

        self.v = np.append(self.v,v)
        self.input_samples.append(sample) 

    def subvec(self,indices):
        sv = []
        for i in indices: 
            sv.append(self.v[i])
        return np.array(sv) 

    def add_CVIS(self,cvis):
        assert type(cvis) == CVecISelector
        self.cvis.append(cvis)

    def cmp(self,output_type=measures.zero_div):
        bvec = [] 
        for i in range(len(self.cvis)):
            ##print("comparator {}".format(i))
            m = self.cmp_by_selector(i,output_type)
            ##print("M")
            ##print(m)
            if type(m) == type(None):
                bvec.append(False) 
                continue 

            b = self.process_cmpmat(m)
            bvec.append(b)
            ##print("--------------")
        return np.array(bvec)

    """
    process comparison matrix
    """
    def process_cmpmat(self,m):
        assert type(m) == np.ndarray
        assert str(m.dtype) in {'bool','int32','float64'},"got type {}".format(m.dtype)

        if m.dtype != 'bool':
            ##print("convert to boolmat")
            m = std_float_matrix_to_bool_matrix(m,metric_div)
        ##print("operating on: {}".format(m))
        ##print()
        return metric_2dboolmat_to_bool(m)

    def cmp_by_selector(self,selector_index:int,output_type): 
        cvis = self.cvis[selector_index]
        prev,ref = cvis.output(len(self.v))
        if type(prev) == type(None) or type(ref) == type(None):
            return None 
        return self.cmp_indexseqs(prev,ref,output_type)

    """
    - description: 
    compares two index-sequences using the function 
    `output_type` (either division or leq). 

    - return:
    |iseq1| x |iseq2| matrix, each (i,j)'th value is 
    output from element iseq1[i] and iseq2[j]. 
    """
    def cmp_indexseqs(self,iseq1,iseq2,output_type):
        assert output_type in {measures.zero_div,np.less_equal}

        sv1 = self.subvec(iseq1)
        sv2 = self.subvec(iseq2) 
        output = np.zeros((len(iseq1),len(iseq2)))

        for i in range(len(iseq1)): 
            for j in range(len(iseq2)):                
                if output_type == np.less_equal:
                    q = output_type(sv2[j],sv1[i])
                else: 
                    q = output_type(sv2[j],sv1[i],1.1)
                ##print("--| ",q)
                output[i][j] = q
                
                ###

        return output

def cmp_seq_with_cvec(cvec,seq):
    assert type(cvec) == CVec
    input_samples = [0 for _ in range(len(seq))]

    cvec.input_samples = input_samples
    cvec.v = seq 

    bvec = cvec.cmp()
    return bvec 

############################ functions below 
############################ used for the k-multiplier search


def kmult_diffvec(V_m,V_f,k):
        return np.abs(V_m - V_f * k)

def CVec__scan_in_range(V_m,V_f,mrange,increment=0.1):
    cvec = CVec(cvis=basic_cvec_iselector_seq())

    q = float(mrange[0])
    min_ed = float('inf')
    best_sample = None
    stat = False
    while q < mrange[1] and not stat:
            V_ = V_f * q
            ed = matrix_methods.euclidean_point_distance(V_,V_m)
            ##print("q: {} d: {}".format(q,ed))
            if ed < min_ed: 
                    min_ed = ed
                    best_sample = q

            cvec.append(ed,q) 
            stat = cvec.cmp(output_type=np.less_equal)
            ##print("STAT: ",stat)
            q += increment
    ##print("DABEST: ",best_sample)
    return best_sample,min_ed 

def CVec__scan__kmult_search(V_m,V_f,depth=5):

        # get the k-vector
        kvec = np.array([measures.zero_div(v2,v1,np.nan) \
                for v2,v1 in zip(V_m,V_f)])
        kvec = np.unique(kvec)

        # sort the elements of kvec
        diffvec = [(k,np.sum(kmult_diffvec(V_m,V_f,k))) for \
                k in kvec]
        diffvec = sorted(diffvec,key=lambda x:x[1])
        ##print("diffvec: ",diffvec)
        
        # best sol'n, best case
        if len(diffvec) == 1: 
                return diffvec[0],0.0 

        dx = 1
        mrange = sorted([diffvec[0][0],diffvec[1][0]])
        ##print("MRANGE: ",mrange)
        bs,ed = None,float('inf')
        while dx <= depth:
                ##print("+ MRANGE: ",mrange)
                best_sample,best_ed = CVec__scan_in_range(\
                        V_m,V_f,mrange,increment=10**-dx)
                ##print("BEST SAMPLE: ",best_sample, best_ed)
                if type(best_sample) == type(None):
                    break 
                if best_ed < ed: 
                        bs,ed = best_sample,best_ed
                        mrange = [best_sample,best_sample + 10 **-dx] 
                else:
                        break

                dx += 1 
        return bs,ed

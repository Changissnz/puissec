from defaults import * 
from morebs2 import matrix_methods,rssi

def top_percentile(v,percentile,is_reverse=False):
    assert matrix_methods.is_vector(v)
    assert type(percentile) == float and percentile >= 0.0 \
        and percentile <= 1.0

    l = int(round(percentile * len(v)))
    return sorted(v,reverse=is_reverse)[:l]

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
        prange = self.prev_range if stype == "p" else self.ref_range
        prev = None
        if type(prange[0]) == int:
            assert prange[0] > self.prange[1]
            assert max(prange) <= 0

            prev = [i for i in range(prange[0],\
                prange[1])]
        else:
            assert prange[0] <= prange[1]
            assert prange[0] >= 0. and prange[1] <= 1.

            start_index = int(round(cvec_length * prange[0]))
            end_index = int(round(cvec_length * prange[1]))
            prev = [i for i in range(start_index,end_index)] 
        return prev 

class CVec:

    def __init__(self,start_vec=np.array([])):
        assert matrix_methods.is_vector(start_vec)
        self.v = start_vec
        self.cvis = []
        return

    def append(self,v):
        assert type(v) in {int,float}
        self.v = np.append(self.v,v)

    def subvec(self,indices):
        sv = []
        for i in indices: 
            sv.append(self.v[i])
        return np.array(sv) 

    def add_CVIS(self,cvis):
        assert type(cvis) == CVecISelector
        self.cvis.append(cvis)

    def cmp(self):
        return -1

    """
    - description: 
    compares two index-sequences using the function 
    `output_type` (either divide or geq). 

    - return:
    |iseq1| x |iseq2| matrix, each (i,j)'th value is 
    output from element iseq1[i] and iseq2[j]. 
    """
    def cmp_indexseqs(self,iseq1,iseq2,output_type):
        assert output_type in {np.divide,np.greater_equal}

        sv1 = self.subvec(iseq1)
        sv2 = self.subvec(iseq2) 

        output = np.zeros((len(iseq1),len(iseq2)))

        for i in range(len(iseq1)): 
            for j in range(len(iseq2)):
                output[i][j] = output_type(iseq1[i],iseq2[j]) 
        return output


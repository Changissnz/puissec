from defaults import * 
from morebs2 import matrix_methods,rssi

def top_percentile(v,percentile,is_reverse=False):
    assert matrix_methods.is_vector(v)
    assert type(percentile) == float and percentile >= 0.0 \
        and percentile <= 1.0

    l = int(round(percentile * len(v)))
    return sorted(v,reverse=is_reverse)[:l]

class CVec:

    def __init__(self,start_vec=np.array([])):
        assert matrix_methods.is_vector(start_vec)
        self.v = start_vec
        return

    def append(self,v):
        assert type(v) in {int,float}
        self.v = np.append(self.v,v)

    def subvec(self,indices):
        sv = []
        for i in indices: 
            sv.append(self.v[i])
        return np.array(sv) 

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


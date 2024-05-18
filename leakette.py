from iso_ring import *

########################### leak functions, uses 
########################### rnd_struct for decision-making

def choose_multiple(decimal,degree):
    mx = all_multiples_decimal(decimal)
    l = len(mx)
    q = min(int(round(l * degree,0)),l-1)
    return mx[q]

def idn_decimal(decimal):
    return decimal 

def subbound_for_decimal(decimal,degree,range2):
    assert range2[1] >= range2[0]

    dr = degree * (range2[1] - range2[0])
    
    s = max([decimal - dr,range2[0]])
    e = min([decimal + dr,range2[1]])
    return (s,e)

LEAKF_MAP = {0:choose_multiple,\
            1:idn_decimal,\
            2:subbound_for_decimal}

# TODO: future. 
def leakf__type_gcd(ir): 
    return -1 

"""
rnd_struct := .
degree := [0] fx := choose multiple -->
                (f1,f2); fx a float 
          [1] fx := idn_decimal -->
                f1; a float. 
          [2] fx := subbound_for_decimal -->
                (f1,(r1,r2)); f1,rx are floats,
                    and (r1,r2) is an ordered 
                    range.
fx := function
"""
def leakf__type_MV(ir:IsoRing,rnd_struct,degree,\
        fx):#degree2):
    assert fx in {choose_multiple,idn_decimal,\
                subbound_for_decimal}
    if fx == choose_multiple:
        assert len(degree) == 2
    elif fx == subbound_for_decimal:
        assert len(degree) == 3
        assert len(degree[2]) == 2 
    else: 
        assert degree >= 0.0 and degree <= 1.0 

    def one_f(i):
        s = seq[i]

        if fx == choose_multiple:
            return fx(s,degree[1])
        elif fx == idn_decimal:
            return fx(s)
        else: 
            return fx(s,degree[1],degree[2])

    def null_value():
        if fx == subbound_for_decimal:
            return np.array([np.nan,np.nan])
        return np.nan

    ti = ir.sec.seq_index()
    seq = ir.sec.optima_points()[ti]

    # choose the indices of the seq to leak. 
    indices = [i for i in range(len(seq))] 
    ls = len(seq)

    dx = None
    if fx == choose_multiple: 
        dx = degree[0]
    elif fx == idn_decimal: 
        dx = degree
    else: 
        dx = degree[0] 
    l = int(round(dx * ls))
    for _ in range(ls - l):
        r = rnd_struct.randrange(0,len(indices))
        indices.pop(r)

    v = []
    for i in range(ls):
        mx = None
        if i in indices:
            mx = one_f(i)
        else: 
            mx = null_value()
        v.append(mx)
    v = np.array(v)
    return v 


# TODO: test. 
"""
container structure holding the (function,value) 
pairs of information from a <Leak> instance onto 
an <IsoRing> w/ identity `ir_idn`. 
"""
class LeakInfo:

    def __init__(self,ir_idn):
        self.ir_idn = ir_idn

        # 0
        self.leak_info = {0:[],1:[],2:[]}

    @staticmethod
    def is_more_potent(self,p1,p2):
        stat1 = p1[0] < p2[0]
        stat2 = len(p1[1]) < len(p2[1])
        return stat1,stat2

    """
    return: 
    - 2-tuple, [0] is summation of potency per index
               [1] is subset of indices w/o leaks
    """
    def potency(self,b,sb):
        q = 0
        l = []
        for i in range(b.shape[0]):
            pl = self.process_leak_value_at_index(i,b,sb)
            if type(pl) != type(None):
                l.append(i)
            else:
                q += pl
        return q,l


    """
    i := index of the vector
    b := np.ndarray, bounds matrix
    sb := np.ndarray, sub-bound of b; the search info. 
    
    return: 
    - the minumum value from the 3 leak functions
      for index i. 
    """
    def process_leak_value_at_index(self,i,b,sb):
        qx0 = self.value_at_findex(0,i)
        qx1 = self.value_at_findex(1,i)
        qx2 = self.value_at_findex(2,i)

        h = [np.inf,np.inf,np.inf]
        if type(qx1) != type(None):
            h[0] = 1
        if type(qx0) != type(None):
            sbx = sb[i]
            h[1] = (sbx[1] - sbx[0]) / qx0
        if type(qx2) != type(None):
            sbx = sb[i]
            bx = b[i]
            q = measures.zero_div(\
                sbx[1] - sbx[0],bx[1] - bx[0],np.inf)
            h[2] = q
        return min(h) 

    def value_at_findex(self,f,i):
        vx = self.valuelist_at_findex(f,i)
        if len(vx) == 0:
            return None

        if i == 0:
            return max(vx)
        
        if i == 1:
            return vx[0]

        if i == 2:
            q = []
            for (i,x) in enumerate(vx):
                q.append((i,x[1] - x[0]))
            j = max(q,key=lambda q_:q_[1])[0]
            return vx[j] 

    def valuelist_at_findex(self,f,i):
        q = self.leak_info[f]
        v = float('inf')

        def not_null(rx): 
            if type(r) != np.ndarray:
                return not np.isnan(r) 
            return (not np.isnan(r[0]) and not np.isnan(r[1]))

        q2 = []
        for qx in q:
            r = qx[i]
            cn = not_null(r)
            if cn:
                q2.append(deepcopy(r)) 
        return q2

    def __add__(self,px):
        assert len(px) == 2
        assert px[0] in {0,1,2}
        assert type(px[1]) in {np.ndarray,list,tuple}
        self.leak_info[px[0]].append(px[1])
        return self

"""
procedure that leaks information about 
an <IsoRing>. To be used w/ intersection.
"""
class Leak: 

    """
    rnd_struct := a random structure, w/ at least 
                  the functions `randrange`.
    fd_seq := list, w/ elements 
            [0] leak* function
            [1] `degree` argument of the leak* function. 
    """
    def __init__(self,rnd_struct,fd_seq,leak_type="stationary"):
        assert len(fd_seq) > 0, "leak cannot be empty"

        for fd in fd_seq:
            assert len(fd) == 2
            assert fd[0] in {choose_multiple,\
                            idn_decimal,\
                            subbound_for_decimal}
            if fd[0] == choose_multiple:
                assert len(fd[1]) == 2
                assert fd[1][0] >= 0.0 and fd[1][0] <= 1.0
                assert fd[1][0] >= 0.0 and fd[1][0] <= 1.0
            elif fd[1] == idn_decimal:
                assert type(fd[1]) in {float,np.float64}
            else: 
                assert len(fd[1]) == 3
        assert leak_type in {"stationary","mobile"}

        self.rnd_struct = rnd_struct
        self.fd_seq = fd_seq 
        return

    @staticmethod
    def generate_Leak():
        return -1

    def leak_info(self,ir:IsoRing):
        i = ir.leak_stage
        if i >= len(self.fd_seq):
            return None

        x = self.fd_seq[i]
        fx = x[0]
        d = x[1] 
        outp = leakf__type_MV(ir,self.rnd_struct,d,fx)
        ir.leak_stage += 1 
        return outp

    def save_mem(self):

        return -1
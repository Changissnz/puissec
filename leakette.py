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

# TODO: future. 
def leakf__type_gcd(ir): 
    return -1 


"""
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
    def __init__(self,rnd_struct,fd_seq):
        assert len(fd_seq) > 0, "leak cannot be empty"

        for fd in fd_seq:
            assert len(fd) == 2
            assert fd[0] in {leakf__type_ivmap,\
                            leakf__type_multiple,\
                            leakf__type_subbound}
            if fd[0] == leakf__type_multiple:
                assert len(fd[1]) == 2
                assert fd[1][0] >= 0.0 and fd[1][0] <= 1.0
            else: 
                assert fd[1] >= 0.0 and fd[1] <= 1.0  

        self.rnd_struct = rnd_struct
        self.fd_seq = fd_seq 
        self.index = 0
        return

    @staticmethod
    def generate_Leak():
        return -1 

    def leak_info(self,ir:IsoRing):
        assert self.index < len(self.fd_seq)

        return -1 
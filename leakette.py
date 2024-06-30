from iso_ring import *

########################### leak functions, uses 
########################### rnd_struct for decision-making

def choose_multiple(decimal,degree):
    mx = all_multiples_decimal(decimal)
    l = len(mx)
    if l == 0:
        return 0.0 
        ##print("DEC: ",decimal)
    q = min(int(round(l * degree,0)),l-1)
    return mx[q]

def idn_decimal(decimal):
    return decimal 

def subbound_for_decimal(decimal,degree,range2):
    assert range2[1] >= range2[0]

    dr = degree * (range2[1] - range2[0])
    s = max([decimal - dr,range2[0]])
    e = min([decimal + dr,range2[1]])
    ##print("D,S,E:",dr,s,e)
    return (s,e)

# TODO: test this
def subbound_for_decimal_with_hop(decimal,degree,\
    range2,h:int,iter_eq=1):
    assert type(h) == int
    assert h > 0
    assert iter_eq >= 0 and iter_eq < h 

    sb = subbound_for_decimal(decimal,degree,range2)

    # calibrate so that the i'th increment in sb
    # equals the decimal
    lx = sb[1] - sb[0]
    q = lx / h * iter_eq
    ns = decimal - q 
    return (ns,ns + lx)

LEAKF_MAP = {0:choose_multiple,\
            1:idn_decimal,\
            2:subbound_for_decimal,\
            3:subbound_for_decimal_with_hop}

DEFAULT_LEAKF_LIST = list(LEAKF_MAP.values())

def leakf_to_index(leakf):

    if LEAKF_MAP[0] == leakf: 
        return 0
    if LEAKF_MAP[1] == leakf: 
        return 1
    if LEAKF_MAP[2] == leakf: 
        return 2
    if LEAKF_MAP[3] == leakf: 
        return 3
    return -1

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
                (f1,f2,(r1,r2)); fx,rx are floats,
                    and (r1,r2) is an ordered 
                    range.
          [3] fx := subbound_for_decimal_with_hop -->
                (f1,f2,(r1,r2),h); fx,rx are floats,
                    and (r1,r2) is an ordered 
                    range. Value h is the hop size 

fx := function
"""
def leakf__type_MV(ir:IsoRing,rnd_struct,degree,\
        fx):
    assert fx in {choose_multiple,idn_decimal,\
                subbound_for_decimal,\
                subbound_for_decimal_with_hop}

    if fx == choose_multiple:
        assert len(degree) == 2
    elif fx == subbound_for_decimal:
        assert len(degree) == 3
        assert len(degree[2]) == 2 
    elif fx == idn_decimal: 
        assert degree >= 0.0 and degree <= 1.0 
    else:
        assert len(degree) == 4
        assert len(degree[2]) == 2 

    def one_f(i):
        s = seq[i]
        if fx == choose_multiple:
            return fx(s,degree[1])
        elif fx == idn_decimal:
            return fx(s)
        elif fx == subbound_for_decimal: 
            return fx(s,degree[1],degree[2])
        else: 
            return fx(s,degree[1],degree[2],degree[3])

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
    #v = np.array(v)
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
        self.leak_info = {0:[],1:[],2:[],3:[]}

    def __str__(self):
        s = "LEAK_INFO:{}\n\n".format(self.ir_idn)
        for (k,v) in self.leak_info.items():
            s += str(k) + "\n" + str(v) + "\n\n"
        return s

    ##########################################
    def __add__(self,px):
        assert len(px) == 2
        assert px[0] in {0,1,2,3}
        assert type(px[1]) in {np.ndarray,list,tuple}
        self.leak_info[px[0]].append(px[1])
        return self

"""
container to hold a <Leak> instance's extracted information
for each <IsoRing>+<Leak> interaction.
"""
class LeakMem:

    def __init__(self):
        # sec idn -> LeakInfo container 
        self.d = {}
        return

    def __str__(self):
        s = "LEAKMEM" + "\n\n\t"
        s += str(self.d)
        return s 

    """

    ir_idn := [0] sec idn tag, 
              [1] function index, 
              [2] leak output
    """
    def __add__(self,ir_idn):
        assert len(ir_idn) == 3
        if ir_idn[0] not in self.d:
            self.d[ir_idn[0]] = LeakInfo(ir_idn[0])
        self.d[ir_idn[0]] = self.d[ir_idn[0]] + (ir_idn[1],ir_idn[2]) 
        return self

    def fetch_LeakInfo(self,sec_idn,f_idn):

        if sec_idn not in self.d:
            return None

        if f_idn not in self.d[sec_idn]:
            return None

        return self.d[sec_idn][f_idn]

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
            assert fd[0] in DEFAULT_LEAKF_LIST
            
            if fd[0] == choose_multiple:
                assert len(fd[1]) == 2
                assert fd[1][0] >= 0.0 and fd[1][0] <= 1.0
                assert fd[1][0] >= 0.0 and fd[1][0] <= 1.0
            elif fd[0] == idn_decimal:
                assert type(fd[1]) in {float,np.float64}
            elif fd[0] == subbound_for_decimal: 
                assert len(fd[1]) == 3
            else: 
                assert len(fd[1]) == 4
        assert leak_type in {"stationary","mobile"}

        self.rnd_struct = rnd_struct
        self.fd_seq = fd_seq 
        self.leakm = LeakMem()
        # (sec idn, F idn,output)
        self.prev_li = None
        return

    """
    pmap := dict, length is reference dimension d
    """
    @staticmethod
    def generate_Leak__type1(num_leakf,rnd_struct):
        leakfs = [] 
        for i in range(num_leakf):
            q = rnd_struct.randrange(0,4) 
            leakf = Leak.generate_leakf_args(q,rnd_struct)
            leakfs.append((LEAKF_MAP[q],leakf)) 
        return Leak(rnd_struct,leakfs)
    
    @staticmethod
    def generate_leakf_args(idn,rnd_struct):
        assert idn in {0,1,2,3}

        if idn == 0: 
            d1 = rnd_struct.uniform(0.,1.)
            d2 = rnd_struct.uniform(d1+1e-10,1.) 
            return np.array((d1,d2)) 
        if idn == 1:
            return rnd_struct.uniform(0.,1.)

        d1 = rnd_struct.uniform(0.,1.)
        d2 = rnd_struct.uniform(0.,1.)
        d3 = (0.,1.)

        if idn == 2: 
            return (d1,d2,d3)

        d4 = rnd_struct.randrange(2,13) 
        return (d1,d2,d3,d4)

    def leak_info(self,ir:IsoRing):
        i = ir.leak_stage
        if i >= len(self.fd_seq):
            self.prev_li = None
            return None

        x = self.fd_seq[i]
        fx = x[0]
        d = x[1] 
        outp = leakf__type_MV(ir,self.rnd_struct,d,fx)
        ir.leak_stage += 1 

        self.save_mem(ir,fx,outp)
        return outp

    def save_mem(self,ir,fx,outp):
        it = ir.sec.idn_tag
        fi = leakf_to_index(fx)
        self.leakm = self.leakm + (it,fi,outp)
        self.prev_li = (it,fi,outp)
        return

def Leak_sample1():
    l1 = (LEAKF_MAP[0],np.array((0.5,0.5)))
    l2 = (LEAKF_MAP[2],(0.2,0.4,(0.0,1.0)))
    l3 = (LEAKF_MAP[0],np.array((0.5,0.5)))

    B = np.array([[0,1.],\
            [0.,1],\
            [0,1.],\
            [0.,1],\
            [0.,1.]])

    SB = np.array([[0,0.5],\
            [0,0.5],\
            [0.,0.5],\
            [0,0.5],\
            [0.,1.]])

    L = [l1,l2,l3]

    random.seed(332)
    return Leak(random,L)

def Leak_sample2(): 

    l1 = (LEAKF_MAP[1],1.0)

    B = np.array([[0,1.],\
                    [0.,1],\
                    [0,1.],\
                    [0.,1],\
                    [0.,1.]])

    SB = np.array([[0,0.5],\
            [0,0.5],\
            [0.,0.5],\
            [0,0.5],\
            [0.,1.]])

    ir = IsoRing_sample_1()
    ir.sec.idn_tag = 12 
    return Leak(random,[l1])

from crackling import * 

class CBridge:

    def __init__(self,crackling,isoring,hs,ssih=5,\
        cidn=None,batch_size=1000,verbose=False):
        assert type(crackling) == Crackling
        assert type(isoring) == IsoRing 
        assert type(hs) == HypStruct

        self.crackling = crackling 
        self.isoring = isoring
        self.hs = hs 
        self.rssi = None 
        self.verbose = verbose 

        self.load_crackf()#ssih)
        self.cidn = cidn
        self.bs = batch_size
        self.batch = None

        self.load_rssi_batch() 
        return

    def load_crackf(self):#,h=5):
        print("LOADING ",self.verbose)
        self.rssi = default_base_RSSI(self.isoring,self.crackling,\
            self.hs,self.verbose)

    def load_rssi_batch(self):
        self.batch = rssi.ResplattingSearchSpaceIterator.\
            iterate_one_batch(self.rssi,self.bs)

    """
    return:
    - (<Crackling> idn,<IsoRing> idn)
    """
    def agent_idns(self):
        return (self.crackling.cidn,self.isoring.sec.idn_tag)

    def __next__(self):
        if self.rssi.terminated: 
            print("DONE")
            return 
        try:
            p = next(self.batch)
        except: 
            print("LOAD BATCH")
            self.load_rssi_batch()
            return 
        

        #p = next(self.batch)
        if type(p) == None: 
            print("FINISHED!!")
            self.load_rssi_batch()
            p = next(self.batch)
            if type(p) == type(None):
                return None

        if self.verbose:
            print('next point on bridge {}:\n{}'.format(\
                self.agent_idns(),p))
        return p

################ 

class IRCLeakDict:

    def __init__(self,num_leakf_range,rnd_struct):
        assert len(num_leakf_range) == 2
        assert num_leakf_range[0] < num_leakf_range[1] and \
            num_leakf_range[0] >= 0  

        self.num_leakf_range = num_leakf_range 
        self.rnd_struct = rnd_struct 
        # sec idn -> opt. dim -> <Leak>
        self.d = defaultdict(defaultdict) 
        return

    """
    """
    def load_default_Leak_for_IsoRing(self,ir:IsoRing):
        q = ir.secdim_seq()
        self.d[ir.sec.idn_tag] = defaultdict(Leak) 
        for q_ in q:
            nl = self.rnd_struct.randrange(self.num_leakf_range[0],\
                self.num_leakf_range[1])
            lx = Leak.generate_Leak__type1(nl,self.rnd_struct)
            self.d[ir.sec.idn_tag][q_] = lx
        return

    def fetch_Leak(self,sec_idn,opt_dim):
        if sec_idn not in self.d: return None
        if opt_dim not in self.d[sec_idn]: return None
        return self.d[sec_idn][opt_dim] 


###################################################################
        
def one_correct_HypStruct_for_IsoRing(ir:IsoRing):
    ti = ir.sec.seq_index()
    seq = ir.sec.optima_points()[ti]
    b = np.array([seq,seq + 10 ** -4]).T
    sbs = [b] 
    sb_pr = np.array([1.0])
    return HypStruct(ir.sec.idn_tag,ti,sbs,sb_pr)

def one_dummy_HypStruct_for_IsoRing(ir:IsoRing):
    ti = ir.sec.seq_index()

    seq = ir.sec.optima_points()[ti]
    b = np.zeros((len(seq),2),dtype=float)
    b[:,1] = 1.0
    sbs = [b]
    sb_pr = np.array([1.0])
    return HypStruct(ir.sec.idn_tag,ti,sbs,sb_pr)


def round_ratio_to_closest_hop(r,h,round_depth=5):
    assert r >= 0.0 and r <= 1.0
    assert h > 0

    q = 1.0 / h 
    x = math.floor(r / q) 
    if type(round_depth) != type(None): 
        return round(x * q,5)
    return x * q

def round_length_to_closest_hop(l,h,round_depth=5):
    assert h > 0
    q = float(l) / h

    if type(round_depth) != type(None):
        return round(q * h,5)
    return q * h 

"""
EX: h = 3.
    b = [0,1]*4
    r = <0.33,0.66, 1.0>

"""
def calibrate_boundloc_values(bound_length,location_ratio,hop_size,lx):
    if type(bound_length) == float:
        bound_length = np.ones((lx,)) * bound_length
    assert len(bound_length) == lx
    bound_length = np.array([round_length_to_closest_hop(bl,hop_size) \
        for bl in bound_length]) 

    if type(location_ratio) == float:
        location_ratio = np.ones((lx,)) * location_ratio
    assert len(location_ratio) == lx
    location_ratio = np.array([round_ratio_to_closest_hop(lr,hop_size) \
        for lr in location_ratio]) 
    return bound_length,location_ratio 

"""
For the <IsoRing> `ir`, declares a <HypStruct> H on 
its `seci`'th <Sec> S. H has one bound B, with each 
j'th dimension of length equal to `bound_length` 
(type float) or `bound_length[j]` (type np.ndarray). 

The `location_ratio` specifies the location of S.seq 
in B. And the hop 

"""
def one_approximate_HypStruct_for_IsoRing(ir:IsoRing,seci:int,\
    bound_length,location_ratio,hop_size):
    assert type(bound_length) in {float,np.ndarray} 
    assert type(location_ratio) in {float,np.ndarray} 

    # make the bounds
    ir.set_isorep(seci) 
    lx = ir.sec.dim()

    bound_length,location_ratio = calibrate_boundloc_values(\
        bound_length,location_ratio,hop_size,lx) 

    v = ir.sec.seq
    v_ = []
    bxs = []
    for (i,v2) in enumerate(v):
        qx = location_ratio[i] % 1.0
        blx = bound_length[i] * qx
        npx1 = round(v2 - blx,5)
        npx2 = npx1 + bound_length[i]
     
        bxs.append([npx1,npx2]) 
    bxs = np.array(bxs) 

    # declare the HypStruct
    ti = ir.sec.seq_index()

    # case: make only one. 
    hsx = HypStruct(ir.sec.idn_tag,ti,[bxs],\
        sb_pr=np.array([1.0]),\
        hs_vec=np.array([hop_size]))

    return hsx 






 






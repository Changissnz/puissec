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

        self.load_crackf(ssih)
        self.cidn = cidn
        self.bs = batch_size
        self.batch = None

        self.load_rssi_batch() 
        return

    def load_crackf(self,h=5):
        print("LOADING")
        self.rssi = default_base_RSSI(self.isoring,self.crackling,\
            self.hs,h,self.verbose)

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
        
        try:
            p = next(self.batch)
        except: 
            print("FINISHED!")
            return 
        
        #p = next(self.batch)

        if type(p) == None: 
            print("FINISHED!!")

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
from leakette import * 
from morebs2 import relevance_functions,hop_pattern
from hype import * 
from tdir import * 

class Crackling:

    def __init__(self,cmp_deg=1.5,cidn=None,\
        cvsz=float('inf')):
        self.cmp_deg = cmp_deg 
        self.cidn = cidn
        self.hs = None
        # <HypStruct> instance from interdiction
        # w/ `hs`+<LeakInfo>.
        self.hsi = None
        self.td = None
        self.rss = None
        self.set_cvec(cvsz)
        # ?target has been cracked?
        self.astat = False
        self.cracked_dict = defaultdict(float)
        self.flagged_pts = []
        # target has been ?captured?
        self.cstat = False
        # target has been ?interdicted?
        self.istat = False 
        # ?<TDirector> path is coordinated by <Cracker>?
        self.is_coordinated = False

    ###################### functions to get/set <HypStruct>,<TDirector> variables

    def load_HypStruct(self,hs):
        assert type(hs) == HypStruct
        self.hs = hs 

    def target_of_HypStruct(self):
        if type(self.hs) == type(None): 
            return None
        return self.hs.seq_idn

    def load_TDirector(self,td): 
        assert type(td) == TDirector
        self.td = td
        return

    def fetch_tdir(self):
        if type(self.td) == type(None):
            return None
        return self.td.td 

    def loc(self):
        if type(self.td) == type(None):
            return None
        return self.td.loc()  

    def set_cvec(self,cvsz):
        ciseq = default_cvec_iselector_seq()
        cvec = CVec(cvis=ciseq,sz_limit=cvsz)
        self.cvec = cvec 

    ###################################### navigation methods

    """
    """
    def td_next(self,timespan=1.0,verbose:bool=False):

        if verbose:
            print("--> TD-NEXT FOR C={},N={}".format(self.cidn,\
                self.td.loc()))
            print("\tOBJSTAT: ",self.td.obj_stat)

        if self.td.obj_stat == "search for target":
            return self.default_TD__search_for_target(timespan,verbose)
        return self.default_TD__capture_target(timespan,verbose)

        #############################################

    # TODO: unused. 
    def default_TD__search_for_target(self,timespan,verbose:bool=False):
        assert self.td.obj_stat == "search for target"

        q = self.td.check_radar()

        if len(q) > 0:
            if verbose: 
                print("switching from [search]->[capture] target.")
            self.td.switch_obj_stat()
            return self.default_TD__capture_target() 

        # check if <TDir> is still active.
        if not self.td.td.active_stat:
            if verbose: 
                print("NOT ACTIVE,SETTING EXTREME")
            self.td.extloc_search__set_TDir(extf=max,rnd_struct=random)
        
        if verbose: 
            print("path is ")
            print(str(self.td.td.node_path))
            print("------")

        v = len(self.td.td.node_path) - 1 
        v = int(round(v / timespan)) 
        self.td.td.velocity = v
        self.td.td.scaled__next__(timespan)
        return v

    # TODO: unused 
    def default_TD__capture_target(self,timespan,verbose:bool):
        assert self.td.obj_stat == "capture target" 

        q = self.td.check_radar()
        if verbose:
            print("RADAR: ",len(q))

        # CASE: no <Crackling>s in sight, switching 
        # objective to "search for target"
        if len(q) == 0:
            if verbose: 
                print("switching from [capture]->[search] target.")
            self.td.switch_obj_stat()
            return self.default_TD__search_for_target(timespan,verbose)

        dpd = self.td.default_crackling_pathdec(\
            predicted_distance=None,rnd_struct=random)        

        # set the velocity equal to the path length
        self.td.load_new_path(dpd)
        v = len(dpd) -1 
        self.td.td.velocity = int(round(v / timespan)) 
        if verbose: 
            print("PATH")
            print(self.td.td.node_path)
            print("============")

        l = self.td.loc()
        self.td.td.scaled__next__(timespan)
        l2 = self.td.loc()
        if verbose: 
            print("-- TRAVEL {}->{}".format(l,l2))
        return v

    """
    used by <Crackling> instance to decide 
    to accept an entry point as a point based 
    on radar status concerning target <IsoRing>.
    """
    def accept_TDirector_at_entry_point(self):
        assert type(self.td) == TDirector
        q = self.td.check_radar()
        return len(q) > 0

    ################################## feedback-response functions 

    """
    p := vector, the attempted point
    q := vector, length equal to number of optima,
            each i'th value is the score of p w.r.t.
            the i'th optimum point. 
    astat := bool, ?is `p` equal to any local optima 
                   of the targeted <IsoRing> instance.
    """
    def register_response(self,p,q,astat:bool): 
        assert type(self.hs) != type(None)
        if not self.astat:
            self.astat = astat 

        ##print("PP: {}\nQQ: {}".format(p,q))
        x = q[self.hs.target_index] 
        self.cvec.append(x,p)
        
        qx1 = self.cvec.cmp(measures.zero_div)
        qx2 = self.cvec.cmp(np.less_equal) 
        ##print("QX1: ",qx1)
        ##print("QX2: ",qx2)

        ##s = [qx1[self.hs.target_index],qx2[self.hs.target_index]] 
        
        ## arg #1
        s = [np.sum(qx1) >= len(qx1) / self.cmp_deg,\
            np.sum(qx2) >= len(qx2) / self.cmp_deg]
        
        d = s[0] or s[1] 
        """
        d = metric_2dboolmat_to_bool(np.array([qx1,qx2]),\
            0.4,True)
        """
        if d: 
            self.flagged_pts.append(np.copy(self.cvec.input_samples[-1]))#np.copy(self.cvec[-1]))
            #self.flagged_pts.append(len(self.cvec) - 1)
        return d

    def register_lo_pr(self,prx):
        assert len(self.cvec.input_samples) > 0
        v = self.cvec.input_samples[-1]
        s = matrix_methods.vector_to_string(v,float)
        self.cracked_dict[s] = prx 
        return

    ################# methods for <Cracker> coordination

    def assign_coordinated_path(self,p):
        assert type(p) == NodePath
        self.td.td.node_path = mp_
        self.td.td.index = 0
        self.is_coordinated = True
        return

    def recv_open_info(self,open_info_type,idn,info):
        assert open_info_type in {1,2}
        self.td.td.open_info_var.append((open_info_type,idn,info)) 
        return

"""
- return: 
a <RChainHead> instance derived from a 
<IsoRing>. 
"""
def IsoRing_and_Crackling_to_base_RChainHead(ir:IsoRing,cracklng:Crackling,verbose):
    assert type(ir) == IsoRing

    def outputf1(p):
        if verbose:
            print("IR registers attempt")
            print(p)

        # NOTE: uncomment in real simulations.
        """
        if cracklng.astat:
            print("TERMINATED")
            return 
        """
        q,stat = ir.register_attempt(p)
        ##
        if verbose: 
            print("-- register")
            print("score: ",q[cracklng.hs.target_index])
            print("stat: ",stat)
        
        ##
        """
        if stat:
            return None
        """

        if verbose:
            print("cracking resp: ")
        
        d = cracklng.register_response(p,q,stat)
        if stat:
            prx = ir.response_to_prompt(cracklng.hs.target_index)
            if verbose: 
                print("pr. output: ",prx)
            cracklng.register_lo_pr(prx) 
        if verbose: 
            print("response: ",d)
            print("----------------------")

        return d

    rch = relevance_functions.RChainHead()

    argseq1 = ['nr',outputf1]
    rch.add_node_at(argseq1)
    return rch

"""
"""
def default_base_RSSI(ir:IsoRing,cracklng:Crackling,\
    hs:HypStruct,verbose=False): 

    #ssih,\verbose=False):
    assert type(ir) == IsoRing
    assert type(hs) == HypStruct

    ##print("\t-- converting to RCH")
    rch = IsoRing_and_Crackling_to_base_RChainHead(ir,cracklng,verbose)
    ##print("\t-- CONVERTED")
    resplattingMode = ("relevance zoom",rch)

    ix = hs.most_probable_subbound_i()
    mpsb = deepcopy(hs.suspected_subbounds[ix])
    ##mpsb = hs.most_probable_subbound()
    ssih = hs.hs_vec[ix] 
    if verbose: 
        print("[X] declarationes de RSSI aufbund")
        print(mpsb)
        print("-------------------------------------")

    print("[X] declarationes de RSSI aufbund")
    print(mpsb)
    print("-------------------------------------")


    start_point = deepcopy(mpsb[:,0])
    rss = rssi.ResplattingSearchSpaceIterator(mpsb,\
        start_point,SSIHop=ssih,resplattingMode = \
            resplattingMode)
    ##print("DECLARED RSS")
    return rss

"""
class w/ no variables; contains basic functions to 
improve a <HypStruct> instance. 
"""
class HypInfer: 

    def __init__(self):
        return

    @staticmethod
    def infer_by_LeakInfo(hypStruct,leak_info,inference_type=1):
        assert type(hypStruct) == HypStruct
        assert type(leak_info) == LeakInfo
        assert inference_type in {1,2}

        ###
        '''
        print("LEAK INFO")
        print(leak_info.leak_info)
        print("-------------------")
        '''
        ###

        ## NOTE: model 1, 
        ## for each leak, process all leaked values from 
        ## it. Order of processing is by map-key order.
        if inference_type == 1:
            for lk1,lk2 in leak_info.leak_info.items():
                for lk2_ in lk2:
                    hypStruct = HypInfer.infer_FX(\
                        hypStruct,lk2_,lk1)
            return hypStruct

        ## NOTE: model 2,
        ## uses a vector of the best calculated leak values for
        ## each index of S.V 
        liv = HypInfer.best_LeakInfo_valuevec(hypStruct,leak_info)
        hypStruct = HypInfer.apply_delta_on_HypStruct_by_best_valuevec(\
            hypStruct,liv)
        return hypStruct

    ############################ model 2

    @staticmethod
    def apply_delta_on_HypStruct_by_best_valuevec(hypStruct,bvv):
        for bv in bvv:
            hs  = HypInfer.infer_fx(hypStruct,bv[1],bv[0])
        return hs


    @staticmethod
    def best_LeakInfo_valuevec(hypStruct,leak_info):
        assert type(leak_info) == LeakInfo
        d = hypStruct.dim() 
        # every elmnt is (leak function identifier,leak value)
        vx = [] 
        for i in range(len(d)):
            qr = leak_info.best_value_at_index(i)
            vx.append(qr) 
        return vx

    #################################### model 1
    @staticmethod
    def infer_FX(hypStruct,leak_value,leak_idn):
        assert leak_idn in {0,1,2,3}

        if leak_idn != 3:
            hs = HypInfer.alter_all_delta(hypStruct,leak_idn,leak_value)
        else:
            hs = HypInfer.replace_all_delta(hypStruct,leak_idn,leak_value)
        return hs 


    ##############################################################

    """
    leak_idn := int, one in [0...3]
    sx := np.ndarray, proper bounds. 
    leak_value := ?, leak value corresponding to `leak_idn`
    """
    @staticmethod
    def exec_fn(leak_idn,sx,leak_value):
        if leak_idn == 0:
            return adjust_bounds__F0(sx,leak_value)
        elif leak_idn == 1:
            return adjust_bounds__F1(sx,leak_value)
        elif leak_idn == 2: 
            return adjust_bounds__F2(sx,leak_value)
        else: 
            return adjust_bounds__F3(sx,leak_value)

    """
    used for leak_idn in [0,1,2]
    """
    @staticmethod
    def alter_all_delta(hypStruct,leak_idn,leak_value):
        assert leak_value != 3
        for (i,sb) in enumerate(hypStruct.suspected_subbounds): 
            outp = HypInfer.exec_fn(leak_idn,sb,leak_value)
            hypStruct.suspected_subbounds[i] = outp 
        return hypStruct 

    """
    used for leak_idn := 3
    """
    @staticmethod
    def replace_all_delta(hypStruct,leak_idn,leak_value):

        hs = hypStruct.new_HypStruct_by_next_subbound()

        outp = HypInfer.exec_fn(leak_idn,hs.suspected_subbounds[0],leak_value)

        hs.suspected_subbounds = outp[0]
        hs.hs_vec = [outp[1] for _ in range(len(hs.suspected_subbounds))]
        hs.sb_pr = [1.0/len(hs.suspected_subbounds) for _ \
            in range(len(hs.suspected_subbounds))]
        return hs

def closest_reference_to_bound_start(b,V_f):

    null_vec = np.zeros((len(V_f,)))
    start_vec = deepcopy(null_vec) 

    for (i,v) in enumerate(V_f):
        if not np.isnan(v): 
            qx = measures.zero_div(b[i,0],v,np.nan)

            if np.isnan(qx):
                null_vec[i] = b[i,0]
                continue

            xq = int(round(qx))
            start_vec[i] = xq * v
        else: 
            null_vec[i] = b[i,0]

    return start_vec + null_vec 

def adjust_bounds__F0(b,V_f,adjust_increment):
    qi = [vf if (not np.isnan(vf)) else 0.0 for vf in V_f] 
    
    nb = deepcopy(b)

    for (i,vf) in enumerate(V_f):
        if np.any(np.isnan(vf)): continue
        nb[i] = [vf,vf * 5] 
    return nb 

def adjust_bounds__F1(b,V_f):
    bx = deepcopy(b)

    for (i,bx_) in enumerate(bx):
        if not np.isnan(V_f[i]):
            rx = np.array([V_f[i],V_f[i] + 0.05]).T
            bx[i] = rx
    return bx

def adjust_bounds__F2(b,V_f):
    ##assert matrix_methods.is_proper_bounds_vector(V_f) 
    ##print("V_F: ",V_f)
    V_f = np.array(V_f)
    assert V_f.shape[1] == 2 and V_f.shape[0] > 0 
    b2 = deepcopy(b)

    for (i,bx) in enumerate(V_f):
        stat1 = np.isnan(bx[0]) 
        stat2 = np.isnan(bx[1])
        stat = stat1 or stat2

        if stat: continue
        b2[i] = deepcopy(bx)
    return b2 

def adjust_range__F3(r,h_prev,h_ref):
    assert len(r) == 2 and r[0] <= r[1]
    hp_prev = hop_pattern.HopPattern(r[0], r[0], r[1], cycleLog = True, DIR = 1.0/h_prev)
    prev_set = set(hop_pattern.cycle_hop_pattern(hp_prev))

    hp_ref = hop_pattern.HopPattern(r[0], r[0], r[1], cycleLog = True, DIR = 1.0/h_ref)
    ref_set = set(hop_pattern.cycle_hop_pattern(hp_ref))

    qx = ref_set.intersection(prev_set)
    prev_set = prev_set - qx
    range_seq = [deepcopy(r)]

    while len(prev_set) > 0:
        mprev = min(prev)
        r1 = np.array([mprev,mprev + (r[1] - r[0])])

        hp_ref = hop_pattern.HopPattern(r1[0], r1[0], r1[1], cycleLog = True, DIR = 1.0/h_ref)
        ref_set = set(hop_pattern.cycle_hop_pattern(hp_ref))

        qx = ref_set.intersection(prev_set)
        prev_set = prev_set - qx
        range_seq.append(r1) 

    return range_seq 

# NOTE: caution, different hop sizes not accounted for.
def adjust_bounds__F3(b,V_f):
    ##print("V_F: ",V_f) 
    b2 = deepcopy(b)
    mlist = []

    hs_counter = Counter() 
    for (i,bx) in enumerate(V_f):
        stat = np.any(np.isnan(bx))
        if stat: continue
        hs_counter[bx[2]] += 1

    if len(hs_counter) == 0:
        return None 

    m = max([(k,v) for k,v in hs_counter.items()],key=lambda x:x[1])
    fmax = m[0] 
    hs_counter.pop(m[0]) 

    b2 = [[b__] for b__ in b]
    for (i,bx) in enumerate(V_f):
        stat = np.any(np.isnan(bx))
        if stat: continue

        if bx[2] == fmax:
            b2[i] = [np.array([bx[0],bx[1]])]
        else:
            rxs = adjust_range__F3([bx[0],bx[1]],bx[2],fmax) 
            b2[i] = rxs

    # select all bounds by permutation
    b2i = [len(b2_) for b2_ in b2] 
    ivp = IndexVecPermuter(b2i)
    all_bounds = []
    while not ivp.finished:
        ix = next(ivp)
        bx = np.array([deepcopy(b2[j][ix_]) for (j,ix_) in enumerate(ix)])
        all_bounds.append(bx)
    return all_bounds,fmax 


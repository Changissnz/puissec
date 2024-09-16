from crackling import * 
from secnet import * 
from leakette import *
from nERg import *  

"""
splits a bound `bs` by a constant 
size c corresponding to the splitting
size `sz`. 
"""
def even_bound_split(bs,sz):
    assert matrix_methods.is_proper_bounds_vector(bs)
    assert sz >= 1 and type(sz) == int

    r = (bs[:,1] - bs[:,0]) / sz 
    bsx = []
    q = deepcopy(bs[:,0])
    for i in range(sz):
        bs_ = q + r
        bsq = deepcopy(np.array([q,bs_]).T)
        bsx.append(bsq)
        q = bs_
    return bsx

DEFAULT_HYPSTRUCT_GENSCHEMES = {0:one_dummy_HypStruct_for_IsoRing,
                1: one_approximate_HypStruct_for_IsoRing,\
                2: one_correct_HypStruct_for_IsoRing}

def generated_HypStruct__int_idn(fx,is_forward:bool=True):
    if is_forward: 
        assert fx in DEFAULT_HYPSTRUCT_GENSCHEMES
        return DEFAULT_HYPSTRUCT_GENSCHEMES[fx]
    
    dh = std_invert_map(DEFAULT_HYPSTRUCT_GENSCHEMES)
    assert fx in dh
    return dh[fx] 
     


class OrderOfCrackn:

    def __init__(self):
        return

    def order_by_depchain_map(self,dcm):
        print("DCMM: ",dcm)
        l = [(k,len(v)) for k,v in dcm.items()]
        l = sorted(l,key=lambda x:x[1])[::-1]
        ks = [l_[0] for l_ in l]
        soln = [{l_[0]} for l_ in l]

        def check_coincidence(k1,k2):
            assert k1 in dcm
            assert k2 in dcm

            v1 = dcm[k1]
            v2 = dcm[k2]
            vx1,vx2 = set(),set()
            ##print("V1: ",v1)
            ##print("V2: ",v2)
            for v1_ in v1:
                vx1 = vx1 | v1_
            for v2_ in v2:
                vx2 = vx2 | v2_

            return len(vx1.intersection(vx2)) > 0

        """
        if depchains do not coincide, merge. 
        """
        def check_set_coincidence(i1,i2):
            s1,s2 = soln[i1],soln[i2]
            for s1_ in s1:
                for s2_ in s2:
                    if check_coincidence(s1_,s2_):
                        return False
            return True

        def sort_swap(key):

            vx = dcm[key][::-1]
            vx_ = []
            ##print("VX: ",vx)
            for vx2 in vx:
                vx_.extend(list(vx2))

            # get the max index of the 
            # dependencies.
            indexio = {}
            for vx2 in vx_:

                for (i,sx) in enumerate(soln):
                    if vx2 in sx:
                        indexio[vx2] = i
                        break

            index = max(indexio.values())

            keyloc = -1
            for (i,sx) in enumerate(soln):
                if key in sx:
                    keyloc = i
                    break
            assert keyloc != -1

            # case: keyloc >= index; do nothing
            if keyloc >= index:
                return

            # case: keyloc < index;
            q = soln.pop(keyloc)
            soln.insert(index,q)

            #   check the coincidence w/ neighbor
            if index + 1 >= len(soln):
                return

            # no coincidences; merge. 
            if check_set_coincidence(index,index+1): 
                s1,s2 = soln[index],soln[index+1]
                s1 = s1 | s2
                soln.pop(index)
                soln.pop(index)
                soln.insert(index,s1) 

            return

        # iterate through each key and locate index
        for ks_ in ks:
            sort_swap(ks_) 
        return soln

class BackgroundInfo:

    """
    """
    def __init__(self,opm,depchain_map,codep_sets,\
        dec_map,leak_map):
        # sec idn -> opt. dim. -> opt. idn -> Pr. 
        self.opm = opm
        # sec idn -> (opt index,other sec idn, other opt index)
        self.dm = depchain_map
        self.cdm = codep_sets
        # sec idn -> index 
        self.dec_map = dec_map 
        # sec idn -> opt. dim -> <LeakInfo>
        self.leak_map = leak_map
        # sec idn -> expected Pr. value 
        self.expected_pr = None
        return

    def load_expected_Pr(self,epr):
        self.expected_pr = epr

    # TODO: test this. 
    """
    i2hm := dict, sec idn -> opt. dim. -> sequence of <HypStruct>
    """
    def apply_leakmap_on_IRC2HypStruct_map(self,i2hm):

        for (k,v) in self.leak_map.items():
            if k not in i2hm: 
                continue
            v2 = i2hm[k]

            for k_,v_ in v.items():
                if k_ in v2:
                    # search for the opt of sec 
                    hs = v2[k_]

                        ### single update 
                    """
                    hsx = None
                    d = self.dec_map[hs]
                    for hs_ in hs: 
                        if hs_.seq_idn == d:
                    """
                        ### all update
                    hsx = []
                    for hs_ in hs:
                        hs2 = HypInfer.infer_by_LeakInfo(hs_,v_)
                        hsx.append(hs2) 
                    v2[k_] = hsx 
                        
            i2hm[k] = v2
        return i2hm 

    ############################### generators: <IsoRingedChain>-to-<HypStruct> map

    # TODO: test this.
    """
    - full_hypseq := bool, True outputs a map that includes 
                         all local optima. False outputs a map
                         with only the correct optimum index. 
    - naive_split := int, positive, partitioning scheme for the super-bounds; 
        (<HypStruct> --> the super-bound of a <IsoRing>'s <Sec>'s `seq` value). 
    - hop_size := int, positive, the number of hops per dimension for an RSSI.
    
    return:
    - sec idn -> sec dim. -> <HypStruct>
    """
    @staticmethod
    def naive_IRC2HypStruct_map(irc,full_hypseq=True,\
        naive_split=2,hop_size=5):
        assert type(irc) == IsoRingedChain

        l = defaultdict(None)
        for irl_ in irc.irl:
            m = BackgroundInfo.naive_IsoRing2HypStruct_map(irl_,\
                full_hypseq,naive_split,hop_size)
            l[irl_.sec.idn_tag] = m 
        return l 

    '''
    return: 
    - dict, sec. dimension -> <HypStruct>
    '''
    @staticmethod 
    def naive_IsoRing2HypStruct_map(ir:IsoRing,full_hypseq,\
        naive_split=1,hop_size=5):

        d = defaultdict(list)  
        for s in ir.sec_cache:
            lx = s.dim()
            z = np.zeros((lx,2),dtype=float)
            z[:,1] = 1.0
            sb = even_bound_split(z,naive_split)
            sb_pr = np.ones((naive_split,)) * 1.0/naive_split
            hs_vec = np.ones((naive_split,)) * hop_size
            if not full_hypseq: 
                si = s.seq_index()
                hs = HypStruct(ir.sec.idn_tag,\
                    si,deepcopy(sb),deepcopy(sb_pr),\
                    deepcopy(hs_vec))
                d[lx] = [hs] 
                continue 

            q = [] 
            for i in range(len(s.opm)):
                hs = HypStruct(ir.sec.idn_tag,\
                    i,deepcopy(sb),deepcopy(sb_pr),deepcopy(hs_vec))
                q.append(hs)
            d[lx] = q 
        return d 

    """

    - arguments: 
    irc := IsoRingedChain
    bound_length := arg. for functions to generate <HypStruct>s
    rd_range := [f0,f1], ordered range for the ratio of dummy <HypStruct>s; f_i::float
    ra_range := [f0,f1], ratio of approximate <HypStruct>s; f_i::float

    - return: 
    Sec idn. -> Sec dim. -> HypStruct
    """
    @staticmethod
    def partially_naive_IRC2HypStruct_map(irc,\
        bound_length,rd_range,ra_range,rnd_struct):
        
        assert type(irc) == IsoRingedChain
        assert min(rd_range) >= 0.0 and max(rd_range) <= 1.0
        assert rd_range[1] >= rd_range[0]
        assert min(ra_range) >= 0.0 and max(ra_range) <= 1.0
        assert ra_range[1] >= ra_range[0]

        bx1,bx2 = {},{}
        for ir in irc.irl:
            r_d = rnd_struct.uniform(rd_range[0],rd_range[1])
            max_diff = 1.0 - r_d
            rg1 = min([ra_range[1],max_diff])

            if rg1 < ra_range[0]:
                mx = 0.0 
            else:
                mx = ra_range[0]                
            r_a = rnd_struct.uniform(mx,rg1) 

            b1,b2 = BackgroundInfo.partially_naive_IsoRing2HypStruct_map(\
                ir,bound_length,r_d,r_a,rnd_struct)

            bx1[ir.sec.idn_tag] = deepcopy(b1)
            bx2[ir.sec.idn_tag] = deepcopy(b2)

        return bx1,bx2 

    """
    `full_hypseq` var from the previous method is set to TRUE in this 
    method.

    - arguments: 
    r_d := ratio of dummy <HypStruct>s
    r_a := ratio of approximate <HypStruct>s

    - return: 
    M1 := map, sec dim -> <HypStruct>
    M2 := map, sec dim -> (HypStruct gen-scheme,*args for scheme)
    """
    @staticmethod 
    def partially_naive_IsoRing2HypStruct_map(ir:IsoRing,\
        bound_length,r_d,r_a,rnd_struct):
        assert min([r_d,r_a]) >= 0.0
        assert r_d + r_a <= 1.0 

        sds__ = ir.secdim_seq()

        def op_fx(fxid,seci):
            fgx = DEFAULT_HYPSTRUCT_GENSCHEMES[fxid]
            if fxid in {0,2}:
                ir.set_isorep(seci)
                return fgx(ir),None

            d = sds__[seci]
            lr = [rnd_struct.random() for d__ in range(d)]
            lr = np.array(lr)
            hs = rnd_struct.randrange(2,13)
            hypStruct = fgx(ir,seci,bound_length,lr,hs)
            return hypStruct,hs 

        # choose the dim 
        qx = [[],[],[]]
        i = 0

        # start w/ dummy
        l0 = int(round(r_d * len(sds__))) 
        l1 = int(round(r_a * len(sds__))) 
        l2 = len(sds__) - (l0 + l1) 
        # dummy, approximate 
        lx = [l0,l1,l2]

        sds = [vxx for vxx in enumerate(sds__)]
        outp1,outp2 = {},{}
        while i < 3:
            for _ in range(lx[i]):
                # choose an index
                ix = rnd_struct.randrange(0,len(sds))
                elmnt = sds.pop(ix)
                hs,hs2 = op_fx(i,elmnt[0]) 
                outp1[elmnt[1]] = [hs]
                outp2[elmnt[1]] = [i,hs2]
            i += 1 
        return outp1,outp2

    # TODO: test 
    @staticmethod
    def naive_populate_IRC2HypStruct_map(irc,dx,populate_ratio_range,rnd_struct=random,\
        scramble:bool=False):
        for ir in irc.irl:
            dx_ = dx[ir.sec.idn_tag]
            f = rnd_struct.uniform(populate_ratio_range[0],populate_ratio_range[1])
            rs = rnd_struct if scramble else None
            dx_ = BackgroundInfo.naive_populate_IsoRing2HypStruct_map(ir,dx_,f,rs)
            dx[ir.sec.idn_tag] = dx_ 
        return dx 
    
    """
    for a map D: sec idn -> sec dim -> [<HypStruct>], fills in 
    each dimension's sequence of <HypStruct> with dummy
    <HypStruct> instances for unaccounted optima indices.
    """
    @staticmethod
    def naive_populate_IsoRing2HypStruct_map(ir,dx,populate_ratio:float,\
        rnd_struct=None):

        assert populate_ratio >= 0.0 and populate_ratio <= 1.0 
        for (i,x) in enumerate(ir.sec_cache):
            d = x.dim()
            assert d in dx, "map has incorrect form."
            
            # gather the remaining local optima to be
            # accounted by <HypStruct>
            q = len(x.opm)
            rem = set([i for i in range(q)])
            
            q2 = dx[d]
            q2 = set([h.target_index for h in q2])
            rem = rem - q2
            ##print("Q2: ",q2)
            ##print("REM: ", rem)

                # apply `populate_ratio` to the number of possible 
                # local optima to be accounted for. 
            remsz = int(round(len(rem) * populate_ratio)) 
            ##print("REMSZ: ",remsz)

            rem = sorted(list(rem))[:remsz] 

            # set the isorep, and generate dummy <HypStruct>s for each
            ir.set_isorep(i)
            for r in rem:
                hs = one_dummy_HypStruct_for_IsoRing(ir,r)
                dx[d].append(hs)

            if type(rnd_struct) != type(None):
                rnd_struct.shuffle(dx[d]) 
        
        return dx  

    ########################### generator using <IsoRingedChain>+<SRefMap>

    # NOTE: default generator. 
    """
    generates a <BackgroundInfo> instance from an <IsoRingedChain> 
    and a corresponding <SRefMap>.

    Details of variable values derived from the arguments: 
    - dependency map (accurate), <Sec> idn. --> [other <Sec> idns required to be processed.]
    - co-dependency map (accurate), <Sec> idn. --> [co-dependent <Sec>]
    - decision map (?accurate?), uses <SRefMap> to produce a decision-map by
        arguments ('cd',max,[0,1])
    - leak map (of leaked info.); uses the function <BackgroundInfo.default_IsoRing_leak>
        on each <IsoRing> for the map M: sec idn. -> opt. dim. -> leaked info.
        The argument `target_ratio` determines the number of opt. dim. per <Sec> idn 
        that are assigned leaks. 
    """
    @staticmethod
    def generate_instance(irc,srm,add_noise_to_Pr:bool=False,\
        ooc_malpermute_degree:float=0.0,target_ratio:float=None,\
        vary_lk_sz:bool=False,srm_decmap_args=('cd',max,[0,1]),\
        rnd_struct=random):

        assert type(irc) == IsoRingedChain
        assert type(srm) == SRefMap

        # sec idn -> sec dim. -> opt. idn -> Pr. 
        dm_pr = defaultdict(None)
        for i in irc.irl:
            it = i.sec.idn_tag
            dm_pr[it] = i.dim_to_opt_pr_map()
        ##dm_pr = srm.collect_prism_points__PrMap('cd',"greedy-dc",1)

        cs = connected_subsets_of_codepmap(srm.cdms)
        keys = srm.dms.keys() 
        dx = {} 
        for k in keys:
            ##print("K: ",k)
            dx[k] = depchain_for_Sec(srm.dms,k)
            ##print(dx[k])

        # fetch the decision map
        dec_map = srm.collect_prism_points__DecMap(*srm_decmap_args)#'cd',max,[0,1])
        stat = set(dec_map.keys()) == set(srm.dms.keys())

        # fill in the map for missing decisions; default is
        # max Pr.
        if not stat:
            qx = srm.directF_proc__best_nodedec_map()
            for k,v in qx.items():
                if k not in dec_map:
                    dec_map[k] = v
        lm = BackgroundInfo.generate_background_leak(irc,random,\
            targets=target_ratio,vary_lk_sz=vary_lk_sz)

        # perform alterations onto <BackgroundInfo> vars here
        if add_noise_to_Pr:
            dm_pr = noise_on_opt_pr_map_SEEDTYPE_PythonANDNumpy(dm_pr,rnd_struct)

        if ooc_malpermute_degree > 0.0:
            dx,cs = BackgroundInfo.permute_depmap_AND_codepsets(dx,cs,\
                ooc_malpermute_degree,rnd_struct)

        return BackgroundInfo(dm_pr,dx,cs,dec_map,lm)

    """
    used for malpermutating a calculated (precise) Order-of-Cracking for
    a <BackgroundInfo> generated from a (<IsoRingedChain>,<SRefMap>). 

    dm := dict, sec idn -> sec dim. -> opt. idn -> Pr. 
    cds := list(set(codependent nodes))
    occ_malpermute_degree := float, 0. <= x <= 1.0
    rnd_struct := ? 
    """
    @staticmethod
    def permute_depmap_AND_codepsets(dm,cds,occ_malpermute_degree,rnd_struct):

        def delta_one(vx):
            ld = np.array([len(v_) for v_ in vx])
            num_total_swaps = np.cumprod(ld)
            if len(num_total_swaps) == 0:
                return None
            num_total_swaps = num_total_swaps[-1]
            num_total_swaps += np.sum(ld) 

            swaps = int(round(num_total_swaps * occ_malpermute_degree))
            vx = permute_setseq(vx,rnd_struct,swaps)[0]
            return vx 

        for k in dm.keys():
            v = dm[k]
            q = delta_one(v)
            if type(q) != type(None):
                dm[k] = deepcopy(q)

        q = delta_one(cds)

        if type(q) != type(None): cds = deepcopy(q)
        return dm,cds 
    
    """
    return:
    - dict,
        [0] sec idn. tag 
        [1] [0] dim. of opt.
            [1]  <LeakInfo> 

    """
    @staticmethod
    def generate_background_leak(irc,rnd_struct,\
        targets=None,vary_lk_sz:bool=False):
        assert type(irc) == IsoRingedChain
        leak_map = {}
        ##print("LENGO: ",len(irc.irl))
        for irl_ in irc.irl:
            ##print("IRL_: ",irl_.sec.idn_tag)
            dx = BackgroundInfo.leak_process_IsoRing(irl_,rnd_struct,\
                targets,vary_lk_sz) 
            leak_map[irl_.sec.idn_tag] = dx
        return leak_map

    """
    return: 
    - dict, 
            [0] dimension of opt.
            [1] <LeakInfo> 
    """
    @staticmethod 
    def leak_process_IsoRing(ir,rnd_struct,targets=None,vary_lk_sz:bool=False):
        if type(targets) == type(None):
            targets = ir.secdim_seq()
        elif type(targets) == float:
            assert targets >= 0.0 and targets <= 1.0
            t = ir.secdim_seq() 
            l = int(round(len(t) * targets))
            targets = []
            while len(targets) < l:
                ix = rnd_struct.randrange(0,len(t))
                targets.append(t.pop(ix)) 
        l = 2 if not vary_lk_sz else \
            rnd_struct.randrange(DEFAULT_LEAKSZ_RANGE[0],\
            DEFAULT_LEAKSZ_RANGE[1] + 1)
                
        lk = BackgroundInfo.default_IsoRing_leak(rnd_struct,l=l)
        d = {}
        ##print("IR SEC: ",ir.sec.idn_tag)
        for i in range(len(ir.sec_cache)):
            ir.set_isorep(i)
            q = len(ir.rep().seq)
            # not a target dimension 
            if q not in targets: continue 

            lk = BackgroundInfo.default_IsoRing_leak(rnd_struct)
            lk.leak_info(ir)
            ir.leak_stage -= 1

            if ir.sec.idn_tag not in lk.leakm.d: 
                continue 
            leakInfo = lk.leakm.d[ir.sec.idn_tag]

            d[q] = leakInfo
        return d

    @staticmethod
    def default_IsoRing_leak(rnd_struct,l=2):
        assert l >= DEFAULT_LEAKSZ_RANGE[0] and l <= DEFAULT_LEAKSZ_RANGE[1]
        return Leak.generate_Leak__type1(l,rnd_struct)

class CrackSoln:

    def __init__(self):
        self.d = defaultdict(None)
        self.prev_soln = defaultdict(list) 
        return

    """
    - t := (sec idn.,opt. idn,expected sec value,pr of score)
    """
    def __add__(self,t:tuple):
        assert type(t) == tuple
        assert len(t) == 4
        qx = deepcopy(self) 
        # case: update soln
        if t[0] in qx.d:
            qx.prev_soln[t[0]].append(qx.d[t[0]])
        qx.d[t[0]] = (t[1],t[2],t[3])
        self = qx 
        return qx

    def __len__(self):
        return len(self.d)

    def match_pr(self):
        return -1

class Cracker:

    def __init__(self,hyp_map,backgroundInfo,\
        crackling_sz:int,radar_radius:int=4,\
        energy=1000.0):
        assert type(hyp_map) in {dict,defaultdict,type(None)}
        assert type(backgroundInfo) == BackgroundInfo
        assert type(crackling_sz) == int and crackling_sz > 0
        
        # <sec idn> -> <sec dim> -> list::<HypStruct> 
        self.hyp_map = hyp_map 
        # holds <HypStruct> instances attempted from `hyp_map`
        self.hyp_map_cache = defaultdict(None)
        self.bi = backgroundInfo
        self.crackling_sz = crackling_sz
        self.radar_radius = radar_radius
        self.energy = NerG(energy)
        self.cracklings = [] 
        self.spent = [] 
        self.cidn_counter = 0

        self.calculate_oop()
        self.csoln = CrackSoln() 

        self.initiated = False 

        self.terminated = self.is_terminated()
        return

    # TODO: test 
    # NOTE: for real-life applications, 
    ## the calculations in this method are to be used
    ## sparingly, since the <generator method> wraps 
    ## around objects that normally should never be 
    ## possessed by a 3rd-party intruder's computer 
    ## program.  
    """
    generates an instance based on arguments to two 
    important components, <BackgroundInfo> and 
    <IRC2HypStruct> maps. 

    bi_args := [0] add_noise_to_Pr:bool,
               [1] ooc_malpermute_degree:float,
               [2] target_ratio:float=None,
               [3] vary_lk_sz:bool,
               [4] srm_decmap_args,

    irc2hs_args := 
               [0] naive|partially_naive
                
                case: naive 
               [1] full_hypseq:bool,
               [2] naive_split:int, > 0 
               [3] hop_size:int, > 0

                case: partially naive 
               [1] bound_length,
               [2] rd_range,
               [3] ra_range 

    return: 
    - <Cracker> 
    """
    @staticmethod 
    def generate_instance_by_engineered_BI(irc,srm,rnd_struct,bi_args,irc2hs_args,\
        crackling_sz:int,radar_radius:int,energy:float):

        assert irc2hs_args[0] in {"naive","partially naive"} 

        bi = BackgroundInfo.generate_instance(irc,srm,\
            bi_args[0],bi_args[1],bi_args[2],bi_args[3],\
            bi_args[4],rnd_struct) 

        irc2hs = None
        if irc2hs_args[0] == "naive":
            irc2hs = BackgroundInfo.naive_IRC2HypStruct_map(\
                irc,irc2hs_args[1],irc2hs_args[2],irc2hs_args[3])
        else: 
            irc2hs = BackgroundInfo.partially_naive_IRC2HypStruct_map(\
                irc,irc2hs_args[1],irc2hs_args[2],irc2hs_args[3],rnd_struct)[0]
        cr = Cracker(irc2hs,bi,crackling_sz,radar_radius,energy) 
        return cr 

#####################################################################

    # NOTE: pickle does not save crackling instances
    def pickle_thyself(self,fp_base):
        if not os.path.isdir(fp_base):
            os.mkdir(fp_base)
        fp1 = fp_base + "/hm"
        fp2 = fp_base + "/bi"
        fp3 = fp_base + "/extra"

        fpx1 = open(fp1,"wb")
        pickle.dump([self.hyp_map,self.hyp_map_cache],fpx1)
        fpx1.close() 

        fpx2 = open(fp2,"wb")
        pickle.dump(self.bi,fpx2)
        fpx2.close()

        fpx3 = open(fp3,"wb")
        pickle.dump([self.crackling_sz,self.radar_radius,\
            self.cidn_counter,self.csoln],fpx3)
        fpx3.close()

    @staticmethod
    def unpickle_thyself(fp_base):
        fpx1 = open(fp_base + "/hm","rb")
        hms = pickle.load(fpx1)
        fpx1.close()

        fpx2 = open(fp_base + "/bi","rb")
        bi = pickle.load(fpx2)
        fpx2.close()

        fpx3 = open(fp_base + "/extra","rb")
        extra = pickle.load(fpx3)
        fpx3.close()

        c = Cracker(hms[0],bi,extra[0],extra[1])

        c.cidn_counter = extra[2]
        c.csoln = extra[3]
        c.hyp_map_cache = hms[1] 
        return c


    '''
    calculate the order of operations based on 
    <BackgroundInfo>.

    return:
    - dict, sec idn. -> integer index in ordered sequence OR 
                        np.nan for any location. 
    '''
    def calculate_oop(self): 
        ooc = OrderOfCrackn()
        soln = ooc.order_by_depchain_map(self.bi.dm)
        self.oop = soln 
        self.oopi = 0
        return

    def fetch_crackling(self,cidn):
        i = self.fetch_crackling_index(cidn)
        if i == -1: return None
        return self.cracklings[i]

    def fetch_crackling_index(self,cidn):
        for (i,c) in enumerate(self.cracklings):
            if c.cidn == cidn:
                return i
        return -1

    ######################## methods for loading 
    ######################## <Crackling> instances

    """
    return:
    - set, idn of <Sec> to target
    """
    def next_target(self):

        if self.oopi >= len(self.oop):
            return None
        return self.oop[self.oopi]

    """
    loads the appropriate number of <Crackling> 
    instances for the target set of <Sec> instances, 
    `targetdim_seq`.

    targetdim_seq := list<(sec_idn,sec_dim)>
    """
    def load_cracklings_for_secset(self,targetdim_seq):
        self.cracklings.clear() 
        for (nt_,d) in targetdim_seq: 
            self.load_crackling(nt_,d)
        return

    # TODO: uses <NerG> 
    """
    load 
    """
    def load_cracklings_for_secset_(self,targetdim_seq):

        return -1

    def load_crackling(self,sec_idn,sec_dim,sz_limit=20):

        hs = self.next_hypstruct(sec_idn,sec_dim)
        if type(hs) == type(None):
            return False

        cr = Crackling(cidn=self.cidn_counter,cvsz=sz_limit)
        self.cidn_counter += 1
        cr.load_HypStruct(hs) 
        self.cracklings.append(cr) 
        return cr

    """
    return:
    - <HypStruct>, the next attempt given the variables
                   `sec_idn`,`sec_dim`. 
    """
    def next_hypstruct(self,sec_idn,sec_dim):
        if sec_idn not in self.hyp_map:
            return None
        print("FAIL")
        if sec_dim not in self.hyp_map[sec_idn]:
            return None

        if len(self.hyp_map[sec_idn][sec_dim]) == 0:
            return None

        hs = self.hyp_map[sec_idn][sec_dim][0].new_HypStruct_by_next_subbound()
        if type(hs) == type(None):
            self.hyp_map[sec_idn][sec_dim].pop(0)
        return hs

    #### TODO: new section; needs to be tested. 
    ########################## C2C communication methods

    ##### new section on C2C comm. 


    def calculate_probable_node_destinations():
        return -1 

    def gauge_velocity():

        return -1

    #### TODO: make a running fit

    def target_secidn_set(self):
        s = set()

        for c in self.cracklings:
            td = c.td
            if type(td) == type(None): continue
            tdx = td.td
            if type(tdx) == type(None): continue
            s = s | {tdx.target_node}
        return s

    def cracklings_by_targetsec(self,i):
        cx = []

        for c in self.cracklings:
            td = c.td
            if type(td) == type(None): continue
            tdx = td.td
            if type(tdx) == type(None): continue

            if tdx.target_node == i:
                cx.append(c)
        return cx

    def crackling_sg_supernodeset(self,targetsec_idn):
        cx = self.cracklings_by_targetsec(targetsec_idn)
        ns = set()
        for c in cx:
            if type(c.td.resource_sg) == type(None):
                continue
            ns = ns | set(c.td.resource_sg.d.keys())
        return ns 

    # TODO: test. 
    """
    Uses a supergraph `G` containing the cumulative nodes
    and edges of <Crackling>s targetting the <IsoRing> with
    identifier `targetsec_idn`.

    return:
    - set, sec idn. that have had their `node_path` updated.
    """
    def supergraph_info_update_to_Cracklings(self,G,targetsec_idn):
        assert type(G) == SNGraphContainer

        cx = self.cracklings_by_targetsec(targetsec_idn)

        # iterate through each <Crackling>; for each one
        # that does not know the location of the target, 
        # update info.
        if targetsec_idn not in G.ring_locs:
            return set() 

        loc = G.ring_locs[targetsec_idn]

        updated = set() 
        for c in cx:
            # check if 
            sgx = c.td.resource_sg 
            assert type(sgx) != type(None)
            tdx = c.td.td
            q = td.search_for_target_node(sgx)

            if type(q) == type(None):
                lx = c.loc()
                dfsc = G.sp[lx]

                if loc not in dfsc.min_paths:
                    continue

                mp = dfsc.min_paths[loc][0] 
                n2 = Cracker.intermediate_node(c,mp)
                mp2 = deepcopy(dfsc.min_paths[n2][0])
                c.assign_coordinated_path(mp2)

                updated = updated | {c.cidn}
        return updated

    """
    calculates the intermediate node given a <NodePath>
    starting from <Crackling> location, in which the 
    <NodePath> may have an endpoint not included by 
    <Crackling>'s <TDirector> variable <resource_sg>.
    """
    @staticmethod
    def intermediate_node(crackling,nodepath):
        q = nodepath.p
        assert q[0] == crackling.loc()
        assert type(crackling.td.resource_sg) != type(None)
        
        sx = crackling.td.resource_sg 
        m = None
        if q[-1] in sx.d:
            m = q[-1]
        else:
            for q_ in q:
                if q_ in sx.d:
                    m_ = q_
                    continue
                break
        return m 

    ########################## performance feedback methods


    ########################## network entry methods
    """
    return:
    - bool; response to accept <TDirector> given the 
            <HypStruct> of <Crackling> `cidn`. 
    """
    def accept_TDirector_at_entry_point(self,cidn,td):
        assert type(td) == TDirector
        assert td.vp() == "C"
        assert td.obj_stat == "search for target"

        c = self.fetch_crackling(cidn)
        assert type(c) != type(None)
        return td.check_obj()

    def load_TDirector(self,cidn,td):
        c = self.fetch_crackling(cidn)
        assert type(c) != type(None)
        c.load_TDirector(td)

    def cstat(self):
        d = {}
        for i in range(len(self.cracklings)):
            q = self.cracklings[i].cidn
            d[q] = self.crackling_stat(i)
        return d

    def crackling_stat(self,index):
        # cracked
        if self.cracklings[index].astat:
            return 2
        
        # cracking
        if self.cracklings[index].cstat:
            return 1

        # interdiction
        if self.cracklings[index].istat:
            return 0

        # continue travelling
        return -1 

    def remove_cracklings__fstat(self):
        for c in self.cracklings:
            if c.fstat == True:
                self.remove_spent_crackling(c.cidn)
        return -1

    def remove_spent_crackling(self,cidn):
        i = self.fetch_crackling_index(cidn)
        if i == -1: return
        q = self.cracklings.pop(i)
        self.spent.append(q) 
        return

    def clear_open_info(self):
        for c in self.cracklings:
            c.td.td.open_info_var.clear() 

    def is_terminated(self):
        self.terminated = self.energy.v <= 0.0 
        return self.terminated 
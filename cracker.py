from crackling import * 
from secnet import * 
from leakette import * 

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

class OrderOfCrackn:

    def __init__(self):
        self.seq = []
        # index of set element -> & (for co-dep.) OR | (for non-duplicate choice in set)
        self.set_operator_map = {}
        return

    # TODO: finish 
    def order_by_depchain_map(self,dcm):
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
        # sec idn -> sec dim. -> opt. idn -> Pr. 
        self.opm = opm
        self.dm = depchain_map
        self.cdm = codep_sets
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
    i2hm := dict, sec idn -> 
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
                        ### all udpate
                    hsx = []
                    for hs_ in hs:
                        hs2 = HypInfer.infer_by_LeakInfo(hs_,v_)
                        hsx.append(hs2) 
                    v2[k_] = hsx 
                        
            i2hm[k] = v2
        return i2hm 

    # TODO: test this.
    """

    full_hypseq := bool, True outputs a map that includes 
                         all local optima.

    return:
    - sec idn -> sec dim. -> <HypStruct>
    """
    @staticmethod
    def naive_IRC2HypStruct_map(irc,full_hypseq=True,\
        naive_split=2):
        assert type(irc) == IsoRingedChain

        l = defaultdict(None)
        for irl_ in irc.irl:
            m = BackgroundInfo.naive_IsoRing2HypStruct_map(irl_,\
                full_hypseq,naive_split)
            l[irl_.sec.idn_tag] = m 
        return l 

    # TODO: test this.
    '''
    outputs a map 
        sec. dimension -> <HypStruct>
    '''
    @staticmethod 
    def naive_IsoRing2HypStruct_map(ir:IsoRing,full_hypseq,\
        naive_split=1):

        d = defaultdict(list)  
        for s in ir.sec_cache:
            lx = s.dim()
            z = np.zeros((lx,2),dtype=float)
            z[:,1] = 1.0
            sb = even_bound_split(z,naive_split)
            sb_pr = np.ones((naive_split,)) * 1.0/naive_split
            if not full_hypseq: 
                si = s.seq_index()
                hs = HypStruct(ir.sec.idn_tag,\
                    si,deepcopy(sb),deepcopy(sb_pr))
                d[lx] = [hs] 
                continue 

            q = [] 
            for i in range(len(s.opm)):
                hs = HypStruct(ir.sec.idn_tag,\
                    i,deepcopy(sb),deepcopy(sb_pr))
                q.append(hs)
            d[lx] = q 
        return d 

    # TODO: test this.
    @staticmethod
    def generate_instance(irc,srm):
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

        dec_map = srm.collect_prism_points__DecMap('cd',max,[0,1])
    
        lm = BackgroundInfo.generate_background_leak(irc,random)
        return BackgroundInfo(dm_pr,dx,cs,dec_map,lm)

    ###

    """
    return:
    - dict,
        [0] sec idn. tag 
        [1] [0] dim. of opt.
            [1]  <LeakInfo> 

    """
    @staticmethod
    def generate_background_leak(irc,rnd_struct):
        assert type(irc) == IsoRingedChain
        leak_map = {}
        ##print("LENGO: ",len(irc.irl))
        for irl_ in irc.irl:
            ##print("IRL_: ",irl_.sec.idn_tag)
            dx = BackgroundInfo.leak_process_IsoRing(irl_,rnd_struct)
            leak_map[irl_.sec.idn_tag] = dx
        return leak_map

    """
    return: 
    - dict, 
            [0] dimension of opt.
            [1] <LeakInfo> 
    """
    @staticmethod 
    def leak_process_IsoRing(ir,rnd_struct,targets=None):
        if type(targets) == type(None):
            targets = ir.secdim_seq() 

        lk = BackgroundInfo.default_IsoRing_leak(rnd_struct)
        d = {}
        ##print("IR SEC: ",ir.sec.idn_tag)
        for i in range(len(ir.sec_cache)):
            ir.set_isorep(i)
            lk = BackgroundInfo.default_IsoRing_leak(rnd_struct)
            lk.leak_info(ir)
            ir.leak_stage -= 1

            #print("LK-LEAKM") 
            #print(lk.leakm.d)
            if ir.sec.idn_tag not in lk.leakm.d: 
                continue 
            leakInfo = lk.leakm.d[ir.sec.idn_tag]
            q = len(ir.rep().seq)
            # not a target dimension 
            if q not in targets: continue 

            d[q] = leakInfo
        return d

    @staticmethod
    def default_IsoRing_leak(rnd_struct):
        #assert type(ir) == IsoRing
        return Leak.generate_Leak__type1(2,rnd_struct)

class CrackSoln:

    def __init__(self):
        self.d = defaultdict(None) 
        return

    """
    - t := (sec idn.,expected sec value,pr of score)
    """
    def __add__(self,t:tuple):
        assert type(t) == tuple
        assert len(t) == 3
        return -1 

    def match_pr(self):
        return -1 

CRACKLING_STAT = {"stationary","mobile"}

class Cracker:

    def __init__(self,hyp_map,backgroundInfo,\
        crackling_sz:int):
        assert type(hyp_map) in {dict,type(None)}
        assert type(backgroundInfo) == BackgroundInfo
        assert type(crackling_sz) == int and crackling_sz > 0
        self.hyp_map = hyp_map 
        self.bi = backgroundInfo
        self.cracklings = [] 

        self.calculate_oop()
        return

    '''
    calculate the order of operations based on 
    <BackgroundInfo>.

    return:
    - dict, sec idn. -> integer index in ordered sequence OR 
                        np.nan for any location. 
    '''
    # TODO: complete
    def calculate_oop(self): 
        return -1 

    # TODO: complete
    def load_crackling(self):

        # declare the HypStruct using <BackgroundInfo> 

        return -1

    # TODO: complete
    def load_crackling_TDirector(self):
        return -1 

    # TODO: complete
    def crackling_stat(self,c_idn):
        return -1 

from crackling import * 
from secnet import * 
from leakette import * 

def even_bound_split(bs,sz):
    assert matrix_methods.is_proper_bounds_vector(bs)
    assert sz >= 1 and type(sz) == int

    print("BS: ",bs)
    r = (bs[:,1] - bs[:,0]) / sz 
    print("BS0: ",bs[:,0])
    bsx = []
    q = deepcopy(bs[:,0])
    print("Q: ",q)
    for i in range(sz):
        bs_ = q + r
        bsq = deepcopy(np.array([q,bs_]).T)
        bsx.append(bsq)
        q = bs_
    return bsx

class BackgroundInfo:

    """
    """
    def __init__(self,opm,depchain_map,codep_sets,\
        dec_map,leak_map):
        self.opm = opm
        self.dm = depchain_map
        self.cdm = codep_sets
        self.dec_map = dec_map 
        # sec idn -> opt. dim -> <LeakInfo>
        self.leak_map = leak_map
        return

    # TODO: test this. 
    """
    """
    def apply_leakmap_on_IRC2HypStruct_map(self,i2hm):

        for (k,v) in self.leak_map.items():
            if k not in i2hm: 
                continue
            v2 = i2hm[k]

            for k_,v_ in v.items():
                if k_ in v2:
                    hs = v2[k_]
                    hs2 = HypInfer.infer_by_LeakInfo(hs,v_)
                    v2[k_] = hs2
            i2hm[k] = v2
        return i2hm 

    # TODO: test this.
    """
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
            print("SB")
            print(sb)
            print("##--##")
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
        dm_pr = srm.collect_prism_points__PrMap('cd',"greedy-dc",1)

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

class Cracker:

    def __init__(self,hyp_map,backgroundInfo):
        assert type(hyp_map) in {dict,type(None)}
        assert type(backgroundInfo) == BackgroundInfo
        self.hyp_map = hyp_map 
        self.bi = backgroundInfo
        self.cracklings = [] 
        return
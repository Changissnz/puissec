from crackling import * 
from secnet import * 
from leakette import * 

class BackgroundInfo:

    """
    """
    def __init__(self,opm,depchain_map,codep_sets,\
        dec_map,leak_map):
        self.opm = opm
        self.dm = depchain_map
        self.cdm = codep_sets
        self.dec_map = dec_map 
        self.leak_map = leak_map
        return

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
    output 
    """
    @staticmethod
    def generate_background_leak(irc,rnd_struct):
        assert type(irc) == IsoRingedChain
        leak_map = {}
        for irl_ in irc.irl:
            dx = BackgroundInfo.leak_process_IsoRing(irl_,rnd_struct)
            leak_map[irl_.sec.idn_tag] = dx 
        return leak_map

    """
    return: 
    - dict, 
            [0] dimension of opt.
            [1] leak
    """
    @staticmethod 
    def leak_process_IsoRing(ir,rnd_struct):
        lk = BackgroundInfo.default_IsoRing_leak(rnd_struct)
        d = {}
        for i in range(len(ir.sec_cache)):
            lk = BackgroundInfo.default_IsoRing_leak(rnd_struct)
            ir.set_isorep(i)
            lk.leak_info(ir)
            ##print("LK-LEAKM") 
            ##print(lk.leakm.d)
            if ir.sec.idn_tag not in lk.leakm.d: 
                continue 
            leakInfo = lk.leakm.d[ir.sec.idn_tag]
            q = len(ir.rep().seq)
            d[q] = leakInfo
        return d

    @staticmethod
    def default_IsoRing_leak(rnd_struct):
        #assert type(ir) == IsoRing
        return Leak.generate_Leak__type1(2,rnd_struct)

class Cracker:

    def __init__(self,backgroundInfo):
        assert type(backgroundInfo) == BackgroundInfo
        self.bi = backgroundInfo
        self.cracklings = [] 
        return
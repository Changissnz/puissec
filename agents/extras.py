"""
code used to generate examples in paper.
"""
import time 
from .psec_env import * 

def example__node_analysis_1(): 
    # declare the SecNet
    sn = SecNet_sample_CallSia()

    irc = sn.irc
    print("EXPLODING")
    irc.explode(300,2)

    srm = sn.srm

    # declare the Cracker 
    bi = BackgroundInfo.generate_instance(irc,srm)

    i2hm = BackgroundInfo.naive_IRC2HypStruct_map(irc,full_hypseq=False,\
        naive_split=2)
    crck = Cracker(i2hm,bi,6) 

    # declare the SecEnv
    random.seed(3752)
    np.random.seed(3752)

    se = SecEnv(sn,crck,rnd_struct=random,\
            ct_ratio=5000,vb=0)

    se.instantiate_cracker_target()
    se.instantiate_td_for_IRC(5,1.0)

    print("LEN CRACKLINGS: ",len(se.crck.cracklings))
    cr0 = se.crck.cracklings[0]


    irx = se.sn.irc.fetch_IsoRing(2)

    td0 = cr0.td
    td1 = irx.td


    print("TARGETNODE ANALYSIS")

    print("-- CRACKLING: ",td0.loc())
    dx1 = td0.targetnode_analysis(np.min)
    for (k,v) in dx1.items():
        print(k,v)

    print("-- ISORING: ",td1.loc())
    dx2 = td1.targetnode_analysis(np.min)
    for (k,v) in dx2.items():
        print(k,v)

    print("DESTNODE ANALYSIS")

def example__CBridge__naive_hyp(num_attempts=5000,\
    hop_size:int=5, default_obj_index=None):

    random.seed(2004)
    np.random.seed(2004)

    ss = SecSeq_sample_4(num_secs=1,\
            singleton_range=DEFAULT_SINGLETON_RANGE,\
            num_conn=1,min_components=1,max_nconn_ratio = 0.3,\
            drange_max=1)

    sc = ss[0] 
    print(sc[0]) 

    bound = [0.,1.]
    irc = IsoRingedChain(sc,bound,random,71)
    for x in irc.irl:
        x.explode_contents()
    q = irc[0]

    # [10, 6, 3, 2, 8, 9]
    seci = 2
    bound_length = np.array([0.4,0.8,0.5])
    location_ratio = np.array([0.1,0.2,0.6])
    hop_size = 5

    #### NOTE: 
    ## case 1
    hs = one_approximate_HypStruct_for_IsoRing(q,seci,\
        bound_length,location_ratio,hop_size)
    ## case 2
    """
    q.set_isorep(seci)
    hs = one_correct_HypStruct_for_IsoRing(q) 
    """


    # run a <CRridge> cracking session
    hs.hs_vec = np.array([hop_size for x in hs.hs_vec])
    cr = Crackling(cidn=2,cvsz=200)
    cr.load_HypStruct(hs) 

    cb  = CBridge(cr,q,cidn=12,batch_size=100,verbose=True)
    i = 0 
    il = 800 
    while i < il:# and not cr.astat: 
        print("i: ",i)
        next(cb)
        i += 1

    l = len(cr.flagged_pts)

    print("STAT")
    print("- c : ",cr.cstat)
    print("- i : ",cr.istat)
    print("- a: ", cr.astat)

    print("FLAGGED POINTS")
    print(l)

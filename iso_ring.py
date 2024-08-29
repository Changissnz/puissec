from obj_func import * 
from secnet_gen import * 
from cvec import * 
from tdir import * 
from nERg import * 
import os 

class BoundedObjFunc:

    """
    takes a sequence of proper bounds `bounds_seq` w/
    each bound corresponding to a different i'th element
    in `corr_obj`. 

    Main method is `output(reference point, current point)`.
    If `current point` is not in any bounds of `bounds_seq`,
    output function is `default_objf`.
    """
    def __init__(self,bounds_seq,corr_objf,default_objf):

        for c in corr_objf: assert type(c) == ObjFunc
        for b in bounds_seq: assert matrix_methods.\
            is_proper_bounds_vector(b) 
        assert len(bounds_seq) == len(corr_objf)
        assert type(default_objf) == ObjFunc

        self.bounds_seq = bounds_seq
        self.corr_objf = corr_objf
        self.dobjf = default_objf

    def __str__(self): 
        s = ""

        for (i,b) in enumerate(self.bounds_seq):
            s += "bounds\n"
            s += str(b)
            s += "\n"
            s += str(self.corr_objf[i])
            s += "\n"
        s += "--------------"
        return s  

    def pickle_thyself(self,fp):
        cobjf = []
        for x in self.corr_objf: 
            cobjf.append((x.obj_type,x.rnd_seed))
        dobjf = (self.dobjf.obj_type,self.dobjf.rnd_seed)

        rx = [deepcopy(self.bounds_seq),\
            cobjf,dobjf]

        fobj = open(fp,"wb")
        pickle.dump(rx,fobj)
        fobj.close()
        return 

    @staticmethod
    def unpickle_thyself(fp):
        fobj = open(fp,"rb")
        q = pickle.load(fobj)
        fobj.close()

        cx = []
        for q_ in q[1]:
            cx.append(ObjFunc(q_[0],q_[1]))
        dx = ObjFunc(q[2][0],q[2][1])
        return BoundedObjFunc(q[0],cx,dx)
    
    @staticmethod
    def one_simple_BoundedObjFunc(superbound,rnd_seed:int):
        bounds_seq = np.array([deepcopy(superbound)])

        corr_objf = [ObjFunc("e.d.",rnd_seed)]
        default_objf = ObjFunc("e.d.",rnd_seed)

        return BoundedObjFunc(bounds_seq,\
            corr_objf,default_objf)

    @staticmethod
    def generate_BoundedObjFunc(superbound,\
        spacing_ratio_range,outlier_pr_ratio,\
        num_bounds,rnd_seed:int):

        assert type(rnd_seed) == int 
        assert type(num_bounds) == int and num_bounds > 0 

        r = {0:"e.d.",1:"s.m.",2:"r.n."}
        random.seed(rnd_seed)

        bounds_seq = generate_bounds_vector_sequence(\
            superbound,spacing_ratio_range,outlier_pr_ratio,num_bounds)

        corr_objf = []
        for i in range(num_bounds):
            # choose a random obj_type
            obj_type = random.randrange(0,3)
            obj_type = r[obj_type]

            s2 = random.randrange(100)
            objf = ObjFunc(obj_type,s2)
            corr_objf.append(objf)
        
        default_objf = deepcopy(corr_objf[0])

        bof = BoundedObjFunc(bounds_seq,corr_objf,\
            default_objf)
        return bof 

    def output(self,ref,x):
        i = self.bounds_index_of_point(x)
        if i == -1:
            cf = self.dobjf
        else: 
            cf = self.corr_objf[i]
        return cf.output(ref,x)

    """
    return: 
    - index of bounds that `x` falls in OR -1 
    """
    def bounds_index_of_point(self,x):
        for (i,b) in enumerate(self.bounds_seq):
            if matrix_methods.point_in_bounds(b,x):
                return i
        return -1 

class IsoRing:

    def __init__(self,sec:Sec,ofunc:BoundedObjFunc,bounds,\
        cvecl=None):
        assert type(sec) == Sec
        assert type(ofunc) == BoundedObjFunc
        assert matrix_methods.is_proper_bounds_vector(bounds)
        if type(cvecl) != type(None):
            for c in cvecl: assert type(c) == CVec
            assert len(cvecl) == len(sec.opm)

        self.sec = sec
        self.sec_cache = [self.sec] 

        self.ofunc = ofunc 
        self.ofunc_cache = [self.ofunc]

        self.bounds = bounds 
        # bloom stat 
        self.bstat = True
        # index of repr in `sec_cache`
        self.repi = 0

        # cracked stat 
        self.cstat = False 

        # TODO: 
        self.td = None
        self.sec_script = None
        self.leak_stage = 0

    def bounds_of_dim(self,d):
        bds = np.empty((d,2),dtype=float)
        q0 = min(self.bounds[:,0])
        q1 = max(self.bounds[:,1])
        bds[:,0] = q0
        bds[:,1] = q1
        return bds

    def loc(self):
        if type(self.td) == type(None):
            return None
        return self.td.loc() 

    def sec_of_dim(self,d):
        for x in self.sec_cache:
            if x.dim() == d: 
                return x
        return None 

    def dim_to_opt_pr_map(self):
        od = {} 
        for x in self.sec_cache:
            k = x.dim()
            od[k] = x.optima_points_to_index_pr_map()
        return od 

    """
    pickles instance through values
    [F1] sec  
    [F2] list::(sec cache)
    [F3] ofunc 

    [F4] bounds
         bstat 
         repi
         leak_stage
    """
    def pickle_thyself(self,fp_base):

        if not os.path.isdir(fp_base):
            os.mkdir(fp_base)
        fp1 = fp_base + "/__SEC0"
        fp2 = fp_base + "/__SECL"
        fp3 = fp_base + "/__OFUNC"
        fp4 = fp_base + "/__REM"

        # [F1]
        self.sec.pickle_thyself(fp1)

        # [F2]
        qx = []
        for i in range(len(self.sec_cache)):
            q = self.sec_cache[i].to_pickle_list()
            qx.append(q)
        fobj = open(fp2,"wb")
        pickle.dump(qx,fobj)
        fobj.close()

        # [F3]
        self.ofunc.pickle_thyself(fp3)

        # [F4]
        rx = [deepcopy(self.bounds),\
            deepcopy(self.bstat),\
            deepcopy(self.repi),\
            deepcopy(self.leak_stage)]
        fobj = open(fp4,"wb")
        pickle.dump(rx,fobj)
        fobj.close()
        return 

    @staticmethod
    def unpickle_thyself(fp_base):
        fp1 = fp_base + "/__SEC0"
        fp2 = fp_base + "/__SECL"
        fp3 = fp_base + "/__OFUNC"
        fp4 = fp_base + "/__REM"

        # [0]
        s = Sec.unpickle_thyself(fp1)

        # [1] 
        s1 = Sec.unpickle_thyselves(fp2)

        # [2]
        s2 = BoundedObjFunc.unpickle_thyself(fp3) 

        # [3]
        fobj = open(fp4,"rb")
        s3 = pickle.load(fobj)
        fobj.close()


        ir = IsoRing(s,s2,s3[0])
        ir.bstat = s3[1]
        ir.repi = s3[2] 
        ir.leak_stage = s3[3]
        ir.sec_cache = s1 
        return ir 

    def explode_contents(self,optima_size_limit=1000,\
        num_blooms=float('inf'),rnd_struct=random): 

        s = len(self.sec_cache[-1].opm)
        ##print("starting length for {}: {}".format\
        ##    (len(self.sec_cache[-1].seq),s))
        i = 0
        while s < optima_size_limit and i < num_blooms:
            sc = self.sec_cache.pop(-1)
            ##?
            sc.process_one_bloomiso(optima_size_limit - s)
            i += 1 
            s2 = sc.generate_next_Sec()
            s2[0].idn_tag = sc.idn_tag 
            ts2 = s2[2]
            print("transition: {}->{}".format(s2[1],s2[2]))
            s2 = s2[0]
            s = len(s2.opm)
            self.sec_cache.append(sc)
            self.sec_cache.append(s2)

            bof = self.next_bof(s2.dim(),rnd_struct)
            self.ofunc_cache.append(bof)
            if type(ts2) == type(None):
                break
        self.bstat = False
        return

    def next_bof(self,d,rnd_struct=random):
        qx = self.bounds_of_dim(d)

        spacing_ratio_range = rnd_struct.uniform(0.,1.)
        spacing_ratio_range = [spacing_ratio_range,\
            rnd_struct.uniform(spacing_ratio_range,1.)]
        outlier_pr_ratio = rnd_struct.random()
        num_bounds = rnd_struct.randrange(1,6)

        rs = rnd_struct.randrange(200)
        bof = BoundedObjFunc.generate_BoundedObjFunc(qx,\
        spacing_ratio_range,outlier_pr_ratio,\
        num_bounds,rs)
        return bof 


    def secdim_seq(self):
        l = []
        for s in self.sec_cache:
            q = len(s.seq) 
            l.append(q) 
        return l 

    def register_attempt(self,p):
        ##print("-- RA")
        q = self.register_attempt_(p)
        stat = matrix_methods.equal_iterables(p,self.sec.seq)
        outstat = []
        ops = self.rep().optima_points()
        for ps in ops:
            stat2 = matrix_methods.equal_iterables(p,ps)
            outstat.append(stat2)
        ##print("got outstat")
        ##print(outstat)
        outstat = np.any(np.array(outstat) == True) 
        self.cstat = stat
        ##print("outstat")
        ##print(outstat)
        return q,outstat 

    def register_attempt_(self,p):
        assert matrix_methods.is_vector(p) 
        r = self.rep()
        ops = r.optima_points()
        if len(p) != ops.shape[1]:
            print("attempt in wrong dim {}, want {}".format(len(p),\
                ops.shape[1]))
            return None 
        
        # get scores for each of the optima 
        q = list(map(lambda x: self.ofunc.output(x,p),ops))
        return np.array(q)

    def response_to_prompt(self,index):
        if type(self.sec_script) in {defaultdict,dict}:
            return self.sec_script[index]
        m = self.sec.optima_points_to_index_pr_map()
        return m[index]

    def set_isorep(self,i):
        assert len(self.sec_cache) > i and i > -1
        self.repi = i 
        self.sec = self.sec_cache[self.repi]
        self.ofunc = self.ofunc_cache[self.repi]
        return

    def rep(self):
        return self.sec_cache[self.repi]

    #######################################

    def default_TDirector_instance(self,l,radius:int,\
        tdts=DEFAULT_TDIRECTOR_TIMESTAMP_SIZE):
        assert radius > 0
        td = TDirector(l,l,"I",self.sec.idn_tag,\
            radius,1,tdts=tdts)
        self.td = td

    def fetch_td(self):
        if type(self.td) != TDirector:
            return None
        return self.td.td

    #########################################

    # TODO: test
    def default_secproc(self,timespan:float,\
        rnd_struct,verbose=False):
        assert timespan >= 0.0
        assert type(self.td) == TDirector
        
        stat = self.td.check_obj()
        if verbose:
            print("\t\t----------DEFSECPROC,\nI={},OBJ={},STAT={},LOC={}".format(self.sec.idn_tag,\
                self.td.obj_stat,stat,self.td.loc()))
            

        # switch obj.
        if stat:
            s1 = self.td.obj_stat
            self.td.switch_obj_stat()
            s2 = self.td.obj_stat 
            if verbose:
                print("SWITCH {}: {}->{}".format(\
                    self.sec.idn_tag,s1,s2))
                
        # decide how to proceed

        ## case: change from `null radar`; fetch a 
        ##       new and active path.
        if stat and self.td.obj_stat == "avoid target":
            pt = self.td.defsec_path(rnd_struct)

            if verbose:
                print("SET NEW PATH 2SEC @ LOC=",self.td.loc())
                print(pt)

            if type(pt) == type(None):
                return 
            self.td.load_new_path(pt)

        # set active stat of <TDir> to False
        if stat and self.td.obj_stat == "null radar":
            self.td.td.active_stat = False
        
        # travel by <TDir>
        q1 = self.td.loc() 
        self.td.td.scaled__next__(timespan)
        q2 = self.td.loc()

        if verbose: 
            print("MOVELOC {}-->{}\n============".format(q1,q2)) 
        return

################################################

### an example of an <IsoRing> instantiation
def IsoRing_sample_1():
    secs = Sec_list_sample1()
    sndg = SecNetDepGen(secs,random,2,0.75,[1,4])
    sndg.assign_conn()
    sq = sndg.sq

    superbound = np.ones((5,2)) * np.array([0.,1.])
    spacing_ratio_range = [0.,0.2]
    outlier_pr_ratio = 0.4#1.0
    num_bounds = 8
    sb = deepcopy(superbound)

    obf = BoundedObjFunc.generate_BoundedObjFunc(\
            superbound,spacing_ratio_range,\
            outlier_pr_ratio,num_bounds,3) 

    return IsoRing(sq[0],obf,sb)

def SecSeq_sample_1(num_components=1):
    s = Sec_list_sample2()
    sndg = SecNetDepGen(s,random,num_components,0.8,[1,4])
    sndg.assign_conn(2500)
    ss = SecSeq(sndg.sq)
    return ss 

def SecSeq_sample_2(num_secs=80,num_conn=5000,\
    min_components=4,drange_max=4):
    s = Sec_list_sample2(num_secs=num_secs)
    sndg = SecNetDepGen(s,random,min_components,0.5,[1,drange_max])
    print("assigning conn")
    sndg.assign_conn(num_conn)
    
    ss = SecSeq(sndg.sq)
    print("declared <SecSeq>")
    return ss 

def duplicate_Sec_list(sec_list,indices_to_dup,\
    dup_iterations):
    print("duplications remaining: {}".format(dup_iterations))

    if dup_iterations <= 0: return sec_list 

    for i in indices_to_dup:
        l = len(sec_list)
        s = sec_list[i].deepcopy(new_idn_tag=l)
        sec_list.append(s) 
    
    return duplicate_Sec_list(sec_list,\
        indices_to_dup,dup_iterations-1)

def SecSeq_sample_3():
    s = Sec_list_sample2(num_secs=40) 
    i = [i_ for i_ in range(40)] 
    s = duplicate_Sec_list(s,i,10)

    random.seed(22)
    np.random.seed(222)
    sndg = SecNetDepGen(s,random,4,0.5,[1,4],\
        conn_candidate_size=5000) 
    print("assigning conn")
    stat,i = True,0
    while stat and i < 500:  
        print("dconn {}".format(i))
        stat = sndg.make_dep_conn()
        i += 1 
    sndg.write_conn_to_Sec() 
    ss = SecSeq(sndg.sq) 
    return ss,sndg

"""
pre-requisite method for the following:
- SecNet_sample_C3
- SecNet_sample_TDirNv1
"""
def SecSeq_sample_4(num_secs=80,singleton_range=DEFAULT_SINGLETON_RANGE,\
    num_conn=5000,min_components=4,max_nconn_ratio=0.4,drange_max=4):

    s = Sec_list_sample3(num_secs,singleton_range,random)
    print("SECLIST: ",len(s))
    sndg = SecNetDepGen(s,random,min_components,\
            max_nconn_ratio,[1,drange_max])
    sndg.assign_conn(num_conn)

    ss = SecSeq(sndg.sq)
    return ss,sndg 

"""
virtually identical to <SecSeq_sample_4> besides from 
additional argument `depconn_ratio`. 
"""
def SecSeq_sample_5(num_secs=80,singleton_range=DEFAULT_SINGLETON_RANGE,\
    num_conn=5000,min_components=4,max_nconn_ratio=0.4,drange_max=4,\
    depconn_ratio=0.3,conn_types=[1,2,3]):

    s = Sec_list_sample3(num_secs,singleton_range,random)
    sndg = SecNetDepGen(s,random,min_components,\
        max_nconn_ratio,[1,drange_max],depconn_ratio)
    sndg.assign_conn(num_conn,l=conn_types)
    
    ss = SecSeq(sndg.sq)
    return ss,sndg
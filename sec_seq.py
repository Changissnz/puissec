from srefmap import * 

############################################

"""
class helps with transforming optima points 
of a Sec from one dimension to another. 
"""
class OptimaBloomFuncSecRep:

    def __init__(self,obf,sz_map,immutable_dim:int,\
        corr_type="e"):#,bloom_sz_limit:int=1000):
        assert type(obf) == OptimaBloomFunc
        assert type(immutable_dim) == int and immutable_dim >= 0
        #assert type(bloom_sz_limit) == int and bloom_sz_limit > 0

        self.obf = obf
        self.sz_map = sz_map

        self.tdim = [immutable_dim]
        self.corr_type = corr_type

        self.sz = self.obf.oseeds.shape[1]
        self.bpoints = {} 
        self.bpoints[self.sz] = deepcopy(\
            self.obf.oseeds)
        # finished blooming
        self.fstat = False
        # terminate at dim `sz` 
        self.tstat = False
        self.set_sz()

        self.set_dpm() 

    def set_dpm(self):
        self.dpm = DerivatorPrMap() 

    def set_sz(self,attempts=10):
        self.tdim.append(self.sz)
        # get dim of OptimaBloomFunc
        self.sz = self.obf.oseeds.shape[1]
        ##self.tdim.append(self.sz) 

        ##print("SZ[0]: ",self.sz)
        self.sz = self.sz_map(self.sz)
        ##print("SZ[1]: ",self.sz)

        if self.sz in self.tdim:
            ##print("VIOLA")
            if attempts > 0:
                return self.set_sz(attempts - 1)
            self.sz = None 
            return None 

        ##print("SZ IS: {}".format(self.sz))
        return self.sz

    # NOTE: wrong. 
    def iso_appear(self,i):
        assert i in self.bpoints
        return deepcopy(self.bpoints[i])

    def __next__(self):
        self.obf.d = self.sz

        #if self.sz in self.tdim:
        #    return None

        d = next(self.obf)
        if type(d) == type(None):
            ##print("[!!] TSTAT")
            self.tstat = True
            self.tdim.append(self.sz)
            return None

        ##print("D: ", d)
        if self.sz not in self.bpoints:
            ##print("D: ",d)
            self.bpoints[self.sz] = np.array([d])
        else:
            v = self.bpoints[self.sz]
            self.bpoints[self.sz] = np.vstack((v,d))

        self.assign_next_value_to_pred(d)
        ##
        """
        print("BPOINTS")
        print(self.bpoints)
        print("------------")
        """
        ##
        return d

    """
    records derivative indices used to generate 
    `next_val` into `dpm`. 
    """
    def assign_next_value_to_pred(self,next_val):
        px1 = self.obf.prev_i1
        px2 = self.obf.prev_i2
        self.dpm.process_index_pair(px1,px2)
        return

    def finalize_dcount(self,pred_opt2pr_map:defaultdict):
        # get the dimension of predecessor
        ##print("SZ IS: ",self.sz)
        pr_vec,exact_corr = self.dpm.fin_count(self.corr_type,\
            self.sz,pred_opt2pr_map)
        self.set_dpm()
        return pr_vec,exact_corr#,dep_map 

    def sample_size_at_i(self,i):
        if i not in self.bpoints: return 0
        s = self.bpoints[i].shape 
        return s[0] * s[1] 

    def check_fstat(self,dimset):
        self.fstat = set(self.tdim) == dimset 

"""
function to generate an instance of 
<OptimaBloomFuncSecRep>.
"""
def OptimaBloomFuncSecRep__type_1(opt_points,singleton_range,\
        rnd_struct):

    abf = default_AltBaseFunc_for_IsoRing()
    bloom_func = abf.load_AltBaseFunc_function(abf)

    # generate the size-delta function for
    # the <IsoRingBloomer>
    multiplier = rnd_struct.uniform(1.1,20.8)
    szm = generate_default_szmult_mod(multiplier)

    bounds = np.array([deepcopy(singleton_range) for _ in \
        range(opt_points.shape[1])])

    split_sz = rnd_struct.randrange(4,14)
    print("D: ",opt_points.shape)
    obf = OptimaBloomFunc(opt_points,bounds,\
            None,bloom_func,opt_points.shape[1],split_sz) 

    obfsr = OptimaBloomFuncSecRep(obf,szm,\
        opt_points.shape[1])
    return obfsr 

class Sec:

    def __init__(self,sequence,singleton_range,\
        optima_pr_map,dep_map,codep_map,obfsr=None):

        assert matrix_methods.is_vector(sequence)
        assert np.array(singleton_range).ndim == 1 and \
            singleton_range[0] <= singleton_range[1] 

        assert type(optima_pr_map) == type(dep_map) and \
            type(dep_map) == type(codep_map)
        assert type(optima_pr_map) == defaultdict
        assert matrix_methods.vector_to_string(sequence,float) in optima_pr_map

        self.seq = sequence
        self.singleton_range = singleton_range
        self.opm = optima_pr_map
        self.dm = dep_map
        self.cdm = codep_map

        self.obfsr = obfsr
        if type(self.obfsr) == type(None):
            self.declare_obf()
        
        # D2 -> D1; used for storing direct
        # correlations from 
        self.child_exact_opt2opt_map = None

        self.declared_dir = []
        self.next_sz = None 

        self.dimso = None 

        self.idn_tag = None 

    def dim(self): 
        return len(self.seq) 

    def deepcopy(self,new_idn_tag="NULL",\
        transfer_obfsr=False):
        s = deepcopy(self.seq)
        sr = deepcopy(self.singleton_range)
        opm_ = deepcopy(self.opm) 
        dm_ = deepcopy(self.dm)
        cdm_ = deepcopy(self.cdm) 
        obfsr = None 
        if transfer_obfsr: 
            obfsr = self.obfsr
            self.obfsr = None
        sc = Sec(s,sr,opm_,dm_,cdm_,obfsr)
        sc.idn_tag = new_idn_tag
        return sc

    def pickle_thyself(self,fp):
        fobj = open(fp,"wb")
        q = self.to_pickle_list()
        pickle.dump(q,fobj)
        fobj.close()
        return

    def to_pickle_list(self):
        return (self.seq,self.singleton_range,\
            self.opm,self.dm,self.cdm,self.idn_tag)

    @staticmethod
    def unpickle_thyself(f): 
        fobj = open(f,"rb")
        q = pickle.load(fobj)
        fobj.close()
        return Sec.varl_to_Sec(q)
    
    @staticmethod 
    def varl_to_Sec(q): 
        assert len(q) == 6
        s = Sec(q[0],q[1],q[2],q[3],q[4])
        s.idn_tag = q[5]
        return s 

    @staticmethod
    def unpickle_thyselves(fx):
        rx_ = open(fx,"rb")
        rx = pickle.load(rx_)


        s = []
        for i in range(len(rx)):
            q = rx[i]
            sec = Sec.varl_to_Sec(q)
            s.append(sec)
        rx_.close()
        return s

    def __str__(self):
        s = "** sequence\n"
        s += matrix_methods.vector_to_string(self.seq,float)
        s += "\n" + "** optima pr." + "\n"
        s += str(self.opm)
        s += "\n" + "** dep. map" + "\n"
        s += str(self.dm)
        s += "\n" + "** co-dep. map" + "\n"
        s += str(self.cdm)

        return s + "\n"

    """
    bare instances do not have any dep. or co-dep. 

    singleton_range :=
    dimension :=
    num_optima := 
    optima_countermeasure := 
    rnd_struct := random number generator
    """
    @staticmethod
    def generate_bare_instance(singleton_range,dimension,num_optima,\
        optima_countermeasure,rnd_struct=np.random):
        # select the optima points 
        bds = np.array([deepcopy(singleton_range) for _ in range(dimension)])
        ps = generate_pointgen_stdrand(bds,num_optima,rnd_struct)
        ##print(ps)
        
        # select a random optima
        qi = random.randrange(0,len(ps))
        seq = deepcopy(ps[qi])

        # declare a Pr map for the optima
        opm = generate_pr_dist_for_seq(ps,qi,optima_countermeasure,\
            rnd_struct)
        dep_map = defaultdict(float)
        codep_map = deepcopy(dep_map)
        return Sec(seq,singleton_range,\
            opm,dep_map,codep_map,obfsr=None)

    def seq_index(self):
        ops = self.optima_points()

        for (i,o) in enumerate(ops): 
            stat = matrix_methods.equal_iterables(o,self.seq)
            if stat: return i 
        return -1

    """
    return: 
    - np.array,rows ordered by alphanumeric order. 
    """
    def optima_points(self):
        ks = sorted(list(self.opm.keys()))
        optima_points = [matrix_methods.string_to_vector(v,float) for v in \
                    ks]
        optima_points = np.array(optima_points)
        return optima_points

    """
    converts opm map 
        key := stringized point
        value := Pr. value
    to a map w/ 
        key := index of ordering for stringized point
        value := Pr. value 
    """
    def optima_points_to_index_pr_map(self):

        ks = sorted(list(self.opm.keys()))
        ks = [(i,self.opm[k]) for (i,k) in enumerate(ks)]

        d = defaultdict(float,ks) 
        return d

    def declare_obf(self):

        abf = default_AltBaseFunc_for_IsoRing()
        bloom_func = abf.load_AltBaseFunc_function(abf)

        # generate the size-delta function for
        # the <IsoRingBloomer>
        multiplier = 3.13
        szm = generate_default_szmult_mod(multiplier)

        # initialize the <OptimaBloomFunc>
        optima_points = self.optima_points()
        bounds = np.array([deepcopy(self.singleton_range) for _ in \
            range(len(self.seq))])

        obf = OptimaBloomFunc(optima_points,bounds,\
            None,bloom_func,optima_points.shape[0],8)
        self.obfsr = OptimaBloomFuncSecRep(obf,szm,\
            optima_points.shape[1])
        return

    """
    return: 
    - bloom value, [index for element 1, index for element 2]
    """
    def __next__(self):
        q = next(self.obfsr)
        ##print("NEXT: ",q) 
        if type(q) == type(None):
            return None,None

        q2 = np.array([self.obfsr.obf.prev_i1,self.obfsr.obf.prev_i2])
        self.dimso = len(q)
        ##print("DIMSO: ",self.dimso)
        return q, q2

    def process_one_bloomiso(self,sz_limit=float('inf')):
        stat = True 
        c = 0 
        while stat: 
            ##print("one bloom pt.")
            # next value, indices of derivators
            b1,b2 = self.__next__()
            
            ##print("bloom index: ",c)
            
            c += 1
            
            if c >= sz_limit: 
                self.obfsr.tstat = True 
                break 
            
            ##
            """
            print("SEC NEXT")
            print(b1)
            print()
            print(b2)
            print("-----")
            """
            ##
            stat = not (type(b1) == type(None))
            if not stat:
                continue

        
        return

    def lone_pr_vec_for_bloom(self):
        x = self.optima_points_to_index_pr_map()
        return self.obfsr.finalize_dcount(x)
        
        #####

    def generate_next_Sec(self):
        assert not self.obfsr.fstat
        assert self.obfsr.tstat     
        ##print("-- blooming for next sec")    
        q = self.lone_pr_vec_for_bloom()

        # get the optima points for the `sz`
        # dimension
        bps = deepcopy(self.obfsr.bpoints[self.obfsr.sz])
        assert bps.shape[0] == len(q[0])

        optima_pr_map = defaultdict(float)
        sm = sum(q[0])
        sm = 1.0 if sm == 0.0 else sm 
        for i,q_ in enumerate(q[0]):
            vs = matrix_methods.vector_to_string(np.round(bps[i],5),float)
            optima_pr_map[vs] = np.round(q_ / sm,5)
        ##print("-- fetching correlation maps")

        # get the dep.,codep. map
        dep_map = exact_correlation_DMap_key_delta(self.dm,deepcopy(q[1]))
        codep_map = exact_correlation_DMap_key_delta(self.cdm,deepcopy(q[1]))

        # store the pred_exact_corrmap
        self.child_exact_opt2opt_map = deepcopy(q[1])

        ##print("-- fetching new seq.")
        # get the new seq, which has the highest
        # Pr. in optima_pr_map
        qx = max([(k,v) for (k,v) in optima_pr_map.items()],key=lambda x:x[1])
        nu_pt = matrix_methods.string_to_vector(qx[0],float)
        nu_pt = np.round(nu_pt,5)
        ##print("-- size of opt. map: ",len(optima_pr_map))
        s = Sec(nu_pt,deepcopy(self.singleton_range),\
            optima_pr_map,dep_map,codep_map,obfsr=self.obfsr)
        sz0 = s.obfsr.sz

        """
        print("-- RESETTING OSEEDS")
        print("\t-- ORIGINAL")
        print(self.obfsr.obf.oseeds)
        print("\t-- NEW")
        print(bps)
        print()
        """
        ##
        sr_mat = np.array([deepcopy(self.singleton_range) for _ in range(sz0)])
        self.obfsr.obf.reset_oseeds(bps,sr_mat,True)
        self.obfsr = None
        s.obfsr.set_sz()
        s.obfsr.tstat = False 
        return s,sz0,s.obfsr.sz

class SecSeq:

    def __init__(self,sequence):
        for (i,s) in enumerate(sequence): 
            assert type(s) == Sec
            s.idn_tag = i 
        self.sequence = sequence
        return

    def pickle_thyself(self,folder):
        os.mkdir(folder)
        for s in sequence: 
            fx = str(s.idn_tag)
            s.pickle_thyself(folder + "/" + fx) 
        return

    @staticmethod 
    def unpickle_thyself(fp):
        s = [] 
        l = os.listdir(fp)

        for l_ in l: 
            sec = Sec.unpickle_thyself(fp + "/" + l_)
            s.append(sec)
        return SecSeq(s) 
        
    def __getitem__(self,i):
        assert i <= len(self.sequence)
        assert i >= 0 
        return self.sequence[i] 

    def __len__(self):
        return len(self.sequence) 

    def sec_instances_to_supermap(self,map_type):
        assert map_type in {'l','d','c'}

        def selectm(i):
            if map_type == 'l':
                return self.sequence[i].optima_points_to_index_pr_map()
            if map_type == 'd':
                return deepcopy(self.sequence[i].dm)
            return deepcopy(self.sequence[i].cdm)

        dx = defaultdict(defaultdict)
        for i in range(len(self.sequence)):
            dx[i] = selectm(i)
        return dx 

    def soln_map(self):
        d = {}
        for s in self.sequence:
            d[s.idn_tag] = s.seq_index()
        return d 

    # TODO: 
    def construct_SRefMap(self):
        dmsm = self.sec_instances_to_supermap('d')
        return -1 
from defaults import *
from morebs2 import aprng_gauge
from lcg2h import *

"""
The `Euler's number type q` function has 
the form: 

f(1 + c1 * e ** e1) / g(1 + c2 * e ** e2);
e1 and e2 are input variables to output
function `F`. 
"""
def generate_efunc_type_q(c1:float,\
    c2:float):#,f,g):

    base = lambda e1: 1 + c1 * math.e ** e1
    undrbase = lambda e2: 1 + c2 * math.e ** e2
    q = lambda e1,e2: np.round(base(e1) / undrbase(e2),5)
    return q

"""
Kustom Random numerical generator

Has these methods named the same 
as default Python `random` library: 
- uniform
- seed
- randrange
- choice
- shuffle
"""
class KRnd:

    def __init__(self):
        return -1

    def initialize_from(self,args):
        if type(args) == int:
            return -1 
        return -1 

    def random(self): 
        return -1

    def uniform(self):
        return -1

    def choice(self):
        return -1
    
    def shuffle(self):
        return -1 

    def randrange(self):
        return -1

##################################################
########## AltBaseFunc

"""
lightweight database to determine 
terminated base-pair functions belonging
to an <AltBaseFunc> instance. 
"""
class BaseFuncTDetectDB:

    def __init__(self,base_pairs_size:int,tcond1,tcond2,\
        score2bool_func,check_size:int=100):
        self.bps = base_pairs_size
        # function F(boolseq) -> bool::is_terminated
        self.tcond1 = tcond1
        # function G(current output, prev output) -> bool
        self.tcond2 = tcond2
        # function H(coverage score, normed. uwpd score) -> bool
        self.score2bool_func = score2bool_func
        # minimum size of `bp_container[i]` for i'th 
        # base-pair function to be checked 
        self.check_size = check_size 

        # data containers for each base-pair function
        self.bp_container = [[] for _ in range(self.bps)]
        # status log for each base-pair function 
        self.bp_log =  [[] for _ in range(self.bps)]

        # corresponding APRNGGauge's for base-pair function
        self.gauge_structs = []
        for i in range(self.bps):
            ap = aprng_gauge.APRNGGauge(None,\
                DEFAULT_SINGLETON_RANGE,\
                (DEFAULT_SINGLETON_RANGE[1] - DEFAULT_SINGLETON_RANGE[0]) / 25.0)
            self.gauge_structs.append(ap)

        self.prev_vstat = [False for _ in range(self.bps)]
        self.bstat = None

    def bool_stat(self):
        c = self.prev_vstat.count(True)
        self.bstat = c == len(self.prev_vstat)

    """
    outputs a terminating value for every i'th element in `gauge_structs`. 
    """
    def vstat(self):
        self.prev_vstat = []
        for i in range(self.bps): 
            stat = self.stat_at_i(i)
            self.prev_vstat.append(stat) 
        return deepcopy(self.prev_vstat)

    def stat_at_i(self,f_index):
        x = self.bp_log[f_index]

        ## used for quality testing
        """
        if f_index == 2:
            print("bp log @ {}: {}".format(f_index,x))
            print("bp container")
            print(self.bp_container[f_index])
        """
        ##

        return self.tcond1(x) 

    def load_into_cycle_var(self,f_index,info):
        if type(self.prev_vstat) != type(None):
            if self.prev_vstat[f_index] == True:
                return 

        self.bp_container[f_index].extend(info) 
        self.check(f_index)
        self.vstat()
        self.bool_stat() 


    """
    check rule:
    - measures the cycle of `bp_container` at 
      `f_index`; converts the 
      (coverage,normed. uwpd) score to a 
      boolean `b`, logs `b` to `bp_log`.
    """
    def check(self,f_index):
        cycle = np.array(self.bp_container[f_index])
        l = len(cycle)
        if l >= self.check_size:
            # reload the cycle
            self.gauge_structs[f_index].cycle = None
            self.gauge_structs[f_index].assign_cycle(cycle)

            # measure cycle
            mc = self.gauge_structs[f_index].measure_cycle(term_func=self.tcond2)

            # convert measurement to bool
            b = self.score2bool_func(mc) 

            # log bool
            self.bp_log[f_index].append(b)

            # clear the bp_container 
            self.bp_container[f_index].clear() 
            return True
        return False

"""
decision structure to output a float f
given two inputs a,b, by randomly 
selecting a base-pair function B in 
`base_pairs`,
B(a,b) -> f. 

If the class function 
`attach_base_func_tdetect_db`
is called, then AltBaseFunc can
determine which functions of 
`base_pairs` have terminated by a 
pseudo-randomness measure, otherwise 
a default `BaseFuncTDetectDB`.
"""
class AltBaseFunc:

    def __init__(self,base_pairs,rnd_struct,\
        replace_inf=DEFAULT_REPLACE_INF):
        self.base_pairs = base_pairs
        self.rnd_struct = rnd_struct 
        self.replace_inf = replace_inf
        self.bftd = None
        self.bftd_count = 0
        self.bftd_tmap = None
        return

    def attach_base_func_tdetect_db(self,bps:int):
        self.bftd = default_BaseFuncTDetectDB(bps)
        self.bftd_tmap = {}
        for i in range(bps): self.bftd_tmap[i] = -1
        return

    ######################### functions for outputting next value

    def catch_value_error(self,x,y):
        stat1 = type(x) == np.ndarray
        stat2 = type(y) == np.ndarray 
        
        if stat1 or stat2:
            x_ = x if stat1 else y
            if x_.ndim in {1,2}:
                q = np.zeros(x_.shape)
            else: 
                raise ValueError
        else:
            q = 0
        print("CAUGHT VALUE ",q)
        return q

    def replace_inf_(self,output):
        if type(output) != np.ndarray:
            return self.replace_inf(output)
        
        assert output.ndim in {1,2}

        m = lambda x: np.array(list(map(self.replace_inf,x)))
        output_ = None
        if output.ndim == 1:
            output_ = m(output)
        else:
            output_ = np.array([m(ot) for ot in output])
        return np.round(output_,5) 

    def output(self,x,y):
        i = self.rnd_struct.randrange(0,len(self.base_pairs))
        
        def f():
            q = None
            try:
                q = self.base_pairs[i](x,y)
            except:
                q = self.catch_value_error(x,y)
            ans = self.replace_inf_(q) 
            ans = ans % (DEFAULT_SINGLETON_RANGE[1] - DEFAULT_SINGLETON_RANGE[0])\
                + DEFAULT_SINGLETON_RANGE[0] 

            if type(self.bftd) != type(None):
                # append the ans to bftd 
                ans_ = [deepcopy(ans)]
                self.bftd.load_into_cycle_var(i,ans_)
                stat = self.bftd.prev_vstat[i]
                if self.bftd_tmap[i] == -1 and stat:
                    self.bftd_tmap[i] = self.bftd_count
            return ans

        self.bftd_count += 1        
        return f()

    def t_stat(self):
        return self.bftd.prev_vstat

    @staticmethod
    def load_AltBaseFunc_function(abf):
        return abf.output

####################################################
####### default basic generator functions
####### for the classes and functions pertinent
####### to this file. 


"""
Generates a function F that takes a boolean sequence 
B as input and outputs its termination status.
Elements of B are True if they fail metrics of <APRNGGauge>.

seq_size := the number of False elements in a row; considers the 
            last [-seq_size:] subsequence of input to F. 
ratio := (optional) / ratio that must be satisfied
"""
def generate_term_cond_for_boolseq(seq_size=None,ratio:float=None):
    if type(ratio) != type(None):
        return lambda bseq: bseq[-seq_size:].count(True) == seq_size \
            and round(bseq[-seq_size:].count(True) / seq_size,5) >= ratio
    return lambda bseq: bseq[-seq_size:].count(True) == seq_size

def generate_cov_uwpd_score_to_bool_function(cov_min:float,uwpd_min:float):
    return lambda x: x[0] <= cov_min and x[1] <= uwpd_min 


def term_cond_for_boolseq_sample_1():
    seq_size = 5
    ratio = 0.75
    return generate_term_cond_for_boolseq(seq_size,ratio)

def cov_uwpd_score_to_bool_function_sample_1():
    cov_min = 0.65
    uwpd_min = 0.10
    return generate_cov_uwpd_score_to_bool_function(cov_min,uwpd_min)

# 0.12, 0.22 

def default_BaseFuncTDetectDB(bps,check_size=5):
    
    bftd = BaseFuncTDetectDB(bps,\
        tcond1=term_cond_for_boolseq_sample_1(),\
        tcond2=DEFAULT_TERMINATE_DETECT_FUNC,\
        score2bool_func=cov_uwpd_score_to_bool_function_sample_1(),\
        check_size=check_size)
    return bftd 
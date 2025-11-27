from base.custom_rnd import * 
import unittest

### lone file test 
"""
python3 -m tests.test_custom_rnd
"""
###

"""
iws := sequence comprised of (iteration size, wanted stat:: bool)
"""
def assert_bftd_stat(seq,bftd,iterations,bfunc_index,iws): 
    
    def q(iterations,stat):
        for _ in range(iterations):
            bftd.load_into_cycle_var(bfunc_index,deepcopy(seq))
            assert bftd.prev_vstat[bfunc_index] == stat
        return
    
    inw_index = 0
    index = 0

    while index < iterations and inw_index < len(iws):
        s = iws[inw_index]
        s[0] = min([s[0], iterations - index])
        q(s[0] ,s[1])
        index = index + s[0] 
        inw_index += 1

class BaseFuncTDetectDBClass(unittest.TestCase):

    def test__BaseFuncTDetectDB_load_into_cycle_var__case1(self):
        # case 1
        bftd = default_BaseFuncTDetectDB(5)
        qx = [0.2 if i % 2 else 0.21 for i in range(5)]

        iterations = 12
        bfunc_index = 1
        iws = [[4,False],[8,True]]
        assert_bftd_stat(qx,bftd,iterations,bfunc_index,iws)


        # case 2 
        bftd = default_BaseFuncTDetectDB(5)
        qx = [1,1,1,1] + [1]

        iterations = 12
        bfunc_index = 1
        iws = [[4,False],[8,True]]
        assert_bftd_stat(qx,bftd,iterations,bfunc_index,iws)


        # case 3
        bftd = default_BaseFuncTDetectDB(5)
        qx = [i % 2 for i in range(5)] 

        iterations = 12
        bfunc_index = 1
        iws = [[12,False]]
        assert_bftd_stat(qx,bftd,iterations,bfunc_index,iws)

        # case 4 
        bftd = default_BaseFuncTDetectDB(5)
        qx = [0.2 if i % 2 else 0.35 for i in range(5)]

        iterations = 12
        bfunc_index = 1
        iws = [[4,False],[8,True]]
        assert_bftd_stat(qx,bftd,iterations,bfunc_index,iws)

        # case 5 
        bftd = default_BaseFuncTDetectDB(5)
        qx = [0.2 if i % 2 else 0.40 for i in range(5)]

        iterations = 12
        bfunc_index = 1
        iws = [[12,False]]
        assert_bftd_stat(qx,bftd,iterations,bfunc_index,iws)

        assert True         
        return

class AltBaseFuncClass(unittest.TestCase):

    def test__AltBaseFunc_load_into_cycle_var__case1(self):
        xx = [[1,2],[1,3]]
        base_pairs = [DEFAULT_PAIRWISE_VEC_FUNCS[3]]
        abf = AltBaseFunc(base_pairs,random)
        abf.attach_base_func_tdetect_db(bps=len(base_pairs))

        for i in range(50):
            qi = xx[i % 2] 
            abf.output(qi[0],qi[1])
            if abf.bftd_count >= 25:
                assert abf.t_stat() == [True]

        assert abf.bftd_tmap == {0:25}

    def test__AltBaseFunc_load_into_cycle_var__case2(self):
        xx = [[0,0],[1,1]]
        #xx = [[1,2],[1,3]]

        q = generate_efunc_type_q(5,10)#,idn_fx,sqrt_fx)
        q2 = generate_efunc_type_q(15,1.4)#,idn_fx,sqrt_fx)
        base_pairs = [DEFAULT_PAIRWISE_VEC_FUNCS[0],\
                DEFAULT_PAIRWISE_VEC_FUNCS[1],\
                DEFAULT_PAIRWISE_VEC_FUNCS[3],\
                q,\
                q2]

        # case 1                 
        abf = AltBaseFunc(base_pairs,random)
        abf.attach_base_func_tdetect_db(bps=len(base_pairs))

        for i in range(500):
            qi = xx[i % 2] 
            abf.output(qi[0],qi[1])
        ltm = list(abf.bftd_tmap.values())
        assert min(ltm) != -1 and max(ltm) < 200, "got {}".format(abf.bftd_tmap)

        # case 2
        xx = [[1,2],[1,3]]
        abf = AltBaseFunc(base_pairs,random)
        abf.attach_base_func_tdetect_db(bps=len(base_pairs))

        for i in range(500):
            qi = xx[i % 2] 
            abf.output(qi[0],qi[1])
        assert abf.bftd_tmap[4] == 278, "got {}, want {}".format(\
            abf.bftd_tmap[4],278)
        
if __name__ == '__main__':
    unittest.main()
from secnet import * 
import unittest

### lone file test 
"""
python3 -m tests.test_srefmap
"""
###

class SRefMapClass(unittest.TestCase):

    def test__SRefMap__fc_proc__best_nodedec_map__case1(self):
        srm = SRefMap_sample1() 
        nm = srm.fc_proc__best_nodedec_map(min,[0,1])
        assert len(nm) == 12

    def test__SRefMap__prmap_for_nodedec__case1(self):
        sn = SecNet_sample1()

        srm = sn.srm
        srm.reprocess('d')

        opm = srm.opmn[0]
        prmap = srm.prmap_for_nodedec(0,6,0,'greedy-d')

        assert opm[0] > prmap[0]
        assert prmap[6] > opm[6] 

        prmap = srm.prmap_for_nodedec(0,0,0,'greedy-d')
        assert opm[0] == prmap[0]
        assert prmap[6] == opm[6] 

    def test__SRefMap__collect_prism_points__DecMap__case1(self):
        ## use <SRefMap> on sample 3
        sn = SecNet_sample1() 
        srm = sn.srm

        dm = srm.collect_prism_points__DecMap('c',max,[0])  
        dm2 = srm.collect_prism_points__DecMap('d',max,[0])  
        dm3 = srm.collect_prism_points__DecMap('cd',max,[0])  

        dm4 = srm.collect_prism_points__DecMap('c',max,[1])
        dm5 = srm.collect_prism_points__DecMap('d',max,[1])
        dm6 = srm.collect_prism_points__DecMap('cd',max,[1])

        assert dm != dm2 
        #assert dm == dm3, "GOT \n{}\n\n{}".format(dm,dm3)
        assert dm4 != dm5
        #assert dm4 == dm6
        assert dm2 == dm5
        return
    
    # NOTE: dummy test; checks for crash-free exec.
    def test__SRefMap__collect_prism_points__PrMap__case1(self):
        sn = SecNet_sample1() 
        srm = sn.srm
        
        dx = srm.collect_prism_points__PrMap('c',"greedy-lone",0)
        dx2 = srm.collect_prism_points__PrMap('c',"greedy-d",0)
        dx3 = srm.collect_prism_points__PrMap('c',"greedy-c",0)
        assert True 

    # NOTE: shallow test; checks only for correct 
    #       number of prism values for each
    def test__SRefMap__build_prism__case1(self):
        sn = SecNet_sample1() 
        srm = sn.srm
        srm.build_prism() 

        assert len(srm.prism_typePr) == 24
        assert len(srm.prism_typeDec) == 18 

    def test__SRefMap__binarycmp_prism_points__typePr__case1(self):
        fobj = open("sampleX","rb")
        srm = pickle.load(fobj) 
        fobj.close()
        assert type(srm) == SRefMap 

        qx = list(srm.prism_typePr.keys())
        prkey1 = qx[0] 
        prkey2 = qx[0] 

        prmap = srm.binarycmp_prism_points__typePr(prkey1,\
                prkey2,0,0,5)

        assert sum(prmap.values()) == 0.0 

    def test__SRefMap__binarycmp_prism_points__typePr__case2(self):
        fobj = open("sampleX","rb")
        srm = pickle.load(fobj) 
        fobj.close()
        assert type(srm) == SRefMap 

        sec_idn = 0
        num_opt = len(srm.opmn[sec_idn])

        prkey1 = 'cd,greedy-dc,0'
        cx_,cx = 0,0
        for i in range(num_opt - 1): 
            for j in range(i+1,num_opt): 
                prmap_ = srm.binarycmp_prism_points__typePr(prkey1,\
                        prkey1,sec_idn,i,j)
                #print("sum for {},{}: {}".format(i,j,sum(prmap_.values())))
                s = sum(prmap_.values())
                if s != 0.0:
                    cx_ += 1
                cx += 1

        assert cx_ == 51 and cx == 66 

    # NOTE: dummy test; does not check for correct sol'n, only 
    #       their existence
    def test__SRefMap__binarycmp_prism_points__typeDec__case1(self):
        fobj = open("sampleX","rb")
        srm = pickle.load(fobj) 
        assert type(srm) == SRefMap 
        fobj.close()

        ## demonstrate the Pr-prism. 
        qx = list(srm.prism_typeDec.keys())
        for i in range(len(qx) - 1):
            for j in range(i+1,len(qx)):
                bpp = srm.binarycmp_prism_points__typeDec(qx[i],qx[j])
                assert len(bpp) > 0 

    def test__SRefMap__cd_density_measure__case1(self):
        fobj = open("sampleX","rb")
        srm = pickle.load(fobj) 
        assert type(srm) == SRefMap 
        fobj.close()

        srm.cd_density_measure(True)
        srm.cd_density_measure(False)

        ##### measure density of dep.

        # count the number of non-null
        def count_active(d):
                return len([v_ for v_ in d.values() if v_[0] != 0.0])

        keys = srm.density_cd.keys()
        dx = {}
        for k in keys:
                c = count_active(srm.density_cd[k]) 
                dx[k] = c

        ans = {0: 1,\
                1: 2,\
                2: 1,\
                3: 0,\
                4: 3,\
                5: 1,\
                6: 2,\
                7: 2,\
                8: 1,\
                9: 0,\
                10: 2,\
                11: 0}
        assert dx == ans 

        ##### measure density of co-dep.

        ql = []
        for k,v in srm.density_d.items():
                for k2,v2 in v.items():
                        ##print("{}:{}:{}".format(k,k2,v2[0]))
                        #assert v2[0] >= 0.63
                        ql.append(v2[0])

        assert min(ql) > 0.09
        assert max(ql) == 1.0 

        mn = sum(ql) / len(ql)
        assert round(abs(mn - 0.7929292929292922),12) == 0.0 


if __name__ == '__main__':
    unittest.main() 
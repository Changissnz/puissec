from secnet_gen import * 
import unittest

### lone file test 
"""
python3 -m tests.test_secnet_gen
"""
###

def Sec_list_sample_1():

    singleton_range = [0.,10.] 
    dimension = 5
    num_optima = 2
    countermeasure = (0.6,0.5) 

    sec = Sec.generate_bare_instance(singleton_range,dimension,num_optima,\
            countermeasure,rnd_struct=np.random)

    sec2 = Sec.generate_bare_instance(singleton_range,dimension,num_optima,\
            countermeasure,rnd_struct=np.random)

    sec3 = Sec.generate_bare_instance(singleton_range,dimension,num_optima,\
            countermeasure,rnd_struct=np.random)

    sec4 = Sec.generate_bare_instance(singleton_range,dimension,num_optima,\
            countermeasure,rnd_struct=np.random)

    return [sec,sec2,sec3,sec4]

class SecNetGenSchemeClass(unittest.TestCase):

    # NOTE: weak assertion in testing
    def test__SecNetGenScheme_next_conn__spineframe(self):
        snv = [0,4,10,12,24]
        snv = [0,4,10,12,24,33,56,112,200]
        nsnv = [1,2,3,5]
        sngs = SecNetGenScheme("spine frame",5)
        snfg = SecNetFrameGen(snv,nsnv,sngs)

        snfg.construct_frame() 
        ##
        """
        print(snfg.node_components)
        print()
        print(snfg.d)
        print("distances")
        for k,v in snfg.node2sec_distance_map.items():
            print("K: ",k)
            print(v)
            print() 
        """
        ##
        assert len(snfg.node_components) == 1 

    def test__SecNetGenScheme_next_conn__pairingframe(self):

        snv = [0,4,10,12,24]
        nsnv = [1,2,3,5]
        sngs = SecNetGenScheme("pairing frame",5)
        snfg = SecNetFrameGen(snv,nsnv,sngs)

        # make the first four connections
        # check each connection b/t a secnode and nsecnode
        for i in range(4):
            snfg.next_conn() 

        for n in snv:
            ci = snfg.node_to_component_index(n)
            q = snfg.node_components[ci]

            stat1 = len(q) == 1
            stat2 = False
            for x in q: 
                if x in snfg.nsec_nodevec:
                    stat2 = True 
                    break 
            assert stat1 or stat2 

        snfg.construct_frame()
        assert len(snfg.node_components) == 1 
        return

    def test__SecNetGenScheme_next_conn__pseudorandom(self):

        snv = [0,4,10,12,24,33,56,112,200]
        nsnv = [1,2,3,5]
        sngs = SecNetGenScheme("pseudo random",5,0.5)
        snfg = SecNetFrameGen(snv,nsnv,sngs)
        snfg.construct_frame() 
        v = snfg.node2sec_distance_map[1] 
        d = defaultdict(int,{3: 1, 4: 1, 0: 1, 33: 1, 1: 0, 2: 2, 5: 2,\
            200: 2, 10: 2, 12: 2, 112: 2, 24: 2, 56: 2})
        assert v == d 

    def test__SecNetGenScheme_generate_graph__type_prop(self):

        irc_sz = 13
        p_s,p_n = 1.0,1.3

        g = SecNetFrameGen.generate_graph__type_prop(irc_sz,p_s,p_n,random)

        assert len(g[0]) >= int(round(13 * 2.3)) 
        assert len(g[1]) == 13 

        for k,v in g[0].items():
            assert len(v) > 0 



class SecNetDepGenClass(unittest.TestCase):

    def test__SecNetDepGen__make_conn__case1(self):

        secs = Sec_list_sample_1()
        sndg = SecNetDepGen(secs,random,2,0.75,[1,4])
        stat = True
        i = 0 
        while stat and i < 100: 
            stat = sndg.make_conn([2,3])
            i += 1 
            '''
            print("stat: ",stat)
            print("dep map")
            print(sndg.dep_map)
            print("---------------------")
            '''

        #print(i)
        assert i <= 25

    def test__SecNetDepGen__make_conn__case2(self):
        secs = Sec_list_sample_1()
        sndg = SecNetDepGen(secs,random,2,0.75,[1,4])

        stat = True
        i = 0 
        #for i in range(10):
        #    print("ITER: ",i)
        while stat and i < 100: 
            stat = sndg.make_conn([1])
            i += 1 
            '''
            print("stat: ",stat)
            print("dep map")
            print(sndg.dep_map)
            print("---------------------")
            '''
        #print(i)
        assert i <= 7, "want {},got {}".format(5,i)

    def test__SecNetDepGen__make_conn__case3(self):

        s = Sec_list_sample2(num_secs=4)
        sndg = SecNetDepGen(s,random,4,0.0,[1,4],depconn_ratio=0.5)
        sndg.assign_conn(500)
        ss = SecSeq(sndg.sq)

        ans = {0:(31,0),1:(64,0),2:(0,0),3:(32,0)}
        for s in ss.sequence:
            q = ans[s.idn_tag] 
            assert q[0] == len(s.dm)
            assert q[1] == len(s.cdm) 

if __name__ == '__main__':
    unittest.main()
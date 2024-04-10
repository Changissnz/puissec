from secnet_gen import * 
import unittest

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

if __name__ == '__main__':
    unittest.main()
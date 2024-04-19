from pr_catgen import * 
import unittest

class PRCatGenFunctionsClass(unittest.TestCase):

    def test__exact_correlation_dep_Pr(self):
    
        random.seed(3)
        d = {0:1,1:3,2:5,3:0,\
            4:1,5:2,6:4,7:2,8:3,\
            9:0,10:2,11:0}
        pred_corrmap = defaultdict(int,d)

        dim = (12,4)
        pr_dim = generate_2d_match_pr(dim)

        pred_opt2pr_map = {0:0.6,1:0.1,2:0.1,3:0.,4:0.2}
        q = exact_correlation_dep_Pr(12,pred_corrmap,pr_dim)#pred_opt2pr_map)

        ans1 = defaultdict(float,\
            {'0,1': 0.71882, '0,2': 0.0,\
            '0,3': 0.75823, '0,4': 0.0,\
            '0,5': 0.0, '0,6': 0.0,\
            '0,7': 0.0, '0,8': 0.71882,\
            '0,9': 0.75823, '0,10': 0.0,\
            '0,11': 0.75823, '1,2': 0.0,\
            '1,3': 0.30127, '1,4': 0.71882,\
            '1,5': 0.87887, '1,6': 0.0,\
            '1,7': 0.87887, '1,8': 0.0,\
            '1,9': 0.30127, '1,10': 0.87887,\
            '1,11': 0.30127,'2,3': 0.0,\
            '2,4': 0.0, '2,5': 0.0,\
            '2,6': 0.0, '2,7': 0.0,\
            '2,8': 0.0, '2,9': 0.0,\
            '2,10': 0.0, '2,11': 0.0,\
            '3,4': 0.75823, '3,5': 0.47275,\
            '3,6': 0.0, '3,7': 0.47275,\
            '3,8': 0.30127, '3,9': 0.8849,\
            '3,10': 0.47275, '3,11': 0.8849,\
            '4,5': 0.0, '4,6': 0.0,\
            '4,7': 0.0, '4,8': 0.71882,\
            '4,9': 0.75823, '4,10': 0.0,\
            '4,11': 0.75823, '5,6': 0.0,\
            '5,7': 0.9211, '5,8': 0.87887,\
            '5,9': 0.47275, '5,10': 0.9211,\
            '5,11': 0.47275, '6,7': 0.0,\
            '6,8': 0.0, '6,9': 0.0,\
            '6,10': 0.0, '6,11': 0.0,\
            '7,8': 0.87887, '7,9': 0.47275,\
            '7,10': 0.9211, '7,11': 0.47275,\
            '8,9': 0.30127, '8,10': 0.87887,\
            '8,11': 0.30127, '9,10': 0.47275,\
            '9,11': 0.8849, '10,11': 0.47275})

        assert q == ans1 

    def test__partial_correlation_dep_Pr(self):
        random.seed(3)
        d = {0:1,1:3,2:5,3:0,\
        4:1,5:2,6:4,7:2,8:3,\
        9:0,10:2,11:0}
        pred_corrmap = defaultdict(int,d)

        dim = (12,4)
        pr_dim = generate_2d_match_pr(dim)

        pred_opt2pr_map = {0:0.6,1:0.1,2:0.1,3:0.,4:0.2}

        random.seed(8)
        pr_dist = generate_pr_dist_D1_to_D2(5,12)

        prx = partial_correlation_dep_Pr(12,pred_corrmap,\
            pr_dist,pred_opt2pr_map)

        ansd = defaultdict(float, {\
            '0,0': 0.0056, '0,1': 0.00535, '0,2': 0.00603,\
            '0,3': 0.00521, '0,4': 0.0056, '0,5': 0.0049,\
            '0,6': 0.00342, '0,7': 0.0049, '0,8': 0.00535,\
            '0,9': 0.00521, '0,10': 0.0049, '0,11': 0.00521,\
            '1,0': 0.00535, '1,1': 0.0121, '1,2': 0.00936,\
            '1,3': 0.00273, '1,4': 0.00535, '1,5': 0.00456,\
            '1,6': 0.00434, '1,7': 0.00456, '1,8': 0.0121,\
            '1,9': 0.00273, '1,10': 0.00456, '1,11': 0.00273,\
            '2,0': 0.00603, '2,1': 0.00936, '2,2': 0.01064,\
            '2,3': 0.00498, '2,4': 0.00603, '2,5': 0.00767,\
            '2,6': 0.00705, '2,7': 0.00767, '2,8': 0.00936,\
            '2,9': 0.00498, '2,10': 0.00767, '2,11': 0.00498,\
            '3,0': 0.00521, '3,1': 0.00273, '3,2': 0.00498,\
            '3,3': 0.00571, '3,4': 0.00521, '3,5': 0.00522,\
            '3,6': 0.00352, '3,7': 0.00522, '3,8': 0.00273,\
            '3,9': 0.00571, '3,10': 0.00522, '3,11': 0.00571,\
            '4,0': 0.0056, '4,1': 0.00535, '4,2': 0.00603,\
            '4,3': 0.00521, '4,4': 0.0056, '4,5': 0.0049,\
            '4,6': 0.00342, '4,7': 0.0049, '4,8': 0.00535,\
            '4,9': 0.00521, '4,10': 0.0049, '4,11': 0.00521,\
            '5,0': 0.0049, '5,1': 0.00456, '5,2': 0.00767,\
            '5,3': 0.00522, '5,4': 0.0049, '5,5': 0.00721,\
            '5,6': 0.00617, '5,7': 0.00721, '5,8': 0.00456,\
            '5,9': 0.00522, '5,10': 0.00721, '5,11': 0.00522,\
            '6,0': 0.00342, '6,1': 0.00434, '6,2': 0.00705,\
            '6,3': 0.00352, '6,4': 0.00342, '6,5': 0.00617,\
            '6,6': 0.00583, '6,7': 0.00617, '6,8': 0.00434,\
            '6,9': 0.00352, '6,10': 0.00617, '6,11': 0.00352,\
            '7,0': 0.0049, '7,1': 0.00456, '7,2': 0.00767,\
            '7,3': 0.00522, '7,4': 0.0049, '7,5': 0.00721,\
            '7,6': 0.00617, '7,7': 0.00721, '7,8': 0.00456,\
            '7,9': 0.00522, '7,10': 0.00721, '7,11': 0.00522,\
            '8,0': 0.00535, '8,1': 0.0121, '8,2': 0.00936,\
            '8,3': 0.00273, '8,4': 0.00535, '8,5': 0.00456,\
            '8,6': 0.00434, '8,7': 0.00456, '8,8': 0.0121,\
            '8,9': 0.00273, '8,10': 0.00456, '8,11': 0.00273,\
            '9,0': 0.00521, '9,1': 0.00273, '9,2': 0.00498,\
            '9,3': 0.00571, '9,4': 0.00521, '9,5': 0.00522,\
            '9,6': 0.00352, '9,7': 0.00522, '9,8': 0.00273,\
            '9,9': 0.00571, '9,10': 0.00522, '9,11': 0.00571,\
            '10,0': 0.0049, '10,1': 0.00456, '10,2': 0.00767,\
            '10,3': 0.00522, '10,4': 0.0049, '10,5': 0.00721,\
            '10,6': 0.00617, '10,7': 0.00721, '10,8': 0.00456,\
            '10,9': 0.00522, '10,10': 0.00721, '10,11': 0.00522,\
            '11,0': 0.00521, '11,1': 0.00273, '11,2': 0.00498,\
            '11,3': 0.00571, '11,4': 0.00521, '11,5': 0.00522,\
            '11,6': 0.00352, '11,7': 0.00522, '11,8': 0.00273,\
            '11,9': 0.00571, '11,10': 0.00522, '11,11': 0.00571})

        assert prx == ansd

    def test__generate_pr_dist_D1_to_D2(self):
        random.seed(8)

        pr_dist = generate_pr_dist_D1_to_D2(5,12)

        ans0 = defaultdict(float,\
        {0: 0.10768, 1: 0.45707,\
        2: 0.06001, 3: 0.33478, 4: 0.04046})
        assert pr_dist[0] == ans0 
        return

    def test__generate_bounds_vector_sequence(self):
        np.random.seed(30)

        superbound = np.array([[0,1.],\
                    [0,1.],\
                    [0,1.]]) 

        spacing_ratio_range = [0.,0.1]
        outlier_pr_ratio = 0.1#1.0
        num_bounds = 3

        bvseq = generate_bounds_vector_sequence(superbound,\
            spacing_ratio_range,outlier_pr_ratio,num_bounds)

        for i in range(1,len(bvseq)):

            b1 = bvseq[i -1]
            b2 = bvseq[i]
            assert np.all(b2[:,0] - b1[:,1])
        return

if __name__ == '__main__':
    unittest.main()
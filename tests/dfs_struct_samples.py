from dfs_struct import *

"""

    0---1
    |\ 
    | \ 
    2  \ 
        3 

"""
def test_dfs_graph_1():
    d = defaultdict(set) 
    d[0] = {1,2,3}
    d[1] = {0}
    d[2] = {0}
    d[3] = {0}
    return d

"""

                _______
     ____  ____/_______\ 
    /    \/   /   \     \ 
    0--1--2--3--4--5--6--7--8
"""
def test_dfs_graph_2():
    d = defaultdict(set) 
    d[0] = {1,2}
    d[1] = {0,2}
    d[2] = {0,1,3,5,7}
    d[3] = {2,4,7}
    d[4] = {3,5}
    d[5] = {2,4,6}
    d[6] = {5,7}
    d[7] = {2,3,6,8}
    d[8] = {7}
    return d

"""

    0---1
    | \/| 
    | /\|
    2---3

"""
def test_dfs_graph_3():
    d = defaultdict(set)
    d[0] = {1,2,3}
    d[1] = {0,2,3}
    d[2] = {0,1,3}
    d[3] = {0,1,2}
    return d

"""

    0--1--2--3--4
"""
def test_dfs_graph_4():
    d = defaultdict(set)
    d[0] = {1}
    d[1] = {0,2}
    d[2] = {1,3}
    d[3] = {2,4}
    d[4] = {3}
    return d
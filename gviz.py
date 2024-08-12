# some basic methods for visualizing a graph 
from defaults import * 
import networkx as nx 
import matplotlib.pyplot as plt

def graph_to_viz(g,save_fig:str=""):
    assert type(g) in {defaultdict,dict}
    G = nx.Graph()

    lk = list(g.keys())
    ex = []
    for k,v in g.items():
        for v_ in v:
            ex.append((k,v_))

    G.add_nodes_from(lk)
    G.add_edges_from(ex) 

    nx.draw(G)

    if save_fig != "":
        plt.savefig(save_fig)
    plt.show()
    
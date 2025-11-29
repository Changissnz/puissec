# some basic methods for visualizing a graph 
from base.defaults import * 
import networkx as nx 
import matplotlib.pyplot as plt

def default_Puissec_node_color_map(g,secnodes,cr_loc,ir_loc):
    color_map = [] 
    
    nodes = list(g.keys())

    for n in nodes: 
        # case: Crackling and IsoRing are on the same node 
        if n in cr_loc and n in ir_loc: 
            x2 = (n, {"color": "black"}) 
        # case: Crackling location 
        elif n in cr_loc: 
            x2 = (n,{"color":"green"}) 
        # case: IsoRing location 
        elif n in ir_loc: 
            x2 = (n,{"color":"orange"}) 
        # case: Sec node 
        elif n in secnodes: 
            x2 = (n,{"color":"red"}) 
        else: 
            x2 = (n,{"color":"blue"}) 
        color_map.append(x2) 
    return color_map 

def SecEnv_data_to_viz(g,secnodes,cr_loc=set(),ir_loc=set(),save_fig:str="",no_show:bool=False):
    assert type(g) in {defaultdict,dict}

    color_map = default_Puissec_node_color_map(g,secnodes,cr_loc,ir_loc)
    G = nx.Graph() 

    nodes = list(g.keys())
    ex = []
    for k,v in g.items():
        for v_ in v:
            ex.append((k,v_))

    G.add_nodes_from(nodes)
    G.add_edges_from(ex) 
    pos = nx.spring_layout(G) 

    ncs = [cm[1]["color"] for cm in color_map]
    lsx = dict([(cm[0],str(cm[0])) for cm in color_map]) 

    nx.draw_networkx_nodes(G, pos, node_color = ncs)
    nx.draw_networkx_labels(G, pos, labels = lsx, font_size = 12)
    nx.draw_networkx_edges(G, pos, edge_color = 'black')

    if save_fig != "":
        plt.savefig(save_fig)

    if no_show: 
        return G 

    plt.show()
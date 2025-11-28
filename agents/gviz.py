# some basic methods for visualizing a graph 
from base.defaults import * 
import networkx as nx 
import matplotlib.pyplot as plt

def get_node_colors(G):
    colors = []
    for i in range(len(G.nodes)):
        try:
            colors.append(f"tab:{G.nodes[i]['data']['color']}")
        except KeyError as e:
            colors.append('tab:red')
    print()
    print('Node colors:')
    print(colors)
    return colors   

def SecNet_graph_to_viz(g,secnodes,save_fig:str=""):
    assert type(g) in {defaultdict,dict}
    G = nx.Graph()

    lk = list(g.keys())
    lk_ = []

    for lk2 in lk:
        cx = "red" if lk2 in secnodes else "blue"
        x2 = (lk2, {"color": cx})
        print("LK2: ",x2)
        lk_.append(x2)
    lk = lk_ 

    print("LK")
    print(lk)

    ex = []
    for k,v in g.items():
        for v_ in v:
            ex.append((k,v_))

    G.add_nodes_from(lk)
    G.add_edges_from(ex) 
    pos = nx.spring_layout(G)

    ncs = [lk__[1]["color"] for lk__ in lk]
    lsx = dict([(lk__[0],str(lk__[0])) for lk__ in lk]) 

    nx.draw_networkx_nodes(G, pos, node_color = ncs)
    nx.draw_networkx_labels(G, pos, labels = lsx, font_size = 12)
    nx.draw_networkx_edges(G, pos, edge_color = 'black')
    
    if save_fig != "":
        plt.savefig(save_fig)
    plt.show()

    
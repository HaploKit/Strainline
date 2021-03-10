import sys
import networkx as nx


def is_connected(gfa_file):
    links = []
    with open(gfa_file) as  fr:
        for line in fr:
            if line.startswith('S'):
                pass
            elif line.startswith('L'):
                a = line.split()
                links.append((a[1], a[3]))
    G = nx.Graph()
    G.add_edges_from(links)
    print('Number connected components:{}'.format(nx.number_connected_components(G)))
    print('Number of nodes: {}'.format(G.number_of_nodes()))
    print('The size of connected components is:{}'.
          format([len(c) for c in sorted(nx.connected_components(G), key=len, reverse=True)]))
    return nx.is_connected(G)

# gfa='/Users/xiaoluo/Documents/CWI/project/vg/bandage/out.pruned2.compressed.gfa'
gfa,=sys.argv[1:]
print("The graph is connected: {}".format(is_connected(gfa)))

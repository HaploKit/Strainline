import networkx as nx
import sys

dot, gfa, min_cov = sys.argv[1:]

min_cov = int(min_cov)
G = nx.DiGraph(nx.drawing.nx_pydot.read_dot(dot))

print('reading dot is finished...')

fw = open(gfa, 'w')
fw.write('H\tVN:Z:1.0\n')

link_infos = []
reliable_nodes = {}
reliable_edges =[]
for edge, label in nx.get_edge_attributes(G, 'label').items():
    cov = int(label.strip('"'))/2 #weight of edge = 2*cov
    # print('cov:{}'.format(cov))
    if cov >= min_cov:
        link_infos.append('L\t' + edge[0] + '\t+\t' + edge[1] + '\t+\t' + '0M\n')
        reliable_edges.append((edge[0],edge[1]))
        reliable_nodes[edge[0]] = 1
        reliable_nodes[edge[1]] = 1

G0 = nx.Graph()
G0.add_edges_from(reliable_edges)
print('The size of connected components is:{}'.
      format([len(c) for c in sorted(nx.connected_components(G0), key=len, reverse=True)]))
G1 = G0.subgraph(sorted(nx.connected_components(G0), key=len, reverse=True)[0])
print('The largest component is connected:{}, size:{}'.format(nx.is_connected(G1), len(G1)))

G1_nodes = {v for v in G1.nodes()}
G1_edges = {(e1, e2) for e1, e2 in G1.edges()}


#only output the largest component
for node, label in nx.get_node_attributes(G, 'label').items():
    if node not in G1_nodes:
        continue
    seq = label.strip('"').split()[-1]
    fw.write('S\t' + node + '\t' + seq + '\n')

for link_info in link_infos:
    node1, node2 = link_info.split()[1],link_info.split()[3]
    if ((node1, node2) in G1.edges) or ((node2, node1) in G1_edges):
        fw.write(link_info)

fw.close()

## View dot file
# from graphviz import Source
# file = open('xx/reads.top10.sorted.poa.dot', 'r')#READING DOT FILE
# text=file.read()
# Source(text) #very slow for large graph

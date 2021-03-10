import networkx as nx
import sys
import copy
import numpy as np


def remove_edges(graph, link2cov, min_ratio, min_cov, mid_min_cov, end_cutoff, max_node_id):
    '''
    prune edges in the De Bruijn graph by min coverage and its relative ratio
    :param graph: input graph
    :param link2cov: {(node1,node2):2,(node4,node8):4}
    :param min_ratio: remove this edge if its coverage ratio < min_ratio
    :param min_cov: remove this edge if its coverage < min_cov
    :param mid_min_cov: remove this edge (not in both ends of graph) if its coverage < mid_min_cov
    :param end_cutoff: define how many nodes in the end belong to 'end' (for VG, not DBG)
    :param max_node_id: the max node id in its corresponding variation graph
    :return: pruned graph
    '''
    G = copy.deepcopy(graph)
    del_links = {}
    # num_nodes = len(G)
    for node in G.nodes():  # node: 12+
        # remove bad in-edges
        cov_hash = {}
        for pre_node in G.predecessors(node):
            if (pre_node, node) in link2cov:
                cov_hash[pre_node] = link2cov[(pre_node, node)]
            else:
                cov_hash[pre_node] = 0
        cov_sum = sum(cov_hash.values())
        for pre_node, cov in cov_hash.items():
            if (cov_sum == 0) or (cov * 1.0 / cov_sum < min_ratio) or (cov < min_cov):
                del_links[(pre_node, node)] = 1
            elif (int(node.split('+')[0]) > end_cutoff) or (int(node.split('+')[-1]) < (max_node_id - end_cutoff)):
                # use more strigent threshold in the middle of genome
                # nodes range from end_node_len to (num_nodes-end_node_len) will be regarded as middle nodes
                if cov < mid_min_cov:
                    del_links[(pre_node, node)] = 1

        # remove bad out-edges
        cov_hash = {}
        for suc_node in G.successors(node):
            if (node, suc_node) in link2cov:
                cov_hash[suc_node] = link2cov[(node, suc_node)]
            else:
                cov_hash[suc_node] = 0
        cov_sum = sum(cov_hash.values())
        for suc_node, cov in cov_hash.items():
            if (cov_sum == 0) or (cov * 1.0 / cov_sum < min_ratio) or (cov < min_cov):
                del_links[(node, suc_node)] = 1
            elif (int(node.split('+')[0]) > end_cutoff) or (int(node.split('+')[-1]) < (max_node_id - end_cutoff)):
                # use more strigent threshold in the middle of genome
                if cov < mid_min_cov:
                    del_links[(node, suc_node)] = 1

    G.remove_edges_from(del_links.keys())
    return G


def get_candidate_paths(vg_gfa, min_kmer_count, min_ratio, k=3, end_cutoff=0, min_path_len=0):
    '''construct a general De Bruijn graph from the read variation graph,
    node in dbg is kmer which is the concatination of k nodes in vg: 1+2+4,
    there is an edge between two consecutive kmers,
    finally, enumerate all candidate paths which are supposed to be haplotypes.
    '''
    kmer2count = {}
    kmer_links = {}  # {('1+2+4','2+4+5'):1}
    with open(vg_gfa) as  fr:
        for line in fr:
            if line.startswith('P'):
                a = line.split()
                nodes = a[2].split(',')
                if nodes[0][-1] == '-':
                    raise Exception("Check path direction, only + is allowed.")
                pre_kmer = ''
                for i in range(len(nodes) - k):
                    kmer = ''.join(nodes[i:(i + k)])[:-1]  # 1+2+4
                    if kmer in kmer2count:
                        kmer2count[kmer] += 1
                    else:
                        kmer2count[kmer] = 1
                    if i == 0:
                        pre_kmer = kmer
                    else:
                        if (pre_kmer, kmer) in kmer_links:
                            kmer_links[(pre_kmer, kmer)] += 1  # coverage of edges
                        else:
                            kmer_links[(pre_kmer, kmer)] = 1
                        pre_kmer = kmer
            else:
                continue
    with open('kmer.txt', 'w') as fw:
        fw.write('\n'.join([k + '\t' + str(v) for k, v in kmer2count.items()]))

    bad_kmers = [kmer for kmer in kmer2count if kmer2count[kmer] < min_kmer_count]
    print('\nNumber of all kmers: {}'.format(len(kmer2count)))
    print('\nNumber of removed kmers: {}'.format(len(bad_kmers)))

    max_node_id = max([max(np.array(kmer.split('+'), dtype='int')) for kmer in kmer2count.keys()])
    print('\nThe max node id in the variation graph: {}\n'.format(max_node_id))

    G = nx.DiGraph()
    G.add_edges_from(kmer_links.keys())
    if len(list(nx.simple_cycles(G))) == 0:
        print('No cycle in De Bruijn graph.')
    else:
        print('Cycles are found in De Bruijn graph:{}'.format(len(list(nx.simple_cycles(G)))))

    # remove bad nodes
    G.remove_nodes_from(bad_kmers)

    # remove bad edges

    # min coverage for k+1 mer, which is the same with the edge weight of two nodes(kmer),
    # so use less strigent cutoff TODO
    min_cov = int(min_kmer_count * 0.9)

    mid_min_cov = 1.5 * min_kmer_count  # TODO
    G = remove_edges(G, kmer_links, min_ratio, min_cov, mid_min_cov, end_cutoff, max_node_id)

    componets = sorted(nx.connected_components(nx.to_undirected(G)), key=len, reverse=True)
    print('\nNumber connected components:{}'.format(len(componets)))
    print('\nNumber of nodes: {}'.format(G.number_of_nodes()))
    print('\nThe size of connected components is:{}'.format([len(c) for c in componets]))

    G = G.subgraph(componets[0])  # only keep the largest component
    print("\nNumber of edges in the largest component:{}\n".format(len(G.edges())))

    # remove tips? TODO

    if len(list(nx.simple_cycles(G))) == 0:
        print('No cycle in pruned De Bruijn graph.')
    else:
        print('Cycles are found in pruned De Bruijn graph:{}'.format(len(list(nx.simple_cycles(G)))))


    # add a source and sink node? not necessary
    source_nodes = []
    sink_nodes = []
    num_tips = 0
    graph_size = len(G)
    for node in G.nodes():
        if G.in_degree(node) == 0:
            pos = int(node.split('+')[0])
            if pos < end_cutoff:
                source_nodes.append(node)
            else:
                num_tips += 1
        elif G.out_degree(node) == 0:
            pos = int(node.split('+')[-1])
            # print('pos:{}, max_node_id:{}'.format(pos, max_node_id))
            if pos > (max_node_id - end_cutoff):
                sink_nodes.append(node)
            else:
                num_tips += 1
    print('Number of tips in De Bruijn graph (end cutoff = {}): {}'.format(end_cutoff, num_tips))

    # write DBG to gfa
    with open('dbg.gfa', 'w') as fw:
        fw.write('H\tVN:Z:1.0\n')
        fw.write('\n'.join(['S\t' + v + '\tA' for v in G.nodes()]) + '\n')  # TODO
        fw.write('\n'.join(['L\t' + v1 + '\t+\t' + v2 + '\t+\t' + str(k - 1) + 'M' for v1, v2 in G.edges()]))

    # enumerate all simple paths from sources to sinks, which can be regarded as candidate haplotypes
    candidate_paths = []
    i = 0
    print('\n\n')
    print('source nodes:{}\n'.format(source_nodes))
    print('sink nodes:{}\n'.format(sink_nodes))
    # print("Number of path:{}".format(len(list(nx.all_simple_paths(G, source=source_nodes[0], target=sink_nodes[0])))))

    for source_node in source_nodes:
        print('\nprcoessing source node:{}\n'.format(source_node))
        for path in nx.all_simple_paths(G, source=source_node, target=sink_nodes):
            # joint_path = []
            # for i, kmer in enumerate(path):
            #     if i == 0:
            #         joint_path = kmer.split('+')
            #     else:
            #         joint_path.append(kmer.split('+')[-1])
            # if len(joint_path) >= min_path_len:
            #     candidate_paths.append('+'.join(joint_path))

            i += 1
            if i % 1000 == 0:
                print("The {} path...".format(i))
            candidate_paths.append(path)

        break

    print('Number of total candidate paths:{}'.format(len(candidate_paths)))
    with open('candidate_paths.txt', 'w') as f:
        # f.write('\n'.join(candidate_paths))
        f.write('\n'.join([','.join(path) for path in candidate_paths]))
    return candidate_paths


if __name__ == '__main__':
    vg_gfa, min_kmer_count, min_ratio, k, end_cutoff, min_path_len = sys.argv[1:]
    get_candidate_paths(vg_gfa, int(min_kmer_count), float(min_ratio), int(k), int(end_cutoff), int(min_path_len))

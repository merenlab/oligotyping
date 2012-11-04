#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (C) 2010 - 2011, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

import os
import sys
import math
import networkx as nx
import matplotlib.pyplot as plt
from networkx import graphviz_layout

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
from utils.utils import pretty_print

def topology_graph(topology_file, match_levels = False):
    G = nx.MultiDiGraph()
    
    datafile = open(topology_file)
   
    nodes = {}
    levels = []
    for line in [l.strip('\n') for l in datafile.readlines()]:
        node_id, node_size, node_parent, node_level, node_children = line.split('\t')
        nodes[node_id] = {'size': node_size, 'parent': node_parent, 'level': int(node_level),
                          'children': node_children.split(',') if node_children else [], 'type': 'node'}
        levels.append(int(node_level))
    
    # match levels
    if match_levels:
        max_level = max(levels)
        new_nodes = {}
        for node_id in nodes:
            node = nodes[node_id]
            if node['level'] < max_level and not node['children']:
                levels_to_cover = range(node['level'] + 1, max_level + 1)
                for level in levels_to_cover:
                    if levels_to_cover.index(level) == 0:
                        new_nodes[node_id + ':l%d' % level] = {'size': node['size'], 'parent': node_id, 
                                                               'level': level, 'children': [], 'type': 'spam'}
                    else:
                        new_nodes[node_id + ':l%d' % level] = {'size': node['size'], 'parent': node_id + ':l%d' % (level - 1),
                                                               'level': level, 'children': [], 'type': 'spam'}

        for node_id in new_nodes:
            nodes[node_id] = new_nodes[node_id]

    for node_id in nodes:
        node = nodes[node_id]
        
        if node['parent']:
            if node['type'] == 'node':
                label = node_id
                while 1:
                    if label[0] == '0':
                        label = label[1:]
                    else:
                        break
                G.add_edge(node_id, node['parent'], size = int(node['size']), label = label)
            else:
                G.add_edge(node_id, node['parent'], size = int(node['size']), label = '')
   
    for node_id in nodes['root']['children']:
        G.add_edge('root', node_id, size = int(nodes['root']['size']), label = 'root')
   

    

    return G


def topology(topology_file, output_file = None, title = None):
    G = topology_graph(topology_file)

    number_of_edges = G.number_of_edges()
    number_of_nodes = G.number_of_nodes()

    print("Loaded %d edges and %d nodes." % (number_of_edges, number_of_nodes))

    plt.figure(figsize=(16, 16))
    
    # use graphviz to find radial layout
    # twopi, gvcolor, wc, ccomps, tred, sccmap, fdp, circo, neato, acyclic, nop, gvpr, dot
    pos=nx.graphviz_layout(G, prog="twopi")

    # node size is proportional to number of reads went into it
    sizes = dict.fromkeys(G.nodes(), 0.0)
    for (u, v, d) in G.edges(data=True):
        sizes[u] = d['size']
    max_size = max(sizes.values())
    k = 10000.0 / max_size
    for node in sizes:
        sizes[node] = sizes[node] * k if sizes[node] * k > 500 else 500
 
    shapes = dict.fromkeys(G.nodes(), 0.0)
    for (u, v, d) in G.edges(data=True):
        shapes[u] = 'o' if d['size'] > 1 else ''
 
 
    # edge width, not in use at this moment
    edgewidth = []
    for (u, v, d) in G.edges(data = True):
        edgewidth.append(2) #len(G.get_edge_data(u,v)))


    network_nodes = nx.draw_networkx_nodes(G, pos, node_shape = 'o', node_size = [sizes[i] for i in G], node_color=['#FFFFFF' for (u, v, d) in G.edges(data=True)])
    network_nodes.set_edgecolor('#FFFFFF')
    nx.draw_networkx_edges(G, pos, alpha=0.4, node_size=10, width = 1, edge_color='#808080')
    nx.draw_networkx_labels(G, pos, font_size=12, font_weight = 'bold', labels = dict([(u, '%s\n(%s)' % (d['label'], pretty_print(d['size']))) for u, v, d in G.edges(data=True)]))
    
    # adjust the plot limits
    xmax = 1.02 * max(x for x, y in pos.values())
    ymax = 1.02 * max(y for x, y in pos.values())
    plt.xlim(0, xmax)
    plt.ylim(0, ymax)
    plt.xticks([])
    plt.yticks([])
    plt.axis('equal')

    plt.subplots_adjust(hspace = 0, wspace = 0, right = 0.995, left = 0.005, top = 0.995, bottom = 0.005)

    plt.text(0.5, 0.97, "Topology",
             horizontalalignment='center',
             transform=plt.gca().transAxes)

    #plt.axis('off')
    plt.show()


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Graph Representation of Minimum Entropy Decomposition Topology')
    parser.add_argument('topology_file', metavar = 'TOPOLOGY',\
                        help = 'Description of the topology in text format')
    parser.add_argument('--output-file', default = None, metavar = 'OUTPUT_FILE',\
                        help = 'File name for the figure to be stored. File name\
                                must end with "png", "jpg", or "tiff".')
    parser.add_argument('--title', default = None, metavar = 'TITLE',\
                        help = 'Title for the figure (project name would be appropriate).')


    args = parser.parse_args()

    topology(args.topology_file, args.output_file, args.title)

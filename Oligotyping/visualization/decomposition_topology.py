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

import cPickle
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

from Oligotyping.utils.utils import pretty_print

final_nodes = []
parent_nodes = []

def topology_graph(topology_dict_path, match_levels = False):
    G = nx.MultiDiGraph()
    
    topology = cPickle.load(open(topology_dict_path))
   
    nodes = {}
    levels = []
    
    for node_id in topology:
        node = topology[node_id]
        
        if node.killed:
            continue

        if not node.children:
            final_nodes.append(node_id)
        else:
            parent_nodes.append(node_id)

        nodes[node_id] = {'size': node.size, 'parent': node.parent, 'level': node.level,
                          'children': [child_node_id for child_node_id in node.children if topology.has_key(child_node_id) and not topology[child_node_id].killed], 'type': 'node'}
        if node.freq_curve_img_path:
            nodes[node_id]['freq_curve_img_path'] = node.freq_curve_img_path
        levels.append(int(node.level))
    
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
                                                               'level': level, 'children': [], 'type': 'spam', 'freq_curve_img_path': None}
                    else:
                        new_nodes[node_id + ':l%d' % level] = {'size': node['size'], 'parent': node_id + ':l%d' % (level - 1),
                                                               'level': level, 'children': [], 'type': 'spam', 'freq_curve_img_path': None}

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
                G.add_edge(node_id, node['parent'], size = int(node['size']), label = label,\
                           image = node['freq_curve_img_path'] if node.has_key('freq_curve_img_path') else None,\
                           final_node = True if not node['children'] else False)
            else:
                G.add_edge(node_id, node['parent'], size = int(node['size']), label = '',\
                           image = node['freq_curve_img_path'] if node.has_key('freq_curve_img_path') else None,\
                           final_node = True if not node['children'] else False)
   
    for node_id in nodes['root']['children']:
        node = nodes['root']
        G.add_edge('root', node_id, size = int(nodes['root']['size']), label = 'root',\
                           image = node['freq_curve_img_path'] if node.has_key('freq_curve_img_path') else None,\
                           final_node = False)
   
    return (G, nodes)


def topology(topology_dict_path, output_file = None, title = None):
    G, nodes_dict = topology_graph(topology_dict_path)

    number_of_edges = G.number_of_edges()
    number_of_nodes = G.number_of_nodes()

    print("Loaded %d edges and %d nodes." % (number_of_edges, number_of_nodes))

    plt.figure(figsize=(24, 16))
    
    # use graphviz to find radial layout
    # twopi, gvcolor, wc, ccomps, tred, sccmap, fdp, circo, neato, acyclic, nop, gvpr, dot
    pos=nx.graphviz_layout(G, prog="fdp")

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

    parent_nodes_network = nx.draw_networkx_nodes(G, pos, nodelist = parent_nodes, node_shape = 'o', node_size = [sizes[i] for i in parent_nodes], node_color = '#EFEFEF')
    final_nodes_network = nx.draw_networkx_nodes(G, pos, nodelist = final_nodes, node_shape = 'o', node_size = [sizes[i] for i in final_nodes], node_color = '#FAFFFA')
    parent_nodes_network.set_edgecolor('#888888')
    final_nodes_network.set_edgecolor('#88BB00')
    nx.draw_networkx_edges(G, pos, alpha=0.4, node_size=10, width = 1, edge_color='#808080')
    nx.draw_networkx_labels(G, pos, font_size=8, font_weight = 'bold', labels = dict([(u, '%s\n(%s)' % (d['label'], pretty_print(d['size']))) for u, v, d in G.edges(data=True)]))
    
    # adjust the plot limits
    xmax = 1.02 * max(x for x, y in pos.values())
    ymax = 1.02 * max(y for x, y in pos.values())
    plt.xlim(0, xmax)
    plt.ylim(0, ymax)
    plt.xticks([])
    plt.yticks([])

    plt.subplots_adjust(hspace = 0, wspace = 0, right = 0.995, left = 0.005, top = 0.995, bottom = 0.005)

    plt.text(0.03, 0.97, title or "Topology", fontsize='xx-large',
             fontname="Arial", fontweight="bold", transform=plt.gca().transAxes)

    ax=plt.gca()
    plt.setp(ax, frame_on=False)
    #plt.axis('off')

    if nodes_dict['root'].has_key('freq_curve_img_path'):
        AX=plt.gca()
        f=plt.gcf()

        for node in nodes_dict.keys():
            (x, y) = pos[node]
            xt,yt = AX.transData.transform((x, y)) # figure coordinates
            xf, yf = f.transFigure.inverted().transform((xt, yt)) # axes coordinates
            print xf, yf
            if node == 'root':
                imsize = 0.04
            else:
                imsize = 0.025
            img =  mpimg.imread(nodes_dict[node]['freq_curve_img_path'])
            a = plt.axes([xf - imsize / 2.0, yf - imsize / 2.0, imsize, imsize ])
            a.imshow(img)
            a.axis('off')

    if output_file:
        plt.savefig(output_file)
    else:
        plt.show()


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Graph Representation of Minimum Entropy Decomposition Topology')
    parser.add_argument('topology_dict', metavar = 'TOPOLOGY',\
                        help = 'Serialized topology dictionary')
    parser.add_argument('--output-file', default = None, metavar = 'OUTPUT_FILE',\
                        help = 'File name for the figure to be stored. File name\
                                must end with "png", "jpg", or "tiff".')
    parser.add_argument('--title', default = None, metavar = 'TITLE',\
                        help = 'Title for the figure (project name would be appropriate).')

    args = parser.parse_args()

    topology(args.topology_dict, args.output_file, args.title)

# -*- coding: utf-8 -*-

# Copyright (C) 2010 - 2012, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

import copy
import random
import matplotlib.cm as cm

#
# all available colormaps in matplotlib can be seen via
# http://wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps
#

def getColor(name, n):
    return cm.get_cmap(name, lut=n+2)

def get_hex_color(rgba_color):
    hex_color = '#'
    for t in rgba_color[0:3]:
        h = str(hex(int(t * 255)))[2:]
        hex_color += '00' if h == '0' else h
    return hex_color + '0' * (7 - len(hex_color))

def random_colors(oligotypes, output_file_path = None, colormap = 'Paired'):

    oligotypes_shuffled = copy.deepcopy(oligotypes)
    random.shuffle(oligotypes_shuffled)

    color_dict = {}

    colors = getColor(colormap, len(oligotypes_shuffled))
    for i in range(0, len(oligotypes_shuffled)):
        color_dict[oligotypes_shuffled[i]] = get_hex_color(colors(i))

    if output_file_path:
        output_file = open(output_file_path, 'w')
        for oligotype in oligotypes:
            output_file.write('\n'.join(['%s\t%s\n' % (oligotype, color_dict[oligotype])]))
    
    return color_dict

def get_list_of_colors(number_of_colors, colormap = 'OrRd'):
    color_list = getColor(colormap, number_of_colors)
    return [get_hex_color(color_list(i)) for i in range(0, number_of_colors)] 

def get_color_shade_dict_for_list_of_values(values, colormap = 'OrRd'):
    color_shade_dict = {}

    colors = getColor(colormap, len(values) * 1000)

    max_val = max(values) if max(values) > 1 else 1

    for val in values:
        i = float(val) / max_val
        color_shade_dict[val] = get_hex_color(colors(i))

    return color_shade_dict

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Generate Random Colors')
    parser.add_argument('oligotypes', metavar = 'FASTA_FILE', help = 'Oligotype sequences in FASTA format')
    parser.add_argument('--output-file', default = None, metavar = 'OUTPUT_FILE',\
                        help = 'File name to store random colors')
    parser.add_argument('--colormap', default = 'Accent', metavar = 'MATPLOTLOB_COLORMAP',\
                        help = 'see http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps')


    args = parser.parse_args()

    oligotypes = [line.strip() for line in open(parser.parse_args().oligotypes) if not line.startswith('>')]
    
    colors_dict = random_colors(oligotypes, output_file_path = args.output_file, colormap = args.colormap) 

    if not args.output_file:
        for oligo in colors_dict:
            print('%s: %s' % (oligo, colors_dict[oligo]))

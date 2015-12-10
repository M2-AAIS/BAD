#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
import argparse
from os import path

parser = argparse.ArgumentParser(description='Plot data from output of the black hole simulation.')
parser.add_argument('--index', type=int, default=1,
                    help='Index to plot (default: %(default)s)')
parser.add_argument('--map-dir', type=str, default='maps', dest='map_dir',
                    help='Directory containing paths (default: %(default)s)')
args = parser.parse_args()


if __name__ == '__main__':
    filename = path.join(args.map_dir, 'map_{:0>5}.dat'.format(args.index))
    with open(filename, 'r') as f:
        f.readline()
        T = [float(i) for i in f.readline().split()]
        f.readline()
        Sigma = [float(i) for i in f.readline().split()]
        lines = ''.join(f.readlines())

    # Read T, Sigma, Q and tau (one liner yeah!)
    Q, tau = [pd.DataFrame([fline.split() for fline in l.split('\n')[1:-1]],
                                     dtype=float)
                                     for l in lines.split('#')][1:]

    extents = Sigma[0], Sigma[1], T[0], T[1]

    #fig, ax = plt.subplots()

    for element in ({'data': Q, 'title': 'Q^+ - Q^-', 'threshold': 5e13},
                    {'data': tau, 'title': '\\tau', 'threshold': 1e-5}):
        plt.figure()
        plt.title(u"${}$".format(element['title']))
        plt.imshow(element['data'], interpolation='none', cmap='viridis', extent=extents, aspect='auto',
                   norm=matplotlib.colors.SymLogNorm(element['threshold']))
        plt.colorbar()
        plt.xlabel('$\Sigma$')
        plt.ylabel('$T$')

    plt.show()

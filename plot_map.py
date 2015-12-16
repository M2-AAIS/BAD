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
        nb_T = int(f.readline().split()[0])
        f.readline()
        nb_S = int(f.readline().split()[0])
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

    X = np.logspace(np.log10(Sigma[0]), np.log10(Sigma[1]), nb_S)
    Y = np.logspace(np.log10(T[1]), np.log10(T[0]), nb_T)

    #fig, ax = plt.subplots()

    for element in ({'var': 'Q', 'data': Q, 'title': u"Chauffage spécifique $Q^+ - Q^-$ (erg·g$^{-1}$·s$^{-1}$)", 'threshold': 1e15, 'ticks': [-1e20,0,1e20]},
                    {'var': 'tau', 'data': tau, 'title': u"Opacité $\\tau$", 'threshold': 1e-5, 'ticks': [0.06, 1, 10]}):
        plt.figure()
        plt.title(element['title'])
        plt.xscale('log')
        plt.yscale('log')
        CSQ = plt.contour(X, Y, Q, origin='image', levels=[0], colors=('w'), extent=extents)
        plt.clabel(CSQ, inline=1, fontsize=10)
        CStau = plt.contour(X, Y, tau, origin='image', levels=[0.06,1, 10], colors=('r','m','k'), extent=extents)
        plt.clabel(CStau, inline=1, fontsize=10)
        plt.imshow(element['data'], interpolation='none', cmap='viridis', extent=extents, aspect='equal',
                   norm=matplotlib.colors.SymLogNorm(element['threshold']))
        plt.colorbar(ticks=element['ticks'],shrink=1)
        plt.xlabel(u"Densité surfacique $\Sigma$ (g·cm$^{-2}$)")
        plt.ylabel(u"Température $T$ (K)")
        plt.tight_layout()
        plt.savefig('maps/'+element['var']+'_map.pdf', transparent=True, dpi=300, bbox_inches='tight', pad_inches=0)

    plt.show()

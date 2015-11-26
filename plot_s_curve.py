#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from numpy import array, log
import sys
import os

import matplotlib.animation as animation

fig = plt.figure()

inpath = sys.argv[1]

if os.path.isfile(inpath):
    print('Visiting {}'.format(inpath))
    filenames = [inpath]
else:
    _filenames = os.listdir(inpath)
    _filenames.sort()
    filenames = [inpath + '/' + fname for fname in _filenames if '_tot.dat' in fname]
    
    print('Visiting all files of {}'.format(inpath))

axline, = plt.plot(0, 0, 'o')

def draw_once(filename):
    x = []
    y = []
    if not 'tot.dat' in filename:
        return ([0], [0])
    else:
        print('Visiting {}'.format(filename))
        outfile = filename.replace('.dat', '.png')
        
    for line in open(filename):
        data = line.replace('\n', '').split()
        try :
            print (data)
            xData = float(data[0])
            yData = float(data[1])
            x.append(xData)
            y.append(yData)
        except ValueError:
            pass

    axline.set_xdata(x)
    axline.set_ydata(y)

    return axline,

def init():
    print('Initialisation')
    plt.ylabel('$\log T$')
    plt.xlabel('$\log \Sigma$')
    plt.xlim(1.8, 3.2)
    plt.ylim(6, 8)
    plt.grid()

if len(filenames) > 1:
    ani = animation.FuncAnimation(fig, draw_once, filenames, init_func=init, interval=10)
    plt.show()
else:
    init()
    draw_once(filenames[0])
    plt.show()
# x, y = draw_once(filenames[2])
# plt.plot(x, y, 'o')


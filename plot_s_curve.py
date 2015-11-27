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

axline, = plt.plot(0, 0, 'o', label='')

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
            xData = float(data[0])
            yData = float(data[1])
            x.append(xData)
            y.append(yData)
        except ValueError:
            pass

    # this is a hack to get the x
    arr = filename.split('/')[-1].split('_')
    number = int(arr[2])

    axline.set_xdata(x)
    axline.set_ydata(y)
    axline.set_label('$(T, \Sigma)_{'+str(number)+'}$')
    plt.legend(loc='upper left')
    return axline,

def init():
    print('Initialisation')
    plt.ylabel('$\log T$')
    plt.xlabel('$\log \Sigma$')
    plt.xlim(1.8, 3.2)
    plt.ylim(6, 8)
    plt.grid()
    plt.legend()

if len(filenames) > 1:
    ani = animation.FuncAnimation(fig, draw_once, filenames, init_func=init, interval=10)
    ani.save('s_curve.mp4', writer='ffmpeg', fps=10, bitrate=10000, dpi=180)
    plt.show()
else:
    init()
    draw_once(filenames[0])
    plt.show()
# x, y = draw_once(filenames[2])
# plt.plot(x, y, 'o')


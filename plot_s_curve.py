#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from numpy import array, log
import sys
import os
import argparse

import matplotlib.animation as animation

fig = plt.figure()

parser = argparse.ArgumentParser(description='Plot data from output of the black hole simulation.')
parser.add_argument('--path', default='s_curves/',
                    help='path to the S curves directory')
parser.add_argument('--index', default=-1, type=int,
                    help='Give a particular index to plot')
parser.add_argument('--video', action='store_true',
                    help='save a video (MUCH slower)')
args = parser.parse_args()

fullpath = os.path.join(args.path,
                        'Temperature_Sigma_{:0>5}_tot.dat'.format(args.index))
print(fullpath)
if os.path.isfile(fullpath):
    print('Visiting {}'.format(fullpath))
    filenames = [fullpath]
else:
    _filenames = os.listdir(args.path)
    _filenames.sort()
    filenames = [os.path.join(args.path,fname) for fname in _filenames if '_tot.dat' in fname]

    print('Visiting all files of {}'.format(args.path))

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

    axline.set_xdata(log(x)/log(10))
    axline.set_ydata(log(y)/log(10))
    axline.set_label('$(T, \Sigma)_{'+str(number)+'}$')
    plt.legend(loc='upper left')
    return axline,

def init():
    print('Initialisation')
    plt.ylabel('$\log T$')
    plt.xlabel('$\log \Sigma$')
    plt.xlim(1.4, 3.6)
    plt.ylim(5.2, 7.5)
    plt.grid()
    plt.legend()

if len(filenames) > 1:
    ani = animation.FuncAnimation(fig, draw_once, filenames, init_func=init, interval=10)
    if args.video:
        ani.save('s_curve.mp4', writer='ffmpeg', fps=10, bitrate=10000, dpi=180)
    plt.show()
else:
    init()
    draw_once(filenames[0])
    plt.show()


#!/usr/bin/python
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.animation as animation
import numpy as np
from itertools import tee

parser = argparse.ArgumentParser(description='Plot data from output of the black hole simulation.')
parser.add_argument('--no-video', action='store_true',
                    help='do not save a video (faster)')
parser.add_argument('--video-file', metavar='file', default='evolution.mp4',
                    help='save the video as (default: %(default)s).')
parser.add_argument('--plot', metavar='dump_id', nargs='*', type=int, default=[],
                    help='dump_ids to plot.')
parser.add_argument('--c-points', metavar='file',
                    help='Critical point file (default: %(default)s).',
                    default='critical_points/file.dat')
parser.add_argument('--i-conditions', metavar='file',
                    help='Initial conditions file (default: %(default)s).', default='CI.dat')
parser.add_argument('--s-curves', metavar='n', nargs='+',
                    help='S curves to plot (default: %(default)s).', default=[1, 10, 100, 200])
parser.add_argument('--s-curves-dir', metavar='dir', nargs=1,
                    help='S curves directory (default: %(default)s).', default='s_curves')
args = parser.parse_args()

fig, ((ax11, ax12), (ax21, ax22)) = plt.subplots(2, 2)

class Data:
    ''' Class that opens filename and reads it. You can directly iterate over it.'''
    def __init__(self, filename):
        self.filename = filename
        self.data = {}
        self.times = {}

        i = 0
        with open(self.filename, 'r') as f:
            for chunk in self.chunkYielder(f):
                dumpNumber, time, data = self.processChunk(chunk)
                self.data[i] = data
                self.times[i] = time
                i += 1

    def chunkYielder(self, f):
        '''Read the file, split it arround '#' and yields it piece by piece. '''
        # arf, not beautiful, should read it more carefully
        lines = [l.replace('\n', '') for l in f]

        i = 0
        
        while i < len(lines):
            niter = int(lines[i].split('#')[1])
            time = float(lines[i+1].split('#')[1])
            headers = lines[i+2].split()
            data = []
            i += 3

            while i < len(lines) and '#' not in lines[i]:
                data.append([float(val) for val in lines[i].split()])
                i += 1
            yield niter, time, headers, data
                    
                    
                
                
    def processChunk(self, chunk):
        ''' Read the lines and get the following data in the lines:
        1: dump number
        2: time
        3: headers
        4+: data'''

        dumpNumber = chunk[0]
        time = chunk[1]
        headers = chunk[2]
        data = pd.DataFrame(data=chunk[3], columns=headers, dtype=float)
        
        return dumpNumber, time, data

    def __iter__(self):
        for key in self.data:
            yield (key, self.data[key])


lines = {'r-T': 0,
         'r-Sigma': 0,
         'Sigma-T': []}

class colorLooper:
    def __init__(self, colors= ['red', 'green', 'blue', 'yellow']):
        self.colors = colors
        self.i = 0

    def __next__(self):
        self.i += 1
        return self.colors[self.i - 1]
    
    def __iter__(self):
        self.i = 0
        while self.i < len(self.colors):
            yield self.__next__()

    def reset(self):
        self.i = 0
        
def init(ic, crit_pts, s_curves, initial_data):
    ''' Plot the initial conditions, the critical points and the s_curves'''
    ax11.plot(ic['r'], ic['T'], '--', label='Initial conditions')
    ax11.plot(ic['r'], crit_pts['Temp_thin'], '--', label='critical')
    
    ax11.set_xlabel('$r\ (cm)$')
    ax11.set_ylabel('$T\ (K)$')

    lines['r-T'] = ax11.plot(initial_data['r'], initial_data['T'], label='iteration 0')[0]
    
    ax11.set_yscale('log')
    ax11.grid()
    ax11.legend()

    

    ax12.set_xlabel('$r\ (cm)$')
    ax12.set_ylabel('$\Sigma\ (g.cm^{-2})$')
    ax12.plot(ic['r'], ic['Sigma'], '--')
    ax12.plot(ic['r'], crit_pts['Sigma_thin'], '--')

    lines['r-Sigma'] = ax12.plot(initial_data['r'], initial_data['Sigma'])[0]
    
    ax12.grid()
    ax12.set_yscale('log')

    
    lines['r-Mdot'] = ax21.plot(initial_data['r'], initial_data['M_dot'])[0]
    
    ax21.set_xlabel('$r\ (cm)$')
    ax21.set_ylabel('$\dot{M}$')

    colorsIter = colorLooper()
    for ind, s_curve in s_curves:
        ax22.plot(s_curve['Surface_density'], s_curve['Temperature'],
                  label='$x_{'+str(ind)+'}$',
                  c=colorsIter.__next__())

    # add the initial data on the s_curve plot
    colorsIter.reset()
    lines['Sigma-T'] = [ (ind, ax22.plot(initial_data['S'][ind],
                                         initial_data['T'][ind], 'o',
                                         c=colorsIter.__next__())[0])
                     for ind, foo in s_curves]
                             

    ax22.set_xlabel('$\log\Sigma\ (g.cm^{-2})$')
    ax22.set_ylabel('$\log T\ (K)$')
    ax22.legend()
    ax22.grid()

def plotData(args):
    ''' Update with the new data. Parameters:
    args: array that contains
      index: the indexes of the datadump
      data : the data
    '''

    index, data = args
    print('Plotting {}'.format(index))
    lines['r-T'].set_ydata(data['T'])
    lines['r-T'].set_label('iteration {}'.format(index))

    lines['r-Mdot'].set_ydata(data['M_dot'])

    lines['r-Sigma'].set_ydata(data['Sigma'])

    for ind, line in lines['Sigma-T']:
        x = data['S'][ind]
        y = data['T'][ind]
        line.set_xdata(x)
        line.set_ydata(y)

    ax11.legend()
    

if __name__ == '__main__':
    print('Reading initial conditions')
    ic = pd.read_csv('CI.dat', delim_whitespace=True)
    
    print('Reading output file')
    simulData = Data('output.dat')

    print('Reading critical points')
    crit_pts = pd.read_csv(args.c_points, delim_whitespace=True)

    print('Reading S_curves')
    s_curves_indexes = [1, 10, 100, 250]
    s_curves = [ (ind, pd.read_csv(args.s_curves_dir + '/Temperature_Sigma_{:0>5}_tot.dat'.format(ind),
                                   delim_whitespace=True, dtype=float, header=0))
                 for ind in args.s_curves]
    

    initFun = lambda: init(ic, crit_pts, s_curves, simulData.data[0])

    if len(simulData.data) == 1 or len(args.plot) == 1:
        initFun()
        plotData((0, simulData.data[args.plot[0]]))
        plt.show()
    else:
        if (len(args.plot) > 0):
            nFrames = len(args.plot)
            data = [(i, d) for i, d in simulData if i in args.plot]
        else:
            nFrames = len(simulData.data)
            data = simulData
            
        print('Go take a coffee, still {} frames to go.'.format(nFrames))

        ani = animation.FuncAnimation(fig, plotData, data, init_func=initFun, interval=10)
        if not args.no_video:
            ani.save(args.video_file, writer='ffmpeg', fps=10, bitrate=10000, dpi=180)

        plt.show()

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
parser.add_argument('--dont-loop', action='store_false', dest='loop', default=True,
                    help='prevent from looping infinitely')
parser.add_argument('--video-file', metavar='file', default='evolution.mp4',
                    help='save the video as (default: %(default)s).')
parser.add_argument('--plot', metavar='dump_id', nargs='*', type=int, default=[],
                    help='dump_ids to plot.')
parser.add_argument('--c-points', metavar='file',
                    help='Critical point file (default: %(default)s).',
                    default='critical_points/file.dat')
parser.add_argument('--i-conditions', metavar='file',
                    help='Initial conditions file (default: %(default)s).', default='CI.dat')
parser.add_argument('--s-curves', metavar='n', nargs='+', type=int,
                    help='S curves to plot (default: %(default)s).', default=[1, 10, 20, 255])
parser.add_argument('--s-curves-dir', metavar='dir', nargs=1,
                    help='S curves directory (default: %(default)s).', default='s_curves')
parser.add_argument('--reread', action='store_true', default=False,
                    help='Reread file when fully read.')
args = parser.parse_args()

fig, ((ax11, ax12), (ax21, ax22)) = plt.subplots(2, 2)

class Data:
    ''' Class that opens filename and reads it. You can directly iterate over it.'''
    def __init__(self, filename):
        self.filename = filename
        self.data = {}
        self.times = {}
        self.restrictTo = [] # list of indexes to restrict to

        self.file = open(self.filename, 'r')
        # read ahead one line
        self.lastLineRead = self.file.readline()

    def getChunk(self):
        '''Get a chunk (niter, time, headers + data)'''
        f = self.file

        if self.lastLineRead == '':
            return -1, -1, [], [], True

        # first line to read == last line read at previous call            
        niter = int(self.lastLineRead.split('#')[1])
        time = float(f.readline().split('#')[1])
        headers = f.readline().split()
        eof = False # flag to set if end of file

        data = []
        line = f.readline()
        # read until reaching a commented line
        while '#' not in line:
            if line == '':
                eof = True
                break
            # remove last character ('\n')
            data.append(line[:-1].split())
            
            line = f.readline()

        pdData = pd.DataFrame(data=data, columns=headers, dtype=float)
        # store the last line read
        self.lastLineRead = line
        return niter, time, headers, pdData, eof

    def __iter__(self):
        ''' Infinite iterator over the whole file'''
        self.restrictTo.sort()
        
        if len(self.restrictTo) > 0:
            noRestriction = False
        else:
            noRestriction = True
            
        i = 0
        niter, time, headers, data, eof = self.getChunk()
        while not eof:
            self.data[niter] = data
            self.times[niter] = time
            if noRestriction or niter in self.restrictTo:                
                yield (niter, time, data)

            # if we have a restriction list and the current iteration
            # is after the last accepted item, stop reading the file
            if not noRestriction and niter > self.restrictTo[-1]:
                break
            niter, time, headers, data, eof = self.getChunk()
            if eof and self.reread:
                self.file.close()
                self.file = open(self.filename, 'r')
                self.lastLineRead = self.file.readline()
                eof = False
                

        # if the loop flag is activated, loop over saved data
        if self.loop:
            if noRestriction:
                dataKeys = list(self.data.keys())
            else:
                dataKeys = [key for key in self.data.keys() if key in self.restrictTo]
                
            dataKeys.sort()
            while True:
                for key in dataKeys:
                    yield (key, self.times[key], self.data[key])


    def get(self):
        ''' Get a chunk of the file without moving the position of the reader in the file '''
        pos = self.file.tell()
        niter, time, headers, data, eof = self.getChunk()
        self.file.seek(pos)
        return niter, time, headers, data
        
    def __del__(self):
        ''' Close the file reader'''
        self.file.close()


lines = {'r-T': 0,
         'r-Sigma': 0,
         'Sigma-T': []}

class colorLooper:
    def __init__(self, colors= ['red', 'green', 'blue', 'grey']):
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
    ## Top left panel
    ax11.plot(ic['r'], ic['T'], '--', label='Initial conditions')
    ax11.plot(ic['r'], crit_pts['Temp_thick'], '--', label='critical')
    
    ax11.set_xlabel('$r\ (cm)$')
    ax11.set_ylabel('$T\ (K)$')

    lines['r-T'] = ax11.plot(initial_data['r'], initial_data['T'])[0]
    
    ax11.set_yscale('log')
    ax11.grid()
    # ax11.legend()

    # add ticks corresponding to the s_curve
    ticks      = list(ax11.get_xticks())
    ticklabels = list(ax11.get_xticklabels())
    for ind, s_curve in s_curves:
        ticks.append(ic['r'][ind])
        ticklabels.append('$r_{'+str(ind)+'}$')

    ax11.set_xticks(ticks)
    ax11.set_xticklabels(ticklabels)


    ## Top right panel
    ax12.set_xlabel('$r\ (\mathrm{cm})$')
    ax12.set_ylabel('$\Sigma\ (\mathrm{g.cm^{-2}})$')
    ax12.plot(ic['r'], ic['Sigma'], '--')
    ax12.plot(ic['r'], crit_pts['Sigma_thick'], '--')

    lines['r-Sigma'] = ax12.plot(initial_data['r'], initial_data['Sigma'])[0]
    
    ax12.grid()
    ax12.set_yscale('log')
    
    lines['r-Mdot'] = ax21.plot(initial_data['r'], initial_data['M_dot'])[0]
    
    # add ticks corresponding to the s_curve
    ticks      = list(ax12.get_xticks())
    ticklabels = list(ax12.get_xticklabels())
    for ind, s_curve in s_curves:
        ticks.append(ic['r'][ind])
        ticklabels.append('$r_{'+str(ind)+'}$')

    ax12.set_xticks(ticks)
    ax12.set_xticklabels(ticklabels)


    ## Bottom left panel
    ax21.set_xlabel('$r\ (\mathrm{cm})$')
    ax21.set_ylabel('$\dot{M}\ (\mathrm{g.s^{-1}})$')
    ax21.grid()
 
    ## Bottom right panel
    colorsIter = colorLooper()
    for ind, s_curve in s_curves:
        ax22.plot(s_curve['Surface_density'], s_curve['Temperature'],
                  label='$r_{'+str(ind)+'}$',
                  c=colorsIter.__next__())

    # add the initial data on the s_curve plot
    colorsIter.reset()
    lines['Sigma-T'] = [ (ind, ax22.plot(initial_data['Sigma'][ind],
                                         initial_data['T'][ind], 'o',
                                         c=colorsIter.__next__())[0])
                     for ind, foo in s_curves]
                             

    ax22.set_xlabel('$\Sigma\ (\mathrm{g.cm^{-2}})$')
    ax22.set_ylabel('$T\ (\mathrm{K})$')
    ax22.legend(loc='upper right')
    ax22.grid()

def plotData(args):
    ''' Update with the new data. Parameters:
    args: array that contains
      index: the indexes of the datadump
      data : the data
    '''

    index, time, data = args
    print('Plotting {}'.format(index))
    lines['r-T'].set_ydata(data['T'])
    # lines['r-T'].set_label('$t = {:.2}s$'.format(time))

    lines['r-Mdot'].set_ydata(data['M_dot'])

    lines['r-Sigma'].set_ydata(data['Sigma'])

    for ind, line in lines['Sigma-T']:
        x = data['Sigma'][ind]
        y = data['T'][ind]
        line.set_xdata(x)
        line.set_ydata(y)

    # ax11.legend()

    fig.suptitle('$t = {:.2f}s$, iteration {}'.format(time, index))
    
if __name__ == '__main__':
    print('Reading initial conditions')
    ic = pd.read_csv('CI.dat', delim_whitespace=True)
    
    print('Reading output file')
    simulData = Data('output.dat')

    print('Reading critical points')
    crit_pts = pd.read_csv(args.c_points, delim_whitespace=True)

    print('Reading S_curves')
    s_curves = [ (ind,
                  pd.read_csv(args.s_curves_dir + '/Temperature_Sigma_{:0>5}_tot.dat'.format(ind),
                              delim_whitespace=True, dtype=float, header=0))
                 for ind in args.s_curves]
    

    data0 = simulData.get()[-1]
    initFun = lambda: init(ic, crit_pts, s_curves, data0)

    if len(args.plot) == 1:
        initFun()
        plotData((0, data0, 0))
        plt.show()
    else:
        simulData.loop = args.loop
        simulData.reRead = args.reread
            
        if (len(args.plot) > 0):
            simulData.restrictTo = args.plot

        data = simulData

        ani = animation.FuncAnimation(fig, plotData, data, init_func=initFun, interval=1)
        if not args.no_video:
            # temporaly deactivate looping for the video
            simulData.loop = False
            ani.save(args.video_file, writer='ffmpeg', fps=10, bitrate=10000, dpi=180)
            simulData.loop = args.loop

        # clear the axes
        for ax in (ax11, ax12, ax21, ax22):
            ax.axes.cla()
        
        plt.show()

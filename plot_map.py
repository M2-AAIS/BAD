import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import argparse
from os import path

parser = argparse.ArgumentParser(description='Plot data from output of the black hole simulation.')
parser.add_argument('--index', type=int, default=1,
                    help='Index to plot (default: %{default}s)')
parser.add_argument('--map-dir', type=str, default='maps', dest='map_dir',
                    help='Directory containing paths (default: %{default}s)')
args = parser.parse_args()


if __name__ == '__main__':
    filename = path.join(args.map_dir, 'map_{:0>5}.dat'.format(args.index))
    with open(filename, 'r') as f:
        lines = ''.join(f.readlines())

    # Read T, Sigma, Q and tau (one liner yeah!)
    T, Sigma, Q, tau = [pd.DataFrame([fline.split() for fline in l.split('\n')[1:-1]],
                                     dtype=float)
                                     for l in lines.split('#')][1:]

    endT = len(T[0]) - 1
    endSigma = len(Sigma[:][0]) - 1
    extents = T[0][0], T[0][endT], Sigma[0][0], Sigma[endSigma][0]

    fig, ax = plt.subplots()
    # FIXME: be able to plot for negative values
    plt.imshow(Q, interpolation='none', cmap='viridis', extent=extents, aspect='auto',
               norm=matplotlib.colors.LogNorm())
    plt.colorbar()
    plt.xlabel('Temperature')
    plt.ylabel('Sigma')

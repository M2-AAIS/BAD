#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import sys

x = []

y = []

infile = open(sys.argv[1])

for line in infile:
    data = line.replace('\n','').split()
    print(data)
    x.append(float(data[0]))
    y.append(float(data[1]))

figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()
plt.xscale('log')
plt.yscale('log')
plt.plot(x,y)
plt.show()

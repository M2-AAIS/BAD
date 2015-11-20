#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from numpy import array, log
import sys

x = []

y = []

infile = open(sys.argv[1])

for line in infile:
    data = line.replace('\n','').split()
    print(data)
    x.append(float(data[0]))
    y.append(float(data[1]))

#x = array(x)
#y = array(y)

figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()
#plt.plot(log(x),log(y))
plt.plot(x,y)
plt.show()

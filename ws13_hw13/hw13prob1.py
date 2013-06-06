#!/usr/bin/env python

import numpy as np
import pylab
import matplotlib.pyplot as plt
from math import pi

np.random.seed(1)

nmax = 100000

npi = np.zeros(nmax)
n = np.zeros(nmax)
ncount = 0


for i in range(0,nmax):
	randx = 2.0*np.random.rand()-1.0
	randy = 2.0*np.random.rand()-1.0
	rad2 = randx**2 + randy**2
	if (rad2 <= 1.0):
		ncount += 1
	npi[i] = (4.0*ncount)/(i+1)
	n[i] = i+1
print npi[nmax-1]
plt.loglog(n, abs(npi - pi), linestyle='-', color = 'b', linewidth=2, label='MC')
plt.loglog(n, 1.0/np.sqrt(n), linestyle = '-', color='m', linewidth=2, label=r'$N^{-1/2}$')
plt.ylim([1.0e-04,1.0])
plt.xlabel('Number of experiments', size=18)
plt.ylabel('|$\pi_{MC} - \pi$|', size=18)
plt.legend()
plt.show()



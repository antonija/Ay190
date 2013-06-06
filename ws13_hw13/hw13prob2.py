#!/usr/bin/env python

import numpy as np
import pylab
import matplotlib.pyplot as plt

def ranwalk():
	loc = 0.0 # initial position
	numstep = 0 # step counter
	while (abs(loc) < 1.0):

		if (np.random.rand() < 0.5):
			step = -lam # step to the left
		else:
			step = lam # step to the right
		loc += step
		numstep += 1
	return numstep

####################################################

np.random.seed(1)

ntot = 10 
num = np.zeros(ntot)
lam = 0.01
for i in range(0, ntot):
    num[i] = 1.0*ranwalk()
print np.mean(num)

nexp = [10, 50, 100, 500, 1000, 5000, 10000, 50000]
nsteps = [7189.8, 7579.08, 8837.76, 9478.592, 9519.028, 10094.1628, 9962.898, 10024.332]
plt.semilogx(nexp, nsteps, 'go', ms=10)
plt.axhline(y=10000.0, color='0.4', linestyle=':', linewidth=1.0)
plt.xlabel('Number of experiments', size=18)
plt.ylabel('$\overline{N}$', size=18)
plt.show()


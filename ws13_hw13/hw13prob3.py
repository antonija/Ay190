#!/usr/bin/env python

import numpy as np
import pylab
import matplotlib.pyplot as plt
from math import sqrt, sin, cos, pi

def ranwalk(lam):
	locx = [0.0] # initial position
	locy = [0.0] # initial position
	dist = 0.0
	numstep = 0 # step counter
	step = np.zeros(2)
	stepcount = 0
	np.random.seed()
	while (dist < 1.0):		
		rand = 8.0*np.random.rand()
		if (rand < 1.0):
			step[0] = lam
			step[1] = 0.0
		elif (rand < 2.0):
			step[0] = lam/1.41421356237
			step[1] = lam/1.41421356237
		elif (rand < 3.0):
			step[0] = 0.0
			step[1] = lam
		elif (rand < 4.0):
			step[0] = -lam/1.41421356237
			step[1] = lam/1.41421356237
		elif (rand < 5.0):	
			step[0] = -lam
			step[1] = 0.0
		elif (rand < 6.0):
			step[0] = -lam/1.41421356237
			step[1] = -lam/1.41421356237
		elif (rand < 7.0):
			step[0] = 0.0
			step[1] = -lam
		else:
			step[0] = lam/1.41421356237
			step[1] = -lam/1.41421356237
		
		stepcount += 1
		locx.append(locx[stepcount-1] + step[0])
		locy.append(locy[stepcount-1] + step[1])
		dist = sqrt(locx[stepcount]**2 + locy[stepcount]**2)
		#print dist
	print stepcount
	plt.plot(locx,locy)

####################################################

l = 0.015
np.random.seed(1)
ranwalk(l)
np.random.seed(2)
ranwalk(l)
np.random.seed(3)
ranwalk(l)
np.random.seed(4)
ranwalk(l)

circx = np.zeros(360) 
circy = np.zeros(360)    
angleStep = pi*2.0/360.0
for a in range(0, 360):
	circx[a] = sin(a*angleStep)
	circy[a] = cos(a*angleStep)
plt.plot(circx, circy, linewidth=3, color='black')
plt.xlim([-1,1])
plt.ylim([-1,1])
plt.axes().set_aspect('equal')
plt.show()

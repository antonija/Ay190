#!/usr/bin/python

import numpy as np
import pylab 
import matplotlib.pyplot as plt
from math import pi,pow, sqrt, sin


# MAIN #############################

# trapezoidal rule:
step = np.zeros(1000)
trapint = np.zeros(1000)
for k in range(1000, 2000):
	num = k/100
	xvar1 = np.linspace(0.0, pi, num)
	h1 = pi/num

	trapint1 = 0.0
	for i in range(0,num-1):
		trapint1 += h1*(sin(xvar1[i+1]) + sin(xvar1[i]))/2.0
	
	trapint[k-1000] = trapint1
	step[k-1000] = num

trapintx = np.zeros(1000)
for k in range(1000, 2000):
	num = k/100
	xvar1 = np.linspace(0.0, pi, num)
	h1 = pi/num

	trapint1 = 0.0
	for i in range(0,num-1):
		trapint1 += h1*(xvar1[i+1]*sin(xvar1[i+1]) + xvar1[i]*sin(xvar1[i]))/2.0
	
	trapintx[k-1000] = trapint1
	step[k-1000] = num

# Simpson's rule:
simpint = np.zeros(1000)
for k in range(1000, 2000):
	num = k/100
	xvar1 = np.linspace(0.0, pi, num)
	h1 = pi/num

	simpint1 = 0.0
	for i in range(0,num-1):
		simpint1 +=  h1*(sin(xvar1[i+1]) + sin(xvar1[i]) + 4.0*sin((xvar1[i+1]+xvar1[i])/2.0))/6.0
	
	simpint[k-1000] = simpint1
	
simpintx = np.zeros(1000)
for k in range(1000, 2000):
	num = k/100
	xvar1 = np.linspace(0.0, pi, num)
	h1 = pi/num

	simpint1 = 0.0
	for i in range(0,num-1):
		simpint1 += h1*(xvar1[i+1]*sin(xvar1[i+1]) + xvar1[i]*sin(xvar1[i]) + 4.0*((xvar1[i+1]+xvar1[i])/2.0)*sin((xvar1[i+1]+xvar1[i])/2.0))/6.0
	
	simpintx[k-1000] = simpint1


plt. plot(step, abs(trapint-2.0), 'g',  linewidth=2, label = 'trapezoidal')
plt. plot(step, abs(simpint-2.0), 'm',  linewidth=2, label = 'Simpson')
plt. plot(step, abs(trapintx-pi), 'g',  linewidth=2, linestyle ='--')
plt. plot(step, abs(simpintx-pi), 'm',  linewidth=2,  linestyle ='--')

plt.xlabel('N (number of subintervals)', fontsize= 17)
plt.ylabel('absolute error', fontsize= 17)
pylab.legend(loc=1)
plt.show()



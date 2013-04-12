#!/usr/bin/python

import numpy as np
import pylab 
import matplotlib.pyplot as plt
from math import pi,pow, sqrt




def flin(a, x):
	b = apmag[a]
	c = (apmag[a+1]-apmag[a])/(time[a+1]-time[a])
	return b + c*(x-time[a])

def fquad(a, x):
	b = (x- time[a+1])*(x - time[a+2])*apmag[a]/((time[a]- time[a+1])*(time[a]-time[a+2]))
	c = (x - time[a])*(x - time[a+2])*apmag[a+1]/((time[a+1]-time[a])*(time[a+1]-time[a+2]))
	d = (x - time[a])*(x - time[a+1])*apmag[a+2]/((time[a+2]-time[a])*(time[a+2]-time[a+1]))
	return b+c+d



# MAIN ##################################

time = [0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0]
apmag = [0.302, 0.185, 0.106, 0.093, 0.24, 0.579, 0.561, 0.468, 0.302]
num = 10

for t in range(0, 8):
	xvar = np.linspace(time[t], time[t+1], num)
	plin = np.zeros(num)
	for k in range (0, num):
		plin[k] = flin(t, xvar[k])
	if t ==0:
		plt. plot(xvar, plin, 'b', linewidth=2, label='Linear')
	else:
		plt. plot(xvar, plin, 'b', linewidth=2)
	
for t in range(0, 7):
	xvar = np.linspace(time[t], time[t+1], num)
	pquad = np.zeros(num)
	for k in range (0, num):
		pquad[k] = fquad(t, xvar[k])
	plt. plot(xvar, pquad, 'g', linewidth=2)

xvar = np.linspace(time[7], time[7+1], num)
pquad = np.zeros(num)
for k in range (0, num):
	pquad[k] = fquad(6, xvar[k])
plt. plot(xvar, pquad, 'g', linewidth=2, label='Quadratic')

plt. plot(time, apmag, 'r.', marker='o', ms= 8)
plt.xlabel('time', fontsize= 17)
plt.ylabel('apparent magnitude', fontsize= 17)
pylab.legend(loc=2)
plt.show()


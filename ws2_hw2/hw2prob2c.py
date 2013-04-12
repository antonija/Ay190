#!/usr/bin/python

import numpy as np
import pylab 
import matplotlib.pyplot as plt
from math import pi,pow, sqrt
import scipy.interpolate

def psi0(z):
	return z*z*(2.0*z-3.0)+1.0

def psi1(z):
	return z*(1.0 -2.0*z + z*z)

def derf(a):
	x1 = apmag[a+1]
	x2 = apmag[a]
	h = time[a+1] - time[a]
	return (x1-x2)/h

def derf2(a):
	x1 = apmag[a]
	x2 = apmag[a-1]
	h = time[a] - time[a-1]
	return (x1-x2)/h

def derf3(a):
	x1 = apmag[a+1]
	x2 = apmag[a-1]
	h = time[a] - time[a-1]
	return (x1-x2)/(2.0*h)

def herm(a, x):
	z = (x-time[a])/(time[a+1]-time[a])
	term1 = apmag[a]*psi0(z)
	term2 = apmag[a+1]*psi0(1.0-z)
	term3 = derf(a)*(time[a+1]-time[a])*psi1(z)
	term4 = derf(a+1)*(time[a+1]-time[a])*psi1(1.0-z)
	return term1+term2+term3-term4

def herm2(a, x):
	z = (x-time[a])/(time[a+1]-time[a])
	term1 = apmag[a]*psi0(z)
	term2 = apmag[a+1]*psi0(1.0-z)
	term3 = derf2(a)*(time[a+1]-time[a])*psi1(z)
	term4 = derf2(a+1)*(time[a+1]-time[a])*psi1(1.0-z)
	return term1+term2+term3-term4

def herm3(a, x):
	z = (x-time[a])/(time[a+1]-time[a])
	term1 = apmag[a]*psi0(z)
	term2 = apmag[a+1]*psi0(1.0-z)
	term3 = derf3(a)*(time[a+1]-time[a])*psi1(z)
	term4 = derf3(a+1)*(time[a+1]-time[a])*psi1(1.0-z)
	return term1+term2+term3-term4

# MAIN ##################################

time = [0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0]
apmag = [0.302, 0.185, 0.106, 0.093, 0.24, 0.579, 0.561, 0.468, 0.302]
num = 10

# Hermite	
for t in range(1, 7):
	xvar = np.linspace(time[t], time[t+1], num)
	hermpol = np.zeros(num)
	for k in range (0, num):
		hermpol[k] = herm3(t, xvar[k])
	plt. plot(xvar, hermpol, 'g', linewidth=2)

xvar = np.linspace(time[0], time[1], num)
hermpol = np.zeros(num)
for k in range (0, num):
	hermpol[k] = herm(0, xvar[k])
plt. plot(xvar, hermpol, 'g', linewidth=2)

xvar = np.linspace(time[7], time[8], num)
hermpol = np.zeros(num)
for k in range (0, num):
	hermpol[k] = herm2(7, xvar[k])
plt. plot(xvar, hermpol, 'g', linewidth=2, label='Hermite')

# Spline
xvar2 = np.linspace(0.0, 1.0, 100)
spl = scipy.interpolate.splrep(time, apmag)
yvar = scipy.interpolate.splev(xvar2,spl)

plt.plot(xvar2, yvar, 'm', linewidth=2, label='Spline')

plt. plot(time, apmag, 'b.', marker='o', ms= 8)
plt.xlabel('time', fontsize= 17)
plt.ylabel('apparent magnitude', fontsize= 17)
pylab.legend(loc=2)
plt.show()


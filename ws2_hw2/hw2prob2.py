#!/usr/bin/python

import numpy as np
import pylab 
import matplotlib.pyplot as plt
from math import pi,pow, sqrt


def L(x, j):
	prod = 1.0
	for i in range(0, 9):
		if (i != j):
			a = prod
			prod = a*(x-time[i])/(time[j]-time[i])
	return prod

def func(x):
	b = 0.0
	for i in range(0,9):
		b += apmag[i]*L(x, i)
	return b

# MAIN ##################################

time = [0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0]
apmag = [0.302, 0.185, 0.106, 0.093, 0.24, 0.579, 0.561, 0.468, 0.302]

num = 500

xvar = np.linspace(0.0, 1.0, num)
p = np.zeros(num)
for k in range (0, num):
	p[k] = func(xvar[k])


plt. plot(xvar, p, 'b', linewidth=2)
plt. plot(time, apmag, 'r.', marker='o', ms= 8)
plt.xlabel('time', fontsize= 17)
plt.ylabel('apparent magnitude', fontsize= 17)

plt.show()


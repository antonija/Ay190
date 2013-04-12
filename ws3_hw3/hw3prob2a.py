#!/usr/bin/python

import numpy as np
import pylab 
import matplotlib.pyplot as plt
from math import pi,pow, sqrt, sin, exp
from scipy import special as sp

def f(x):
	return x*x/(exp(x) + 1.0)

def g(x):
	return f(x)*exp(x)

def lag(n):
	[lagroots, lagweights] = sp.l_roots(n,0)
	integ = 0.0
	for j in range(0,n):
		integ +=lagweights[j]*g(lagroots[j])
	return integ

# MAIN ###################################

step = np.zeros(30)
integral1 = np.zeros(30)
for i in range(10, 40):
	step[i-10] = i
	integral1[i-10] = lag(i)

bestint = integral1[29]
print bestint

plt. plot(step, integral1, 'g',  linewidth=2)
plt.xlabel('N (number of nodes)', fontsize= 17)
plt.ylabel('Integral', fontsize= 17)

plt.show()

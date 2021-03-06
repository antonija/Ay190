#!/usr/bin/python

import numpy as np
import pylab 
import matplotlib.pyplot as plt
from math import pi,pow, sqrt, sin, exp
from scipy import special as sp

def f(x):
	return x*x/(exp(x) + 1.0)


def leg(n, a, b):
	[legroots, legweights] = sp.p_roots(n,0)
	integ = 0.0
	for i in range(0, n):
		troot = ((b-a)*legroots[i]+(a+b))/2.0
		integ += legweights[i]*f(troot)
	return integ*(b-a)/2.0


# MAIN ###################################


energbin = np.linspace(0.0,8.0, 32)
nume = np.zeros(32)
for i in range(0, 31):
	nume[i] = leg(20, energbin[i], energbin[i+1])
energy = energbin*20.0
total = sum(nume)
print total

spec = nume*100.0/total

plt.plot(energy, spec, 'g', linestyle='steps', linewidth=2)
plt.xlabel('energy', fontsize= 17)
plt.ylabel('percentage of electrons', fontsize= 17)

plt.show()



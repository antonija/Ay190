#!/usr/bin/python

import numpy as np
import pylab 
import matplotlib.pyplot as plt
from math import pi,pow, sqrt

def f(x):
	return x*x*x - 5.0*x*x + x

def derf(x):
	return 3.0*x*x - 10.0*x + 1.0
 

def forward(x, h):
	return (f(x) - f(x-h))/h


def central(x, h):
	return (f(x+h) - f(x-h))/(2.0*h)
	

	

# MAIN #############################


num = 10


xvar1 = np.linspace(-2.0, 6.0, num)
h1 = 8.0/num

abserr1f = np.zeros(num)
abserr1c = np.zeros(num)


for i in range(0,num):
	abserr1f[i] = forward(xvar1[i], h1) - derf(xvar1[i])
	abserr1c[i] = central(xvar1[i], h1) - derf(xvar1[i])


plt. plot(xvar1, abserr1f, 'm', linestyle='-', linewidth=2, label='forward, h1')
plt. plot(xvar1, abserr1c, 'g', linestyle='-', linewidth=2, label= 'central, h1')


xvar2 = np.linspace(-2.0, 6.0, 2*num)
h2 = 4.0/num

abserr2f = np.zeros(2*num)
abserr2c = np.zeros(2*num)


for i in range(0, 2*num):
	abserr2f[i] = forward(xvar2[i], h2) - derf(xvar2[i])
	abserr2c[i] = central(xvar2[i], h2) - derf(xvar2[i])


plt. plot(xvar2, abserr2f, 'm', linestyle='--', linewidth=2, label= 'forward, h2')
plt. plot(xvar2, abserr2c, 'g', linestyle='--', linewidth=2, label = 'central, h2')
pylab.legend(loc=3)
pylab.axhline(y=0.0, color='0.4', linestyle=':', linewidth=1.0)
plt.show()



#!/usr/bin/python

import numpy as np
import pylab 
import matplotlib.pyplot as plt
from math import pi,pow, sqrt, sin, exp, cos
from scipy import special as sp

def f(x):
	return 3.0*pow(x,5) + 5.0*pow(x,4) -1.0*pow(x,3)


def findroot(a,b):

	c = (a + b)/2.0
	f1a = f(a)
	f1b = f(b)
	f1c = f(c)
	#print fa, fb, fc
	numit = 0

	while (abs(f1c) > eps):

		if (((f1a > 0.0) and (f1c < 0.0)) or ((f1a < 0.0) and (f1c > 0.0))):
			b = c
		else:
			a = c		
		c = (a + b)/2.0
		f1a = f(a)
		f1b = f(b)
		f1c = f(c)	
			
		numit += 1
	#print numit
	
	return c
# MAIN ###################################

amax = 2.6
bmin = -3.2
eps = 1.0e-12

intervals = np.arange(bmin, amax, 0.1)

num = len(intervals)

for i in range(0, num-1):
	a = intervals[i+1]
	b = intervals[i]
	fa = f(a)
	fb = f(b)
	if (abs(fa) <= eps):
		print a
	elif (abs(fb) <= eps):
		print b
	elif ((fa > 0.0) and (fb < 0.0)) or ((fa < 0.0) and (fb > 0.0)):
		root = findroot(a,b)
		print root

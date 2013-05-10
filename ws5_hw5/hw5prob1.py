#!/usr/bin/python

import numpy as np
import pylab 
import matplotlib.pyplot as plt
from math import pi,pow, sqrt, sin, exp, cos
from scipy import special as sp

def f(x,t):
	return x - omega*t - e*sin(x)


def findroot(t):
	a = 0.0
	b = pi + omega*t
	c = (a + b)/2.0
	fa = f(a, t)
	fb = f(b, t)
	fc = f(c, t)
	#print fa, fb, fc
	numit = 0

	while (abs(fc) > eps):

		if (((fa > 0.0) and (fc < 0.0)) or ((fa < 0.0) and (fc > 0.0))):
			b = c
		else:
			a = c		
		c = (a + b)/2.0
		fa = f(a, t)
		fb = f(b, t)
		fc = f(c, t)	
			
		numit += 1
	print numit
	
	return c
# MAIN ###################################

T = 31558148.64
omega = 2.0*pi/T
majoraxis = 1.496e13
#e = 0.0167
e = 0.9999
minoraxis = majoraxis*sqrt(1.0 - e*e)
eps = 1.0e-10

time = [7862400.0, 15724800.0, 23587200.0]

xpos = np.zeros(3)
ypos = np.zeros(3)

for i in range(3):
	E = findroot(time[i])
	xpos[i] = majoraxis*cos(E)
	ypos[i] = minoraxis*sin(E)
	print 'E= ', E, '  x= ', xpos[i]/majoraxis, '  y= ',  ypos[i]/majoraxis 


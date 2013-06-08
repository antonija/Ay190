#!/usr/bin/env python
import sys,math
from pylab import *
import numpy


# Plotting errors


def analytic(x, time):
	x0prime = x0 + v*time # moving Gaussian peak
	return exp(-((x-x0prime)**2)/(2.0*sigma*sigma))

# Upwind method
def newyUPW(yold,t):
	new=zeros(n)
	new[0] = yold[0]
	for i in range(1,n):
		new[i] = yold[i]-v*dt*(yold[i]-yold[i-1])/dx
	return new


# Downwind method
def newyDNW(yold,t):
	new=zeros(n)
	new[0] = yold[0]
	for i in range(1,n-1):
		new[i] = yold[i]-v*dt*(yold[i+1]-yold[i])/dx
	new[n-1] = new[n-2]
	return new


# FTCS method 
def newyFTCS(yold,t):
	new=zeros(n)
	new[0] = yold[0]
	for i in range(1,n-1):
		new[i] = yold[i]-v*dt*(yold[i+1]-yold[i-1])/dx
	new[n-1] = new[n-2]
	return new


# Lax-Friedrich method 
def newyLF(yold,t):
	new=zeros(n)
	new[0] = yold[0]
	for i in range(1,n-1):
		new[i] = 0.5*(yold[i+1]+yold[i-1])-0.5*v*dt*(yold[i+1]-yold[i-1])/dx
	new[n-1] = new[n-2]
	return new

# leapfrog method 
def newyleap(yold,t):
	new=zeros(n)
	new[0] = yold[0]
	for i in range(1,n-1):
		new[i] =  yold[i]-v*dt*(yold[i+1]-yold[i-1])/dx
	new[n-1] = new[n-2]
	return new

# Lax-Wendroff method 
def newyLW(yold,t):
	new=zeros(n)
	new[0] = yold[0]
	for i in range(1,n-1):
		new[i] = yold[i]-0.5*v*dt*(yold[i+1]-yold[i-1])/dx
		new[i] += 0.5*v*v*dt*dt*(yold[i-1]-2*yold[i]+yold[i+1])/(dx*dx)
	new[n-1] = new[n-2]
	return new
 

# parameters
dx = 0.1
v = 0.4
# set up the grid here. Use a decent number of zones;
# perhaps to get a dx of 0.1
x = arange(0.0,100.0,dx)
n = len(x)
y = zeros(n)
cfl = 0.6
dt = cfl*dx/v
t = 0.0

# for initial data
sigma = sqrt(15.0)
x0 = 30.0

#set up initial conditions
y = analytic(x,t)

yold2 = y
yold = y
ntmax = 100

errarray = zeros(ntmax)
tarray = zeros(ntmax)

for it in range(ntmax):
    t = t + dt
    # save previous and previous previous data
    yold2 = yold
    yold = y

    # get new data; ideally just call a function
    y = newyLW(yold,t)

    # get analytic result for time t
    yana = analytic(x,t)
    # compute error estimage
    err = sum(abs(y - yana))
    errarray[it] = err
    tarray[it] = t

    print "it = ",it,err
 
plot(tarray, errarray, linewidth=3)
xlabel('time', size=18)
ylabel('error',size=18)
show()




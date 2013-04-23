#!/usr/bin/python

import sys
from scipy import *
from pylab import *
from math import pi,pow, sqrt

# global constants
ggrav = 6.67e-8
clite = 3.0e10
msun = 1.99e33


# EOS
# neutron stars:
# polyG = 2.0
# polyK = 100.0 * 5.55e38/6.1755e17**polyG

# EOS for 
# white dwarfs:
polyG = 4.0/3.0
polyK = 1.244e15*0.5**polyG

# central values
rhoc = 1.0e10

# minimum pressure
rhomin = 1.0
min_press = polyK*rhomin**polyG


# grid
rmax = 1.55e8


def set_grid(rmax,nzones):
    # set up the grid and return the
    # radius array and dr
	dr = rmax/nzones
	rad = arange(dr*0.0001,rmax,dr)
	return (rad,dr)

def tov_RHS(r,data):
	rhs = zeros(2)
	mass = data[1]
	press = max(data[0], min_press)
	rho = (press/polyK)**(1.0/polyG)
	rhs[0] = -ggrav*mass*rho/(r*r)
	rhs[1] = 4.0*pi*r*r*rho
	return rhs

def tov_RK2(old_data,r,dr):
	k1 = zeros(2)
	k2 = zeros(2)
	new_data = zeros(2)
	k1 = dr*tov_RHS(r, old_data)
	k2 = dr*tov_RHS(r+0.5*dr, old_data + 0.5*k1)
	new_data = old_data + k2
	return new_data
    
def tov_RK3(old_data,r,dr):
	k1 = zeros(2)
	k2 = zeros(2)
	k3 = zeros(2)
	new_data = zeros(2)
	k1 = dr*tov_RHS(r, old_data)
	k2 = dr*tov_RHS(r + 0.5*dr, old_data + 0.5*k1)
	k3 = dr*tov_RHS(r + dr, old_data - k1 + 2.0*k2)
	new_data = old_data + (k1 + 4.0*k2 + k3)/6.0
	return new_data


def tov_RK4(old_data,r,dr):
	k1 = zeros(2)
	k2 = zeros(2)
	k3 = zeros(2)
	k4 = zeros(2)
	new_data = zeros(2)
	k1 = dr*tov_RHS(r, old_data)
	k2 = dr*tov_RHS(r + 0.5*dr, old_data + 0.5*k1)
	k3 = dr*tov_RHS(r + 0.5*dr, old_data + 0.5*k2)
	k4 = dr*tov_RHS(r + dr, old_data + k3)
	new_data = old_data + (k1 + k4 + 2.0*(k2 + k3))/6.0
	return new_data


def tov_integrate(rmax,nzones):

    # set up grid
    (rad,dr) = set_grid(rmax,nzones)

    # initialize some variables
    tovdata = zeros((nzones,2))
    # 0 -- press
    # 1 -- mbary

    tovout = zeros((nzones,4))
    # 0 -- rho
    # 1 -- press
    # 2 -- eps
    # 3 -- mass
     
    # central values
    tovdata[0,0] = polyK * rhoc**polyG
    tovdata[0,1] = 0.0


    # you will need to track the surface (where press <= press_min)
    isurf = 0
    for i in range(nzones-1):

        # integrate one step using RK2 (can change this to
        # RK3 or RK4)
        tovdata[i+1,:] = tov_RK4(tovdata[i,:],rad[i],dr)

        # check if press below 0
        if(tovdata[i+1,0] <=min_press):
            isurf = i

        # press and mass
        tovout[i+1,1] = tovdata[i+1,0]

        if (i+1 > isurf and isurf > 0):
            tovout[i+1,3] = tovdata[isurf,1]
        else:
            tovout[i+1,3] = tovdata[i+1,1]
        # compute density
        tovout[i+1,0] = (tovout[i+1,1]/polyK)**(1.0/polyG)
        # compute eps
        tovout[i+1,2] = tovout[i+1,1]/((polyG-1.0)*tovout[i+1,0])

    return (tovout,isurf,dr)

# for convergence: 
# number of points
na = array([600, 800, 1000])
# to store masses
masses = zeros(len(na))
# to store the drs
drs = zeros(len(na))

for i in range(len(na)):
    (tov_star,isurf,dr) = tov_integrate(rmax,na[i])

    masses[i] = tov_star[na[i]-1, 3]
    drs[i] = dr

xvec=linspace(0.0,1.0,1000)
semilogx(xvec, tov_star[:,0]/tov_star[1,0], 'b', linewidth=2,  label = r'density, $\rho_{norm} = 10^{10}$ g cm$^{-1}$')
semilogx(xvec, tov_star[:,1]/tov_star[1,1], 'g', linewidth=2, label = r'pressure, $P_{norm} = 10^{28}$ dyn cm$^{-2}$')
semilogx(xvec, tov_star[:,3]/tov_star[999,3], 'r', linewidth=2, label= r'mass, $M_{norm} = 2.9\times 10^{33}$ g $= 1.45 M_{\odot}$')
xlabel('radius / 1550 km', fontsize= 17)
ylabel('normalized values', fontsize= 17)
ylim(-0.45, 1.0)
print tov_star[1,0], tov_star[1,1], tov_star[999,3]
Qs = (masses[2]-masses[1])/(masses[1]-masses[0])
Qs2 = (pow(drs[2],2) - pow(drs[1],2))/(pow(drs[1],2) - pow(drs[0],2))
Qs3 = (pow(drs[2],3) - pow(drs[1],3))/(pow(drs[1],3) - pow(drs[0],3))
Qs4 = (pow(drs[2],4) - pow(drs[1],4))/(pow(drs[1],4) - pow(drs[0],4))
print Qs
print Qs2, Qs3, Qs4
legend(loc=3)
show()

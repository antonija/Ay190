#!/usr/bin/env python


import numpy as np
import pylab
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from math import pi,acos, log, exp, pow, trunc, sin, cos, sqrt
from time import time
import struct

# mass inside radius x
def massin(x):
	return 4.0*pi*rho0*rs*rs*rs*(log((rs+x)/rs) - x/(rs + x))


# NFW profile
def density(x):
	if (x < rexp):
	
		dens1 = rho0*rs/(x*pow(1.0 + x/rs,2))
	else:
		dens1 = rho0*rs*exp(-x/rexp)/(x*pow(1.0 + x/rs,2))
	return dens1


# randomly assign coordinate theta 
#(the angle between radial vector and the z-axis)
def theta():
	u = np.random.uniform(-1.0, 1.0)
	#u = 2.0*np.random.rand() - 1.0
	theta1 = acos(u)
	return theta1

# randomly assign coordinate phi (azimuthal angle)
def phi():
	phi1 = 2.0*pi*np.random.rand()
	return phi1


# d(ln rho)/dr
def dersigma(sigma, x):
	derlnrho = (x/rs - 1.0)/(x*(1.0 +x/rs))
	return -G*massin(x)/(2.0*sigma*x*x) - sigma*derlnrho/2.0



# Plots: #############################################################

# Plot the initial positions in xy plane:
def plotxy():
	plt.plot(coordxyz[:,0], coordxyz[:,1], 'b.')
	plt.grid(True)
	plt.axes().set_aspect('equal')
	plt.xlabel('x')
	plt.ylabel('y')

# Plot the inital distribution of angles theta or phi:
def plotangles():
	f, axarr = plt.subplots(2, sharex=True)
	axarr[0].plot(coord[:,1]/pi, 'g.', ms=2)
	axarr[1].plot(coord[:,2]/pi, 'r.', ms=2)
	axarr[0].set_ylabel('$\Theta/\pi$', size=16)
	axarr[1].set_ylabel('$\phi/\pi$', size=16)
	axarr[1].set_xlabel('particle ID', size=16)
	plt.ylabel(r'$\theta$')

# Plot the number of enclosed particles and compare with
# the prediction from NFW formula
def plotenmass():
	plt.plot(rad/3.086e21, sumperbin/2.0e33, linewidth= 4, label='Generated distribution')
	plt.plot(rad/3.086e21, sumanal/2.0e33, 'r--',linewidth= 3, label='NFW profile')
	plt.xlabel('Radius [kpc]', size =18)
	plt.ylabel('Interior mass of particles [M$_{\odot}$]', size =18)
	plt.xlim(0,200.0)
	plt.legend(loc=4)


# Plot the density profile and compare with the NFW formula
def plotdenspro():
	plt.loglog(rad/3.086e21, nperbin*m1/(4.0*pi*rad*rad*delr), linewidth= 4, label='Generated profile')
	plt.loglog(rad/3.086e21, rho0*rs/(rad*(1.0+rad/rs)*(1.0+rad/rs)), 'r--',linewidth= 3, label='NFW profile')
	plt.xlabel('Radius [kpc]', size=18)
	plt.ylabel(r'$\rho(r)$ [g cm$^{-3}$]', size=18)
	plt.xlim(0,200.0)
	plt.legend(loc=3)

# Plot calculated velocity dispersion as a function of radius
def plotveldisp():
	plt.plot(radinv/3.086e+21, sigmavec*1.0e-5, linewidth= 3, color='black')
	plt.xlabel('radius [kpc]', size=18)
	plt.ylabel(r'$\sigma$ [km/s]', size=18)

# Plot generated velocity dispersion as a function of radius
def plotveldisp2():
	radaux=np.zeros(partid)
	velx=np.zeros(partid)
	vely=np.zeros(partid)
	velz=np.zeros(partid)
	for k in range(0, partid):
		velx[k] = velxyz[k,0]
		vely[k] = velxyz[k,1]
		velz[k] = velxyz[k,2]
		radaux[k] =coord[k,0]/3.086e+21	
	#plt.plot(radaux, velx, 'b.')
	rmean=np.zeros(100)
	vdisperx= np.zeros(100)
	vdispery= np.zeros(100)
	vdisperz= np.zeros(100)
	for j in range(0,100):	
		r=np.zeros(99)
		vx= np.zeros(99)
		vy= np.zeros(99)
		vz= np.zeros(99)
		for i in range(0, 99):
			r[i]=coord[j*90+i,0]
			vx[i]=velx[j*90+i]
			vy[i]=vely[j*90+i]
			vz[i]=velz[j*90+i]
		rmean[j]=np.mean(r)/3.086e+21	
		vdisperx[j]=np.std(vx)
		vdispery[j]=np.std(vy)
		vdisperz[j]=np.std(vz)
	plt.plot(rmean, vdisperx, 'r.')
	plt.plot(rmean, vdispery, 'b.')
	plt.plot(rmean, vdisperz, 'g.')
	plt.xlim(0, 200)

############################################################
# Parameters: ##############################################



rs = 20.0*3.086e+21 
rexp = 200.0*3.086e+21
rend = 250.0*3.086e+21 
rho0 = 4.546e-25
Mtot= 2.0e45
G = 6.67e-08


numpart = 10000 	# total number of particles
numbinspl1 = 1001 	# number of radial bins +1
delr = rend/numbinspl1 	# radial bin size
rad = np.arange(delr/2.0, rend-delr/2.0, delr) 	# radial bins from the center out
radinv = np.arange(rend-delr/2.0, delr/2.0, -delr)	# radial bins from outside in
nbins = len(rad) 	# number of radial bins
m1 = Mtot/numpart 	# mass per particle

nperbin = np.zeros(nbins)	# number of particles in each radial bin
sumperbin = np.zeros(nbins)	# sum of particles up to that bin
sumanal = np.zeros(nbins)	# sum of particles up to that predicted by analytic NFW profile
sigmabin = np.zeros(nbins)	# velocity dispersion at each radial bin



# Main: #####################################################
t_start = time()
np.random.seed()

coord = np.zeros([numpart,3]) 	# spherical coordinates
# angles:
for i in range(0, numpart):
	coord[i,1] = theta()
	coord[i,2] = phi()
# radius:
partid=0
for j in range(0, nbins):
	r = rad[j]
	deln = trunc(4.0*pi*r*r*density(r)*delr/m1)
	nperbin[j] = deln
	if (j == 0):
		sumperbin[0] = deln*m1
		#sumanal[0] = massin(r)
	else:
		sumperbin[j] = sumperbin[j-1] + deln*m1
		#sumanal[j] = sumanal[j-1] + massin(r)
	
	sumanal[j] = massin(r)
	for l in range(0, deln):
		coord[partid, 0] = rad[j] + delr*np.random.uniform(-0.5, 0.5)
		partid = partid + 1

print partid # actual number of particles within 250 kpc halo
pid = np.arange(0,numpart,1)

#for u in range(0,numpart):
#	print u, coord[u,0]

coordxyz = np.zeros([partid,3]) 	# Cartesian coordinates 
for k in range(0, partid):
	coordxyz[k,0] = coord[k,0]*sin(coord[k,1])*cos(coord[k,2])/3.086e+21 #in kpc
	coordxyz[k,1] = coord[k,0]*sin(coord[k,1])*sin(coord[k,2])/3.086e+21
	coordxyz[k,2] = coord[k,0]*cos(coord[k,1])/3.086e+21


# Velocity dispersion:
sigmafunc= odeint(dersigma, 1.0e+1, radinv)
sigmavec= sigmafunc[:,0]

# Initial velocities:
velxyz = np.zeros([partid,3])
numvel=0
for k in range(0, partid):
	for h in range(0, nbins):
		if (abs(coord[k,0]-radinv[h]) < delr/2.0):
			velxyz[k,0] = sigmavec[h]*np.random.randn()*1.0e-5 # in km/s
			velxyz[k,1] = sigmavec[h]*np.random.randn()*1.0e-5 
			velxyz[k,2] = sigmavec[h]*np.random.randn()*1.0e-5 
			numvel+=1
print numvel # Number of particles with assigned velocities


# Plotting #############################################################

plotveldisp()
plt.show()


# Printing initial conditions ###########################################

print 'Writing initial condition file...'
npart           = np.zeros(6,int)
massarr         = np.zeros(6,float)
timestart       = 0.0
redshift        = 0.0
flag_sfr        = 0
flag_feedback   = 0
nall            = np.zeros(6,int)
flag_cooling    = 0
num_files       = 1
box_size        = 0.0
Omega0          = 0.0
OmegaLambda     = 0.0
HubbleParam     = 0.0
flag_age        = 0
flag_metals     = 0
nall_HW         = np.zeros(6,int)
flag_entropy    = 0

npart[1] = partid
nall[1] = partid
massarr[1] = m1/2.0e+43

f = open('dm_halo.dat','wb')
fmt = '=iIIIIIIddddddddiiIIIIIIiiddddiiIIIIIIi60xi'
if struct.calcsize(fmt) != 256+8: print 'Gadget header format string is not 256 bytes'
header = struct.pack(fmt,
                     256,
                     npart[0],npart[1],npart[2],npart[3],npart[4],npart[5],
                     massarr[0],massarr[1],massarr[2],massarr[3],massarr[4],massarr[5],
                     timestart,redshift,flag_sfr,flag_feedback,
                     nall[0],nall[1],nall[2],nall[3],nall[4],nall[5],
                     flag_cooling,num_files,box_size,Omega0,OmegaLambda,HubbleParam,flag_age,flag_metals,
                     nall_HW[0],nall_HW[1],nall_HW[2],nall_HW[3],nall_HW[4],nall_HW[5],
                     flag_entropy,
                     256)
f.write(header)

f.write(struct.pack('=i',int(4*3*len(coordxyz))))
for n in range(partid):
    for i in range(3): f.write(struct.pack('=f',coordxyz[n,i]))
f.write(struct.pack('=i',int(4*3*len(coordxyz))))

f.write(struct.pack('=i',int(4*3*len(velxyz))))
for n in range(partid):
    for i in range(3): f.write(struct.pack('=f',velxyz[n,i]))
f.write(struct.pack('=i',int(4*3*len(velxyz))))

f.write(struct.pack('=i',int(4*len(pid))))
for n in range(partid): f.write(struct.pack('=i',pid[n]))
f.write(struct.pack('=i',int(4*len(pid))))

t_end = time()
print t_end-t_start,'s for '+str(partid)+' particles'



#!/usr/bin/python

import numpy as np
import pylab 
import matplotlib.pyplot as plt

import time
import sys

# Gauss-Jordan elimination code written by Jarno Elonen (from http://elonen.iki.fi/code/misc-notes/python-gaussj/)

def gauss_jordan(m, eps = 1.0/(10**10)):
  """Puts given matrix (2D array) into the Reduced Row Echelon Form.
     Returns True if successful, False if 'm' is singular.
     NOTE: make sure all the matrix items support fractions! Int matrix will NOT work!
     Written by Jarno Elonen in April 2005, released into Public Domain"""
  (h, w) = (len(m), len(m[0]))
  for y in range(0,h):
    maxrow = y
    for y2 in range(y+1, h):    # Find max pivot
      if abs(m[y2][y]) > abs(m[maxrow][y]):
        maxrow = y2
    (m[y], m[maxrow]) = (m[maxrow], m[y])
    if abs(m[y][y]) <= eps:     # Singular?
      return False
    for y2 in range(y+1, h):    # Eliminate column y
      c = m[y2][y] / m[y][y]
      for x in range(y, w):
        m[y2][x] -= m[y][x] * c
  for y in range(h-1, 0-1, -1): # Backsubstitute
    c  = m[y][y]
    for y2 in range(0,y):
      for x in range(w-1, y-1, -1):
        m[y2][x] -=  m[y][x] * m[y2][y] / c
    m[y][y] /= c
    for x in range(h, w):       # Normalize row y
      m[y][x] /= c
  return True

def solveGE(M, b):
  """
  solves M*x = b
  return vector x so that M*x = b
  :param M: a matrix in the form of a list of list
  :param b: a vector in the form of a simple list of scalars
  """
  m2 = [row[:]+[right] for row,right in zip(M,b) ]
  return [row[-1] for row in m2] if gauss_jordan(m2) else None

# end of J. Elonen's code



for i in range(0, 5):
	print ' '
	print 'matrix no. ', i+1
	mat = np.loadtxt('LSE'+str(i+1)+'_m.dat')
	bvec = np.loadtxt('LSE'+str(i+1)+'_bvec.dat')
	sh = np.shape(mat)
	det = np.linalg.det(mat)
	print 'shape: ', sh
	print 'determinant: ', det
	start1 = time.time()
	solveGE(mat, bvec)
	end1 = time.time()
	print 'GE time (in seconds): ', end1-start1
	
	start2 = time.time()
	s2 = np. linalg.solve(mat, bvec)
	end2 = time.time()
	print 'numpy.linalg.solve time (in seconds): ', end2-start2


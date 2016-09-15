"""
Minimal *working* script for parameter estimation of Lotka-Volterra model
"""

import pylab as plt
import numpy as np
import tellurium as te
from scipy.optimize import differential_evolution
import time

from fractions import gcd
from functools import reduce

import os
os.chdir('/Users/d3m629/Google Drive/tellurium_exercises')



# equdist_timepoints()
# A helper function that maps inequidistantly 
# spaced time points onto equally spaced array
#==============================================================================
# input:
#   1) original array
#   2) number of significant digits for time values
# return tuple of 4 items:
#   1) starting time
#   2) end time
#   3) numer of points
#   4) index of points that correspond to original timepoints array

# rounding to a predefined number of significant figures.
def round_sigfig(x, sig=2):
    sigfig = sig - np.floor(np.log10(x)) - 1
    sigfig[sigfig == np.inf] = 0 # no need to round zero
    xOut = [np.around(i[0], int(i[1])) for i in np.column_stack((x, sigfig))]
    return np.array(xOut)

def equdist_timepoints(x, sig):
    x = round_sigfig(x, sig)
    minStep = reduce(gcd, np.ediff1d(x))
    numPoints = int(np.ptp(x)/minStep) + 1
    xNew = np.linspace(np.min(x), np.max(x), num = numPoints)
    idx = [np.where(xNew == i)[0].tolist()[0] for i in x]
    return (min(x), max(x), numPoints, idx)
#==============================================================================





# Previously simulated data. Pretending as it is experimental.
#==============================================================================
ed = np.genfromtxt("lv_nontransformed.txt", skip_header = 1, 
                   missing_values = '', delimiter = '\t',
                   usecols=(0,1,2,3,4))

plt.plot(ed[:,0], ed[:,1:])
# aliases
timepts  = ed[:,0]
a_data = ed[:,1]
x_data = ed[:,2]
y_data = ed[:,3]
b_data = ed[:,4]
#==============================================================================



# The Lotka-Volterra model
#==============================================================================
r = te.loada("""
   A + X -> 2 X; k1*A*X;
   X + Y -> 2 Y; k2*X*Y;
   Y -> B; k3*Y;
   
   A = 1; X = 1e-1; Y = 1e-1; B = 1e-6;
   k1 = 0.3; k2 = 3; k3 = 0.1;
   
""")
#==============================================================================




# objective function
#==============================================================================
# Compute the simulated data and calculate the value of objective function
# This one will take some work since the data is all relative
def objectiveFunc(p):
    r.resetToOrigin()  
    k1, k2, k3, A, X, Y, B = p
    r.model['k1'] = k1
    r.model['k2'] = k2
    r.model['k3'] = k3
    r.model['A'] = A
    r.model['X'] = X 
    r.model['Y'] = Y
    r.model['B'] = B
    # core - the simulation
    try:
        s = r.simulate(timing[0], timing[1], timing[2])
        s = s[timing[3],1:]
    except RuntimeError:
        return 1e6    
    # error residuals
    a = a_data - s[:, 0]
    x = x_data - s[:, 1]
    y = y_data - s[:, 2]
    b = b_data - s[:, 3]        
    return np.nansum(a*a) + np.nansum(x*x) + np.nansum(y*y) + np.nansum(b*b)
#==============================================================================






# testing the timings
#==============================================================================
bounds = [(0.,5.), (0.,5.), (0.,5.), (0.,5.), (0.,5.), (0.,5.), (0.,5.)]
tic = time.time()
timing = equdist_timepoints(timepts, 2) # we'll make it global
result = differential_evolution(objectiveFunc, bounds)
toc = time.time()
print(dict(zip(['k1','k2','k3','A','X','Y','B'], result.x)))
print(toc - tic)
# 77, 34, 31, 31, 60, 34 seconds
#==============================================================================








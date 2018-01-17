from pyIsochrone import pymodels
import numpy as np
import matplotlib.pyplot as plt

# Call the fortran code models.f from python and properly save the data 
V,by,L,Teff,g,M,logM,R,Nout = pymodels(0,830,0.03,4300,'try.dat',3)
V = V[:Nout] # V
by = by[:Nout] # b-y
L = L[:Nout] # log(L/Lsun)
Teff = Teff[:Nout] # log(Teff)
g = g[:Nout] # log g
M = M[:Nout] # M/Msun
logM = logM[:Nout] # log(M/Msun)
R = R[:Nout] # R/Rsun

# plot the results
plt.plot(by,V)
plt.show()

# load stars.lst
stars = np.loadtxt('stars_std.lst')

# remove the blue stragglers and foreground/background stars from the data
## CODE MISSING

# minimize the MSE between data and isochrone in order to find optimum age, metalicity, reddening and distance


# print the final values

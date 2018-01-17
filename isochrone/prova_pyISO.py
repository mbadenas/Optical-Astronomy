from pyIsochrone import pymodels
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import scipy.optimize as opt
from cycler import cycler

# Call the fortran code models.f from python and properly save the data 
def call_pymodels(z,dist,eby,age,model):
    V,by,L,Teff,g,M,logM,R,Nout = pymodels(z,dist,eby,age,model)
    V = V[:Nout] # V
    by = by[:Nout] # b-y
    L = L[:Nout] # log(L/Lsun)
    Teff = Teff[:Nout] # log(Teff)
    g = g[:Nout] # log g
    M = M[:Nout] # M/Msun
    logM = logM[:Nout] # log(M/Msun)
    R = R[:Nout] # R/Rsun
    return V,by,L,Teff,g,M,logM,R

V,by,L,Teff,g,M,logM,R = call_pymodels(0,830,0.03,4300,3)
# plot the results
plt.figure()
plt.subplot(121)
plt.plot(by,V,'+')

# load stars.lst
stars = np.loadtxt('stars_std.lst')

# remove the blue stragglers and foreground/background stars from the data
stars_filtered = stars[stars[:,2]>0.25,:] # remove all stars with b-y<0.25
stars_filtered = stars_filtered[(stars_filtered[:,2]*10+11.5)>stars_filtered[:,0],:] # remove all stars over the line V=10(b-y)+11.5
stars_filtered = stars_filtered[(stars_filtered[:,2]*(-7.9)+14.7)<stars_filtered[:,0],:] # remove all stars below the line V=-7.9(b-y)+14.7
Vs = stars_filtered[:,0]
bys = stars_filtered[:,2]

# Show the filtered data and the real data to visually check the filtering makes sense
plt.plot(bys, Vs,'.',label='Observations',color='k',markersize=2)
plt.subplot(122)
plt.plot(stars[:,2], stars[:,0],'.',label='Observations',color='k',markersize=2)
plt.plot(by,V,'+')
plt.show()


# minimize the MSE between data and isochrone in order to find optimum age, metalicity, reddening and distance
Mones = np.ones((len(stars_filtered[:,2]),2000))
def mse_iso(z,dist,reddening,age,model):
    V,by,L,Teff,g,M,logM,R = call_pymodels(z,dist,reddening,age,model)
    N = len(V)
    Viso_Vs = ((Mones[:,:N]*V).transpose()-Vs)**2
    byiso_bys = ((Mones[:,:N]*by).transpose()-bys)**2
    msemat = Viso_Vs+byiso_bys
    msevec = np.min(msemat,axis=0)
    return sum(msevec)

def mse_weight(z,dist,reddening,age,model):
    V,by,L,Teff,g,M,logM,R = call_pymodels(z,dist,reddening,age,model)
    N = len(V)
    Viso_Vs = ((Mones[:,:N]*V).transpose()-Vs)**2/Vs
    byiso_bys = ((Mones[:,:N]*by).transpose()-bys)**2/bys
    msemat = Viso_Vs+byiso_bys
    msevec = np.min(msemat,axis=0)
    return sum(msevec)



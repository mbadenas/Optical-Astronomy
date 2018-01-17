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

# load stars.lst
stars = np.loadtxt('stars_std.lst')
Vs = stars[:,0]
bys = stars[:,2]

# plot isochrones and data for many combinations
fcolors=['g','b','r','c']
linestyle=['-','--','-.',':',(0,(5, 1.5, 1, 1.5, 1, 1.5))]
zvec = np.array([0,0.01])
ebyvec = np.array([0.03,0.035,0.04])
dvec = np.arange(700,901,100)
dvec = [700,900]
agevec = np.arange(4200,4701,100)
agevec = [4100,4300,4500,4700]
mod_names = ['Cassisi (overshoot)','Cassisi (canonical)','Padova']
leg_models = []
leg_dist = []
for k,name in enumerate(mod_names):
    leg_models.append(mlines.Line2D([], [], linestyle='-', color=fcolors[k],
        linewidth=1,label=name))
for j,dist in enumerate(dvec):
    leg_dist.append(mlines.Line2D([], [], linestyle=linestyle[j], color='k',
                linewidth=1,label='d=%g' %dist))
for z in zvec:
    for eby in ebyvec:
        fig = plt.figure()
        fig.suptitle('z=%.1g  eby=%.2g' %(z,eby))
        for i,age in enumerate(agevec):
            ax = fig.add_subplot(2,2,i+1)
            ax.plot(bys, Vs,'.',label='Observations',color='k',markersize=2)
            ax.set_xlabel('b-y')
            ax.set_ylabel('V')
            ax.set_title('age=%g' %age)
            first_legend=ax.legend(handles=leg_dist, handlelength=3, loc='right')
            axl = ax.add_artist(first_legend)
            ax.legend(handles=leg_models,loc='lower left')
            for j,dist in enumerate(dvec):
                for k in range(3):
                    V,by,L,Teff,g,M,logM,R = call_pymodels(z,dist,eby,age,k+1)
                    ax.plot(by,V,ls=linestyle[j],color=fcolors[k],lw=.8)
            ax.set_ylim([19,9])
            ax.set_xlim([-0.05,1.1])
        fig.tight_layout()
plt.show()







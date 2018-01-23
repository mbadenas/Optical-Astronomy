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
param_names = ['age','d','[Fe/H]','E(b-y)']
key = [3,1,2,0]
params = [4200, 830, 0.0, 0.03]
delta_params = [600, 70, 0.01, 0.02]
params2 = [[params[i]-delta_params[i], params[i]+delta_params[i]] for i in range(len(params))]
agevec = np.arange(4200,4701,100)
agevec = [4100,4300,4500,4700]
mod_names = ['Cassisi (overshoot)','Cassisi (canonical)','Padova']
mod_keys = [3]
leg_models = []
for k,model in enumerate(mod_keys):
    leg_models.append(mlines.Line2D([], [], linestyle='-', color=fcolors[k],
        linewidth=1,label=mod_names[model-1]))

fig=plt.figure(1)
ax3 = fig.add_subplot(2,2,3)
ax1 = fig.add_subplot(2,2,1,sharex=ax3)
ax4 = fig.add_subplot(2,2,4,sharey=ax3)
ax2 = fig.add_subplot(2,2,2,sharey=ax1,sharex=ax4)
ax = [ax1,ax2,ax3,ax4]
for i,el in enumerate(params):
    leg = []
    ax[i].plot(bys, Vs,'.',label='Observations',color='k',markersize=2)
    for j,elem in enumerate(params2[i]):
        leg.append(mlines.Line2D([], [], linestyle=linestyle[j], color='k',
                linewidth=1,label='%s=%g' %(param_names[i],elem)))
        for k,model in enumerate(mod_keys):
            inputs = [params[key[0]],params[key[1]],params[key[2]],params[key[3]],model]
            inputs[key[i]]=elem
            print inputs
            V,by,L,Teff,g,M,logM,R = call_pymodels(*inputs)
            ax[i].plot(by,V,ls=linestyle[j],color=fcolors[k],lw=.8)
    first_legend=ax[i].legend(handles=leg, loc='right')
    axl = ax[i].add_artist(first_legend)
    ax[i].legend(handles=leg_models,loc='lower left')

plt.setp(ax1.get_xticklabels(),visible=False)
ax1.set_ylabel('V')
plt.setp(ax2.get_xticklabels(),visible=False)
plt.setp(ax2.get_yticklabels(),visible=False)
ax3.set_xlabel('b-y')
ax3.set_ylabel('V')
ax4.set_xlabel('b-y')
plt.setp(ax4.get_yticklabels(),visible=False)
ax1.set_ylim([19,9])
ax4.set_xlim([-0.05,1.1])
ax3.set_ylim([19,9])
ax3.set_xlim([-0.05,1.1])
fig.tight_layout()
plt.show()







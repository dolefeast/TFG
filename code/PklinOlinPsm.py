import util_tools
import pandas as pd
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from pathlib import Path
plt.rc('lines', linewidth=2.5)
#matplotlib.use('pgf') #Saves the output as pgf
plt.rcParams['axes.unicode_minus'] = False #Latex format stuff
plt.rcParams.update({
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False
})

#files = Path('/home/santi/TFG/outputs_santi/class_output/').glob('*pk.dat')
pklin = Path('/home/santi/TFG/outputs_santi/linspace_class').glob('*pk*069*')
psmooth = Path('/home/santi/TFG/outputs_santi/linspace_class').glob('*psm*069*')
Olin = Path('/home/santi/TFG/outputs_santi/linspace_class').glob('*Olin*069*')

pklin, params = util_tools.many_files(list(pklin))
psmooth, params = util_tools.many_files(list(psmooth))
Olin, params = util_tools.many_files(list(Olin))

h = 0.676
kmax = 0.51
fontsize = 28

for i, data in enumerate([pklin, psmooth, Olin]):
    fig, ax = plt.subplots()
    ax.spines['top'].set_visible(False), ax.spines['right'].set_visible(False)
    ax.set_xscale('log'), ax.set_yscale('log')
    fig.set_tight_layout(True)
    idx = np.where(data[0]<=kmax)
    x, y = np.array(data[0]), np.array(data[1])
    x = x[idx]
    y = y[idx]

    ax.plot(x, y, color='teal')
    ax.set_xlabel(r'$k [h $Mpc$^{-1}] $', fontsize=fontsize) 
    if i==0:
        ax.set_ylabel(r'$P(k) [$Mpc$^3 h^{-3}]$', fontsize=fontsize)
        plt.savefig('../figs/Pklin.png')
    elif i==1:
        ax.set_ylabel(r'$P_{smooth}(k) [$Mpc$^3 h^{-3}]$', fontsize=fontsize)
        plt.savefig('../figs/Psm.png')
    else:
        ax.set_ylabel(r'$O_{lin}(k)$', fontsize=fontsize)
        plt.savefig('../figs/Olin.png')

#    ax.plot(psmooth[0],psmooth[1], color='teal')
#    ax.plot(Olin[0],Olin[1], color='teal')

    #plt.savefig('../figs/PkOlPsm.pdf')
   

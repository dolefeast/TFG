import util_tools
import pandas as pd
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from pathlib import Path
#plt.style.use('fivethirtyeight')
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
kmin = 0.02
kmax = 0.51
fontsize = 28

for i, data in enumerate([pklin, psmooth, Olin]):
    fig, ax = plt.subplots()
    ax.spines['top'].set_visible(False), ax.spines['right'].set_visible(False)
    ax.set_xscale('log'), ax.set_yscale('log')
    fig.set_tight_layout(True)
    idx = np.where(np.logical_and(data[0]<=kmax, data[0]>=kmin))
    x, y = np.array(data[0]), np.array(data[1])
    x = x[idx]
    y = y[idx]

    ax.plot(x, y, color='teal')
#    ax.axhline(y=min(y), color='black', linewidth=1.3, alpha=0.7)
#    ax.axvline(x=kmin, color='black', linewidth=1.3, alpha=0.7)
    ax.set_xlabel(r'$k [h $Mpc$^{-1}] $', fontsize=fontsize) 
    ax.set_xticks((0.05, 1e-1, 0.45))
    ax.set_xticklabels((0.05, 1e-1, 0.45), fontsize=fontsize)
    yticks = np.linspace(min(y), max(y), 4)

    if i==0:
        yticks = [int(x) for x in np.round(yticks, -2)]
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticks, fontsize=fontsize/1.3)
        ax.set_ylabel(r'$\log_{10} P(k) [$Mpc$^3 h^{-3}]$', fontsize=fontsize)
        plt.savefig('../figs/Pklin.png')
    elif i==1:
        yticks = [int(x) for x in np.round(yticks, -2)]
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticks, fontsize=fontsize/1.3)
        ax.set_ylabel(r'$\log_{10} P_{smooth}(k) [$Mpc$^3 h^{-3}]$', fontsize=fontsize)
        plt.savefig('../figs/Psm.png')
    else:
        ax.set_ylabel(r'$\log_{10} O_{lin}(k)$', fontsize=fontsize)
        plt.savefig('../figs/Olin.png')

#    ax.plot(psmooth[0],psmooth[1], color='teal')
#    ax.plot(Olin[0],Olin[1], color='teal')

    #plt.savefig('../figs/PkOlPsm.pdf')
   

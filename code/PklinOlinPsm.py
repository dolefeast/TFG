import util_tools
import pandas as pd
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from pathlib import Path
plt.style.use('fivethirtyeight')
plt.rc('lines', linewidth=1.7)
import matplotlib
#matplotlib.use('pgf') #Saves the output as pgf
matplotlib.rcParams['axes.unicode_minus'] = False #Latex format stuff
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
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


for i, data in enumerate([pklin, psmooth, Olin]):
    fig, ax = plt.subplots()
    fig.set_tight_layout(True)
    idx = np.where(data[0]<=kmax)
    x, y = np.array(data[0]), np.array(data[1])
    x = x[idx]
    y = y[idx]
    
    

    ax.plot(x, y, color='teal')
    ax.set_xlabel(r'$k [h $Mpc$^{-1}] $') 
    if i==0:
        ax.set_ylabel(r'$P(k) [$Mpc$^3 h^{-3}]$')
        plt.savefig('../figs/Pklin.png')
    elif i==1:
        ax.set_ylabel(r'$P_{smooth}(k) [$Mpc$^3 h^{-3}]$')
        plt.savefig('../figs/Psm.png')
    else:
        ax.set_ylabel(r'$O_{lin}(k)$')
        plt.savefig('../figs/Olin.png')

#    ax.plot(psmooth[0],psmooth[1], color='teal')
#    ax.plot(Olin[0],Olin[1], color='teal')

    #plt.savefig('../figs/PkOlPsm.pdf')
    #plt.show()

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
pklin = list(Path('/home/santi/TFG/outputs_santi/linspace_class').glob('*pk*069*'))
psmooth = list(Path('/home/santi/TFG/outputs_santi/linspace_class').glob('*psm*069*'))
Olin = list(Path('/home/santi/TFG/outputs_santi/linspace_class').glob('*Olin*069*'))


print('Opening files...')
pklin, params = util_tools.many_files(pklin)
psmooth, params = util_tools.many_files(psmooth)
Olin, params = util_tools.many_files(Olin)
h = 0.676

fig, ax = plt.subplots(1, 3)
fig.tight_layout()

ax[0].plot(pklin[0],pklin[1]*Olin[1], color='teal')
[axes.set_xlabel(r'$k [h Mpc^{-1}] $') for axes in ax]
ax[0].set_ylabel(r'$P_{lin}(k) [Mpc^3 h^{-3}]$')
ax[1].set_ylabel(r'$P_{smooth}(k) [Mpc^3 h^{-3}]$')
ax[2].set_ylabel(r'$O_{lin}(k)$')

ax[1].plot(psmooth[0],psmooth[1], color='teal')
ax[2].plot(Olin[0],Olin[1], color='teal')

plt.savefig('../figs/PkOlPsm.pdf')
plt.show()

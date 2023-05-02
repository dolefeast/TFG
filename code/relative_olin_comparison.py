import util_tools
import pandas as pd
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from pathlib import Path
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
files = list(Path('/home/santi/TFG/outputs_santi/linspace_class').glob('*Olin*69*'))
hector_files = list(Path('/home/santi/TFG/lrg_eboss/model/').glob('*Olin*'))
hector_df, _ = util_tools.many_files(hector_files)      #all pk smooth as outputs from class

df, params = util_tools.many_files(files)
h = 0.676
fontsize = 20
color = 'teal'

fig, ax = plt.subplots(1, 1, figsize=(5*4/3, 5))

k_in, pk_in = df[0], df[1]
y = pk_in#util_tools.remove_bao(k_in, pk_in)
kmin, kmax = 0.0, 0.51
idx = np.where(np.logical_and(k_in<=kmax, k_in>=kmin))
x, y = np.array(k_in), np.array(pk_in)
x = x[idx]
y = y[idx]
ax.plot(x, y, color=color, label='Flat $O_{lin}(k)$')

k_hector, pk_hector = hector_df[0], hector_df[1]
ax.plot(k_hector, pk_hector, color='coral', label='Non flat $O_{lin}(k)(k)$')
ax.legend(loc='best')


x1, y1 = k_hector, pk_hector
ymod=np.interp(x, x1, y1)
#ymod = ymod(x2)
ax.set_xlabel(r'k [h Mpc$^{-1}$]', fontsize=fontsize)
ax.set_ylabel(r'$O_{lin}(k)$', fontsize=fontsize)
#    plt.plot(x2, (1-0.665)*y2, label='My data')
#    plt.plot(x1, y1, label='Original Y')
#plt.savefig('../figs/Olin_relative_comparison.pdf')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
yticks = ax.get_yticks()
ylabel = ax.get_ylabel()
ax.set_yticks(yticks)
ax.set_yticklabels([round(tick, 2) for tick in yticks], fontsize=fontsize/1.3)

xticks = ax.get_xticks()
xlabel = ax.get_xlabel()
ax.set_xticks(xticks)
ax.set_xticklabels([round(tick, 2) for tick in xticks], fontsize=fontsize/1.3)
fig.tight_layout()
plt.savefig('../figs/olin_comparison.png')
plt.show()

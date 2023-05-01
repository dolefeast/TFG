import matplotlib.pyplot as plt
import numpy as np
import util_tools
from pathlib import Path
import pandas as pd
#plt.style.use('fivethirtyeight')
plt.rc('lines', linewidth=1.7)
import matplotlib
#matplotlib.use('pgf') #Saves the output as pgf
matplotlib.rcParams['axes.unicode_minus'] = False #Latex format stuff
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
    'figure.facecolor': 'white'
})

files = list(Path('/home/santi/TFG/outputs_santi/rustico/').glob('Power_*NGC*69*'))
#file_name = util_tools.select_files(files)
kmin, kmax = (0., 0.45)

df_list, param_list = util_tools.many_files(files)  

fig, ax = plt.subplots(1, 1, figsize=(16, 9))
ax.set_xscale('log'), ax.set_yscale('log')
ax.spines['top'].set_visible(False), ax.spines['right'].set_visible(False)
fontsize = 29
#ax1.set_xlim((0,0.16)), ax2.set_xlim((0,0.16))

for data, params in zip(df_list, param_list):
    kk, P0, P2 = data[1], data[2], data[3] 
    ax.plot(kk, P0, 'x', markersize=17, color='teal')
    #ax2.plot(kk, kk*P2, '--', linewidth=0.8)

ax.set_ylabel('$P(k)[ h^{-3}$Mpc$^3]$', fontsize=fontsize)#, ax2.set_title('$k*P_2(k)$')
ax.set_xlabel('$k[h$Mpc$^{-1} ]$', fontsize=fontsize)#, ax2.set_title('$k*P_2(k)$')
ax.set_xticks((0.01, 0.05, 1e-1, 0.45))
ax.set_xticklabels((0.01, 0.05, 1e-1, 0.45), fontsize=fontsize)
yticks = (2600, 12500, 25000, 50000, 75000)
ax.set_yticks(yticks)
ax.set_yticklabels(yticks, fontsize=fontsize/1.3)
fig.set_tight_layout(True)
plt.savefig('../figs/Pkrustico.png')



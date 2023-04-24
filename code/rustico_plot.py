import matplotlib.pyplot as plt
import numpy as np
import util_tools
from pathlib import Path
import pandas as pd
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

files = list(Path('/home/santi/TFG/outputs_santi/rustico/').glob('Power_*NGC*'))
#file_name = util_tools.select_files(files)
kmin, kmax = (0., 0.45)

df_list, param_list = util_tools.many_files(files)  

fig, ax = plt.subplots(1, 1, figsize=(16, 9))
ax.set_xscale('log'), ax.set_yscale('log')
ax.spines['top'].set_visible(False), ax.spines['right'].set_visible(False)
#ax1.set_xlim((0,0.16)), ax2.set_xlim((0,0.16))

for data, params in zip(df_list, param_list):
    kk, P0, P2 = data[1], data[2], data[3] 
    ax.plot(kk, P0, 'x')
    #ax2.plot(kk, kk*P2, '--', linewidth=0.8)

ax.set_ylabel('$P(k)[ h^{-3}$Mpc$^3]$')#, ax2.set_title('$k*P_2(k)$')
ax.set_xlabel('$k[h$Mpc$^{-1} ]$')#, ax2.set_title('$k*P_2(k)$')
plt.show()


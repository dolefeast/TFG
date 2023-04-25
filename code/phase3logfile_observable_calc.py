import matplotlib as mpl
from pathlib import Path
import re
import util_tools
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as ct
import scipy as sp
#plt.style.use('fivethirtyeight')
import matplotlib
#matplotlib.use('pgf') #Saves the output as pgf
matplotlib.rcParams['axes.unicode_minus'] = False #Latex format stuff
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
})

zmax = 0.698
H0 = 67.6
rs = 147.784

files = list(Path('/home/santi/TFG/outputs_santi/phase3/logfiles_phase3').glob('*2nd*'))
#calculate_observables(*params[0])
Ok_list = []
a_para = []
a_perp = []

for data in files:
    out = util_tools.alpha_from_logfile(data)
    Ok_list.append(out[0])
    a_para.append(out[1])
    a_perp.append(out[2])

n_points = 500
Ok_min, Ok_max = min(Ok_list), max(Ok_list)
Ok_cont = np.linspace(Ok_min, Ok_max, n_points)
Ok_rang = Ok_max - Ok_min
H = lambda z, Ok, Om=0.31: H0*np.sqrt(Om*(1+z)**3 + Ok*(1+z)**2 + 1-Ok-Om)
DH_fid = lambda z, Ok: ct.c/1000/H(z, Ok)
DC_fid = lambda z, Ok: sp.integrate.quad(DH_fid, 0, z, args=(Ok,))[0] 
#DC_fid = np.array([sp.integrate.quad(DH_fid, 0, zmax, args=(ok,))[0] for ok in Ok_cont])

def DA(z, Ok):
    DC = DC_fid(z, Ok)
    DH = DH_fid(z, Ok)
    if Ok>0:
        k =  DH/np.sqrt(Ok)
        return k*np.sinh(np.sqrt(Ok)*DC/DH)
    elif Ok<0:
        k =  DH/np.sqrt(np.abs(Ok))
        return k*np.sin(np.sqrt(np.abs(Ok))*DC/DH)
    elif not Ok:
        return DH

#---Changing z->d (phase 2)
#DA_fid = np.array([DA(zmax, Ok, DC) for Ok, DC in zip(Ok_cont, DC_fid)])
#---Changing template (phase 3)
DA_fid = np.array([DA(zmax, Ok) for Ok in Ok_cont])

fig, axes = plt.subplots(3, 2, sharex=True, figsize=(10, 7))

for ax in axes.ravel():
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

elinewidth=1.7
capsize=3
capthick=1.5
color = 'teal'
fontsize = 20
for Ok, apara, aperp in zip(Ok_list, a_para, a_perp):
    axes[0,0].errorbar(Ok, apara[0], yerr=apara[1], fmt='x', 
                 elinewidth=elinewidth, capsize=capsize, capthick=capthick,
                color=color)
    axes[0,1].errorbar(Ok, aperp[0], yerr=aperp[1], fmt='x', 
                 elinewidth=elinewidth, capsize=capsize, capthick=capthick,
                color=color)
    axes[2,0].errorbar(Ok, DH_fid(zmax, Ok)*apara[0]/rs,  yerr=DH_fid(zmax, Ok)*apara[1]/rs, fmt='x', 
                 elinewidth=elinewidth, capsize=capsize, capthick=capthick, 
                 color=color) 
    #Changing z->d
    #axes[2,1].errorbar(Ok, DA(zmax, Ok)*aperp[0]/rs,  yerr=DA(zmax, Ok)*aperp[1]/rs, fmt='x', 
    #             elinewidth=elinewidth, capsize=capsize, capthick=capthick, 
    #             color=color) 
    #Changing template 
    axes[2,1].errorbar(Ok, DA(zmax, 0)*aperp[0]/rs,  yerr=DA(zmax, 0)*aperp[1]/rs, fmt='x', 
                 elinewidth=elinewidth, capsize=capsize, capthick=capthick, 
                 color=color) 
#    idx = max(-1+int(n_points*(Ok-Ok_min)/Ok_rang), 0)
#    axes[2,1].errorbar(Ok, DA_fid[idx]*aperp[0]/rs,  yerr=DA_fid[idx]*aperp[1]/rs, fmt='x', 
#                 elinewidth=elinewidth, capsize=capsize, capthick=capthick, 
#                 color=color) 


axes[1,0].plot(*Ok_cont, np.ones(Ok_cont)*DH_fid(zmax, 0)/rs, color=color) #Multiply by 0 is phase3
axes[1,1].plot(*Ok_cont, DA(zmax, 0)/rs, color=color)
axes[0,0].set_ylabel(r'$\alpha_{\parallel}$', fontsize=fontsize), axes[0,1].set_ylabel(r'$\alpha_{\perp}$', fontsize=fontsize)
axes[1,0].set_ylabel(r'$\left[ D_H/r_s\right]_{fid}$', fontsize=fontsize), axes[1,1].set_ylabel(r'$\left[ D_A/r_s\right]_{fid}$', fontsize=fontsize)
axes[2,0].set_ylabel(r'$D_H/r_s$', fontsize=fontsize), axes[2,1].set_ylabel(r'$D_A/r_s$', fontsize=fontsize)
axes[2,0].set_xlabel(r'$\left[ \Omega_k\right]^{fid\, 2}$', fontsize=fontsize), axes[2,1].set_xlabel(r'$\left[ \Omega_k\right]^{fid \,2}$', fontsize=fontsize)
plt.tight_layout()
plt.savefig('/home/santi/TFG/figs/phase3_DA_DH_flat.pdf')
plt.show()

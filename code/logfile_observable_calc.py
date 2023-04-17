import matplotlib as mpl
from pathlib import Path
import re
import util_tools
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as ct
import scipy as sp
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


zmax = 0.698
H0 = 67.6
rs = 147.784

files = list(Path('/home/santi/TFG/outputs_santi/phase2/logfiles_phase2').glob('*2nd*'))
#calculate_observables(*params[0])
Ok_list = []
a_para = []
a_perp = []

for data in files:
    out = util_tools.alpha_from_logfile(data)
    Ok_list.append(out[0])
    a_para.append(out[1])
    a_perp.append(out[2])

fig, ((ax1, ax2),(ax12, ax22),(ax13, ax23)) = plt.subplots(3, 2, sharex=True, figsize=(10, 7))
H = lambda z, Ok, Om=0.31: H0*np.sqrt(Om*(1+z)**3 + Ok*(1+z)**2 + 1-Ok-Om)
DH_fid = lambda z, Ok: ct.c/1000/H(z, Ok)
n_points = 500
Ok_min, Ok_max = min(Ok_list), max(Ok_list)
Ok_cont = np.linspace(Ok_min, Ok_max,n_points)
rang = Ok_max - Ok_min
DC_fid = np.array([sp.integrate.quad(DH_fid, 0, zmax, args=(ok,))[0] for ok in Ok_cont])

def DA(z, Ok, DC, DH):
    k =  DH/np.sqrt(np.abs(Ok))
    if Ok>0:
        return k*np.sinh(np.sqrt(Ok)*DC/DH)
    elif Ok<0:
        return k*np.sin(np.sqrt(np.abs(Ok))*DC/DH)
    elif not Ok:
        return DH

DA_fid = np.array([DA(zmax, Ok, DC, ct.c/H0) for Ok, DC in zip(Ok_cont, DC_fid)])

elinewidth=1.7
capsize=3
capthick=1.5
color = 'teal'
for Ok, apara, aperp in zip(Ok_list, a_para, a_perp):
    ax1.errorbar(Ok, apara[0], yerr=apara[1], fmt='x', 
                 elinewidth=elinewidth, capsize=capsize, capthick=capthick,
                color=color)
    ax2.errorbar(Ok, aperp[0], yerr=aperp[1], fmt='x', 
                 elinewidth=elinewidth, capsize=capsize, capthick=capthick,
                color=color)
    ax13.errorbar(Ok, DH_fid(zmax, Ok)*apara[0]/rs,  yerr=DH_fid(zmax, Ok)*apara[1]/rs, fmt='x', 
                 elinewidth=elinewidth, capsize=capsize, capthick=capthick, 
                 color=color) 
    idx = max(-1+int(n_points*(Ok-Ok_min)/rang), 0)
    ax23.errorbar(Ok, DA_fid[idx]*aperp[0]/rs,  yerr=DA_fid[idx]*aperp[1]/rs, fmt='x', 
                 elinewidth=elinewidth, capsize=capsize, capthick=capthick, 
                 color=color) 


ax12.plot(Ok_cont, DH_fid(zmax, Ok_cont)/rs, color=color)
ax22.plot(Ok_cont, DA_fid/rs, color=color)
ax1.set_ylabel(r'$\alpha_{\parallel}$'), ax2.set_ylabel(r'$\alpha_{\perp}$')
ax12.set_ylabel(r'$\left[ D_H/r_s\right]_{fid}$'), ax22.set_ylabel(r'$\left[ D_A/r_s\right]_{fid}$')
ax13.set_ylabel(r'$D_H/r_s$'), ax23.set_ylabel(r'$D_A/r_s$')
ax13.set_xlabel(r'$\left[ \Omega_k\right]^{fid\, 1}$'), ax23.set_xlabel(r'$\left[ \Omega_k\right]^{fid \,1}$')
plt.tight_layout()
plt.savefig('/home/santi/TFG/figs/phase2_DA_DH_flat.pdf')
plt.show()

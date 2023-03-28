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

files_santi = list(Path('/home/santi/TFG/outputs_santi/phase2/logfiles_phase2').glob('*3rd*'))
files_hector = list(Path('/home/santi/TFG/lrg_eboss/output/').glob('*'))

#calculate_observables(*params[0])
Ok_list_santi = []
a_para_santi = []
a_perp_santi = []
Ok_list_hector = []
a_para_hector = []
a_perp_hector = []

for data in files_santi:
    out = util_tools.alpha_from_logfile(data)
    Ok_list_santi.append(out[0])
    a_para_santi.append(out[1])
    a_perp_santi.append(out[2])

for data in files_hector:
    out = util_tools.alpha_from_logfile(data)
    Ok_list_hector.append(out[0])
    a_para_hector.append(out[1])
    a_perp_hector.append(out[2])


fig, ((ax1, ax2),(ax12, ax22),(ax13, ax23)) = plt.subplots(3, 2, sharex=True, figsize=(10, 7))
H = lambda z, Ok, Om=0.31: H0*np.sqrt(Om*(1+z)**3 + Ok*(1+z)**2 + 1-Ok-Om)
DH_fid = lambda z, Ok: ct.c/1000/H(z, Ok)
n_points = 500
Ok_cont = np.linspace(min(Ok_list_hector),max(Ok_list_hector),n_points)
DA_fid = np.array([sp.integrate.quad(DH_fid, 0, zmax, args=(ok,))[0] for ok in Ok_cont])
elinewidth=1
capsize=3
capthick=1.5
for Ok, apara, aperp in zip(Ok_list_santi, a_para_santi, a_perp_santi):
    ax1.errorbar(Ok, apara[0], yerr=apara[1], fmt='xb', 
                 elinewidth=elinewidth, capsize=capsize, capthick=capthick)
    ax2.errorbar(Ok, aperp[0], yerr=aperp[1], fmt='xb', 
                 elinewidth=elinewidth, capsize=capsize, capthick=capthick)
    ax13.errorbar(Ok, DH_fid(zmax, Ok)*apara[0]/rs,  yerr=DH_fid(zmax, Ok)*apara[1]/rs, fmt='xb', 
                 elinewidth=elinewidth, capsize=capsize, capthick=capthick) 
    idx = max(-1+int(n_points*(Ok+0.15)/0.3), 0)
    ax23.errorbar(Ok, DA_fid[idx]*aperp[0]/rs,  yerr=DA_fid[-1+int(100*(Ok+0.15)/0.3)]*aperp[1]/rs, fmt='xb', 
                 elinewidth=elinewidth, capsize=capsize, capthick=capthick) 
for Ok, apara, aperp in zip(Ok_list_hector, a_para_hector, a_perp_hector):
    ax1.errorbar(Ok, apara[0], yerr=apara[1], fmt='xr', 
                 elinewidth=elinewidth, capsize=capsize, capthick=capthick)
    ax2.errorbar(Ok, aperp[0], yerr=aperp[1], fmt='xr', 
                 elinewidth=elinewidth, capsize=capsize, capthick=capthick)
    ax13.errorbar(Ok, DH_fid(zmax, Ok)*apara[0]/rs,  yerr=DH_fid(zmax, Ok)*apara[1]/rs, fmt='xr', 
                 elinewidth=elinewidth, capsize=capsize, capthick=capthick) 
    idx = max(-1+int(n_points*(Ok+0.15)/0.3), 0)
    ax23.errorbar(Ok, DA_fid[idx]*aperp[0]/rs,  yerr=DA_fid[-1+int(100*(Ok+0.15)/0.3)]*aperp[1]/rs, fmt='xr', 
                 elinewidth=elinewidth, capsize=capsize, capthick=capthick) 



fig.suptitle('Phase 2 comparison plot')
ax12.plot(Ok_cont, DH_fid(zmax, Ok_cont)/rs)
ax22.plot(Ok_cont, DA_fid/rs)
ax1.set_ylabel(r'$\alpha_{\parallel}$'), ax2.set_ylabel(r'$\alpha_{\perp}$')
ax12.set_ylabel(r'$\left[ DH/r_s\right]_{fid}$'), ax22.set_ylabel(r'$\left[ DA/r_s\right]_{fid}$')
ax13.set_ylabel(r'$DH/r_s$'), ax23.set_ylabel(r'$DA/r_s$')
ax13.set_xlabel(r'$\left[ \Omega_k\right]^{fid\, 1}$'), ax23.set_xlabel(r'$\left[ \Omega_k\right]^{fid \,1}$')
plt.tight_layout()
plt.savefig('/home/santi/TFG/figs/comparison_phase2_DA_DH_flat.pdf')
plt.show()



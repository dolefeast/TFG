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

files = list(Path('/home/santi/TFG/outputs_santi/phase2/logfiles_phase2').glob('*3rd*'))
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
Ok_cont = np.linspace(-0.15,0.15,n_points)
DA_fid = np.array([sp.integrate.quad(DH_fid, 0, zmax, args=(ok,))[0] for ok in Ok_cont])
elinewidth=1
capsize=3
capthick=1.5
for Ok, apara, aperp in zip(Ok_list, a_para, a_perp):
    ax1.errorbar(Ok, apara[0], yerr=apara[1], fmt='xb', 
                 elinewidth=elinewidth, capsize=capsize, capthick=capthick)
    ax2.errorbar(Ok, aperp[0], yerr=aperp[1], fmt='xb', 
                 elinewidth=elinewidth, capsize=capsize, capthick=capthick)
    ax13.errorbar(Ok, DH_fid(zmax, Ok)*apara[0]/rs,  yerr=DH_fid(zmax, Ok)*apara[1]/rs, fmt='xb', 
                 elinewidth=elinewidth, capsize=capsize, capthick=capthick) 
    idx = max(-1+int(n_points*(Ok+0.15)/0.3), 0)
    ax23.errorbar(Ok, DA_fid[idx]*aperp[0]/rs,  yerr=DA_fid[-1+int(100*(Ok+0.15)/0.3)]*aperp[1]/rs, fmt='xb', 
                 elinewidth=elinewidth, capsize=capsize, capthick=capthick) 


ax12.plot(Ok_cont, DH_fid(zmax, Ok_cont)/rs)
ax22.plot(Ok_cont, DA_fid/rs)
ax1.set_ylabel(r'$\alpha_{\parallel}$'), ax2.set_ylabel(r'$\alpha_{\perp}$')
ax12.set_ylabel(r'$\left[ DH/r_s\right]_{fid}$'), ax22.set_ylabel(r'$\left[ DA/r_s\right]_{fid}$')
ax13.set_ylabel(r'$DH/r_s$'), ax23.set_ylabel(r'$DA/r_s$')
ax13.set_xlabel(r'$\left[ \Omega_k\right]^{fid\, 1}$'), ax23.set_xlabel(r'$\left[ \Omega_k\right]^{fid \,1}$')
plt.tight_layout()
plt.savefig('/home/santi/TFG/figs/phase2_DA_DH_flat.pdf')
plt.show()



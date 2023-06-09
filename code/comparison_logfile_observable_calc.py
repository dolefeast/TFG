from __init__ import *

files_santi = list(Path('/home/santi/TFG/outputs_santi/class_Om031_OL069/noflat_olin').glob('*'))
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
    if abs(out[0]) == 0.2:
        print('Continued!')
        continue
    Ok_list_santi.append(out[0])
    a_para_santi.append(out[1])
    a_perp_santi.append(out[2])

for data in files_hector:
    out = util_tools.alpha_from_logfile(data)
    Ok_list_hector.append(out[0])
    a_para_hector.append(out[1])
    a_perp_hector.append(out[2])


fig, axes = plt.subplots(3, 2, sharex=True, figsize=(10, 7))
H = lambda z, Ok, Om=0.31: H0*np.sqrt(Om*(1+z)**3 + Ok*(1+z)**2 + 1-Ok-Om)
DH_fid = lambda z, Ok: ct.c/1000/H(z, Ok)
n_points = 500
Ok_cont = np.linspace(min(Ok_list_hector),max(Ok_list_hector),n_points)
DA_fid = np.array([sp.integrate.quad(DH_fid, 0, zmax, args=(ok,))[0] for ok in Ok_cont])
santi_color = 'teal'
hector_color = 'coral'

for Ok, apara, aperp in zip(Ok_list_santi, a_para_santi, a_perp_santi):
    axes[0,0].errorbar(Ok, apara[0], yerr=apara[1], fmt='x', 
                 elinewidth=elinewidth, capsize=capsize, capthick=capthick, color=santi_color)
    axes[0,1].errorbar(Ok, aperp[0], yerr=aperp[1], fmt='x', 
                 elinewidth=elinewidth, capsize=capsize, capthick=capthick, color=santi_color)
    axes[2,0].errorbar(Ok, DH_fid(zmax, Ok)*apara[0]/rs,  yerr=DH_fid(zmax, Ok)*apara[1]/rs, fmt='x', 
                 elinewidth=elinewidth, capsize=capsize, capthick=capthick, color=santi_color) 
    idx = max(-1+int(n_points*(Ok+0.15)/0.3), 0)
    axes[2,1].errorbar(Ok, DA_fid[idx]*aperp[0]/rs,  yerr=DA_fid[idx]*aperp[1]/rs, fmt='x', 
                 elinewidth=elinewidth, capsize=capsize, capthick=capthick, color=santi_color) 

for Ok, apara, aperp in zip(Ok_list_hector, a_para_hector, a_perp_hector):
    axes[0,0].errorbar(Ok, apara[0], yerr=apara[1], fmt='x', 
                 elinewidth=elinewidth, capsize=capsize, capthick=capthick, color=hector_color)
    axes[0,1].errorbar(Ok, aperp[0], yerr=aperp[1], fmt='x', 
                 elinewidth=elinewidth, capsize=capsize, capthick=capthick, color=hector_color)
    axes[2,0].errorbar(Ok, DH_fid(zmax, Ok)*apara[0]/rs,  yerr=DH_fid(zmax, Ok)*apara[1]/rs, fmt='x', 
                 elinewidth=elinewidth, capsize=capsize, capthick=capthick, color=hector_color) 
    idx = max(-1+int(n_points*(Ok+0.15)/0.3), 0)
    axes[2,1].errorbar(Ok, DA_fid[idx]*aperp[0]/rs,  yerr=DA_fid[idx]*aperp[1]/rs, fmt='x', 
                 elinewidth=elinewidth, capsize=capsize, capthick=capthick, color=hector_color) 
#fig.suptitle(r'Comparison of {\color{teal}flat} vs {\color[rgb]{1, 0.5, 0.31}no flat} $O_{lin}(k)$', fontsize='large')

axes[1,0].plot(Ok_cont, DH_fid(zmax, Ok_cont)/rs, color='teal')
axes[1,1].plot(Ok_cont, DA_fid/rs, color='teal')
axes[0,0].set_ylabel(r'$\alpha_{\parallel}$', fontsize=fontsize)
axes[0,1].set_ylabel(r'$\alpha_{\perp}$', fontsize=fontsize)
axes[0,0].plot([], [], 'x', color=santi_color, label='Flat $O_{lin}(k)$')
axes[0,0].plot([], [], 'x', color=hector_color, label='No flat $O_{lin}(k)$')
axes[0,0].legend(loc='best', fontsize=fontsize/1.5)
axes[1,0].set_ylabel(r'$\left[ D_H/r_d\right]^{r}$', fontsize=fontsize)
axes[1,1].set_ylabel(r'$\left[ D_M/r_d\right]^{r}$', fontsize=fontsize)
axes[2,0].set_ylabel(r'$D_H/r_d$', fontsize=fontsize)
axes[2,1].set_ylabel(r'$D_M/r_d$', fontsize=fontsize)
axes[2,0].set_xlabel(r'$\left[ \Omega_k\right]^{r}$', fontsize=fontsize)
axes[2,1].set_xlabel(r'$\left[ \Omega_k\right]^{r}$', fontsize=fontsize)
for ax in axes.ravel():
    ticks = ax.get_yticks()
    label = ax.get_ylabel()
    ax.set_yticks(ticks)
    ax.set_yticklabels([round(tick, 2) for tick in ticks], fontsize=fontsize)
    ticks = ax.get_xticks()
    label = ax.get_ylabel()
    ax.set_xticks(ticks)
    ax.set_xticklabels([round(tick, 2) for tick in ticks], fontsize=fontsize)

plt.tight_layout()
plt.savefig('/home/santi/TFG/figs/flatnoflat_DADH.pdf')
#plt.show()



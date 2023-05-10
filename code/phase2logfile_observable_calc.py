from __init__ import * 

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

n_points = 500
Ok_min, Ok_max = min(Ok_list), max(Ok_list)
Ok_cont = np.linspace(Ok_min, Ok_max, n_points)
Ok_rang = Ok_max - Ok_min
H = lambda z, Ok, Om=0.31: H0*np.sqrt(Om*(1+z)**3 + Ok*(1+z)**2 + 1-Ok-Om)

@util_tools.iterable_output
def DH_fid(z, Ok):
    return ct.c/1000/H(z, Ok)
@util_tools.iterable_output
def DC_fid(z, Ok):
    return sp.integrate.quad(DH_fid, 0, z, args=(Ok,))[0]

#DC_fid = np.array([sp.integrate.quad(DH_fid, 0, zmax, args=(ok,))[0] for ok in Ok_cont])

@util_tools.iterable_output
def DM(z, Ok):
    DC = DC_fid(z, Ok)
    DH = DH_fid(z, Ok)
    if Ok>0:
        k =  DH/np.sqrt(Ok)
        return k*np.sinh(np.sqrt(Ok)*DC/DH)
    elif Ok<0:
        k =  DH/np.sqrt(np.abs(Ok))
        return k*np.sin(np.sqrt(np.abs(Ok))*DC/DH)
    elif not Ok:
        return DC

#---Changing z->d (phase 2)
#DM_fid = np.array([DM(zmax, Ok, DC) for Ok, DC in zip(Ok_cont, DC_fid)])
#---Changing template (phase 3)

fig, axes = plt.subplots(3, 2, sharex=True, figsize=(10, 7))
elinewidth=1.7
capsize=3
capthick=1.5
color = 'teal'
fontsize = 20

DH_list = []
DM_list = []

for Ok, apara, aperp in zip(Ok_list, a_para, a_perp):
    current_DH = DH_fid(zmax, Ok) * apara[0]/rs
    current_DH_std = DH_fid(zmax, Ok) * apara[1]/rs
    current_DM = DM(zmax, Ok) * aperp[0]/rs
    current_DM_std = DM(zmax, Ok) * aperp[1]/rs
    axes[0,0].errorbar(Ok, apara[0], yerr=apara[1], fmt='x', 
                 elinewidth=elinewidth, capsize=capsize, capthick=capthick,
                color=color)
    axes[0,1].errorbar(Ok, aperp[0], yerr=aperp[1], fmt='x', 
                 elinewidth=elinewidth, capsize=capsize, capthick=capthick,
                color=color)
    axes[2,0].errorbar(Ok, current_DH,  yerr=current_DH_std, fmt='x', 
                 elinewidth=elinewidth, capsize=capsize, capthick=capthick, 
                 color=color) 
    axes[2,1].errorbar(Ok, current_DM,  yerr=current_DM_std, fmt='x', 
                 elinewidth=elinewidth, capsize=capsize, capthick=capthick, 
                 color=color) 
    DH_list.append((current_DH, current_DH_std))
    DM_list.append((current_DM, current_DM_std))
    print("""For Ok = {} 
          D_H/r_d = {} \\pm {} 
          D_M/r_d = {} \\pm {} """.format(Ok,*[round(x, 2) for x in (current_DH,current_DH_std,current_DM,current_DM_std)]))

axes[1,0].plot(Ok_cont, DH_fid(zmax, Ok_cont)/rs, color=color, linewidth=2) #Multiply by 0 is phase2
axes[1,1].plot(Ok_cont, DM(zmax, Ok_cont)/rs, color=color, linewidth=2)
axes[0,0].set_ylabel(r'$\alpha_{\parallel}$', fontsize=fontsize)
axes[0,1].set_ylabel(r'$\alpha_{\perp}$', fontsize=fontsize)
axes[1,0].set_ylabel(r'$\left[ D_H/r_d\right]^{r}$', fontsize=fontsize)
axes[1,1].set_ylabel(r'$\left[ D_M/r_d\right]^{r}$', fontsize=fontsize)
axes[2,0].set_ylabel(r'$D_H/r_d$', fontsize=fontsize)
axes[2,1].set_ylabel(r'$D_M/r_d$', fontsize=fontsize)
axes[2,0].set_xlabel(r'$\left[ \Omega_k\right]^{r}$', fontsize=fontsize)
axes[2,1].set_xlabel(r'$\left[ \Omega_k\right]^{r}$', fontsize=fontsize)
axes[2,0].set_xticks(Ok_list)
axes[2,0].set_xticklabels(Ok_list, fontsize=fontsize/1.3)
axes[2,1].set_xticklabels(Ok_list, fontsize=fontsize/1.3)
for ax in axes.ravel():
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ticks = ax.get_yticks()
    label = ax.get_ylabel()
    ax.set_yticks(ticks)
    ax.set_yticklabels([round(tick, 2) for tick in ticks], fontsize=fontsize/1.3)

plt.tight_layout()
plt.savefig('/home/santi/TFG/figs/phase2_DA_DH_flat.pdf')
plt.show()

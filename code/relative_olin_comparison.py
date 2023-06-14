from __init__ import *

#files = Path('/home/santi/TFG/outputs_santi/class_output/').glob('*pk.dat')
files = list(Path('/home/santi/TFG/outputs_santi/linspace_class').glob('*Olin*69*'))
df, params = util_tools.many_files(files)
df = df[0]

hector_files = list(Path('/home/santi/TFG/lrg_eboss/model/').glob('*Olin*'))
hector_df, _ = util_tools.many_files(hector_files)      #all pk smooth as outputs from class
hector_df = hector_df[0]

h = 0.676
fontsize = 20
color = 'teal'

fig, ax = plt.subplots(1, 1, figsize=(5*4/3, 5))

k_in, pk_in = df[0], df[1]
y = pk_in#util_tools.remove_bao(k_in, pk_in)
kmin, kmax = 0.028, 0.51
idx = np.where(np.logical_and(k_in<=kmax, k_in>=kmin))
x, y = np.array(k_in), np.array(pk_in)
x = x[idx]
y = y[idx]
ax.set_xscale('log')
ax.set_yscale('log')
ax.plot(x, y, color=color, label='Flat $O_{lin}(k)$')

k_hector, pk_hector = np.array(hector_df[0]), np.array(hector_df[1])
idx = np.where(np.logical_and(k_hector<=kmax, k_hector>=kmin))
k_hector = k_hector[idx]
pk_hector = pk_hector[idx]
ax.plot(k_hector, pk_hector, color='coral', label='Non flat $O_{lin}(k)$')
ax.legend(loc='best', fontsize=fontsize/1.5)



#x1, y1 = k_hector, pk_hector
#ymod=np.interp(x, x1, y1)
#ymod = ymod(x2)
ax.set_xlabel(r'k [h Mpc$^{-1}$]', fontsize=fontsize)
ax.set_ylabel(r'$O_{lin}(k)$', fontsize=fontsize)
#    plt.plot(x2, (1-0.665)*y2, label='My data')
#    plt.plot(x1, y1, label='Original Y')
#plt.savefig('../figs/Olin_relative_comparison.pdf')
yticks = ax.get_yticks()
ylabel = ax.get_ylabel()
ax.set_yticks(yticks, minor=True)
ax.set_yticklabels([round(tick, 2) for tick in yticks], fontsize=fontsize/1.3)

xticks = ax.get_xticks()
xlabel = ax.get_xlabel()
ax.set_xticks(xticks, minor=True)
ax.set_xticklabels([round(tick, 2) for tick in xticks], fontsize=fontsize/1.3)
ax.set_xlim((kmin,kmax))
ax.set_ylim((min(pk_hector)**1.7, max(pk_hector)**1.7))
fig.tight_layout()
plt.savefig('../figs/olin_comparison.pdf')
plt.show()

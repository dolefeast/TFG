from __init__ import *

#files = Path('/home/santi/TFG/outputs_santi/class_output/').glob('*pk.dat')
files = list(Path('/home/santi/TFG/outputs_santi/linspace_class').glob('*Olin*69*'))
hector_files = list(Path('/home/santi/TFG/lrg_eboss/model/').glob('*Olin*'))
hector_df, _ = util_tools.many_files(hector_files)      #all pk smooth as outputs from class

df, params = util_tools.many_files(files)
h = 0.676
fontsize = 20
color = 'teal'

fig, ax = plt.subplots(1, 2, figsize=(5*16/9, 5))

k_in, pk_in = df[0], df[1]
y = pk_in#util_tools.remove_bao(k_in, pk_in)
kmin, kmax = 0.028, 0.51
idx = np.where(np.logical_and(k_in<=kmax, k_in>=kmin))
x, y = np.array(k_in), np.array(pk_in)
x = x[idx]
y = y[idx]
ax[0].plot(x, y, color=color, label='Flat $O_{lin}(k)$')
ax[0].plot(x, y, color=color, label='Flat $O_{lin}(k)$')


k_hector, pk_hector = hector_df[0], hector_df[1]
ax[0].plot(k_hector, pk_hector, '--k', label='No flat $O_{lin}(k)(k)$')
ax[0].legend(loc='best')


x1, y1 = k_hector, pk_hector
ymod=np.interp(x, x1, y1)
#ymod = ymod(x2)
ax[1].plot(x, -(y-ymod)/ymod*100, color='teal')
ax[0].set_xlabel(r'k [h Mpc$^{-1}$]', fontsize=fontsize)
ax[0].set_ylabel(r'$O_{lin}(k)$', fontsize=fontsize)
ax[1].set_xlabel(r'k [h Mpc$^{-1}$]', fontsize=fontsize)
#    plt.plot(x2, (1-0.665)*y2, label='My data')
#    plt.plot(x1, y1, label='Original Y')
#plt.savefig('../figs/Olin_relative_comparison.pdf')
for axis in ax.ravel():
    yticks = axis.get_yticks()
    ylabel = axis.get_ylabel()
    axis.set_yticks(yticks)
    axis.set_yticklabels([round(tick, 2) for tick in yticks], fontsize=fontsize/1.3)

    xticks = axis.get_xticks()
    xlabel = axis.get_xlabel()
    axis.set_xticks(xticks)
    axis.set_xticklabels([round(tick, 2) for tick in xticks], fontsize=fontsize/1.3)
fig.tight_layout()
plt.show()

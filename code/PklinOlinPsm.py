from __init__ import *

#CLASS output that was linearly spaced and with the BAO removed.  
pklin = Path('/home/santi/TFG/outputs_santi/linspace_class').glob('*pk*069*') #The power spectrum
psmooth = Path('/home/santi/TFG/outputs_santi/linspace_class').glob('*psm*069*') #The smoothed power spectrum
Olin = Path('/home/santi/TFG/outputs_santi/linspace_class').glob('*Olin*069*') #The pure BAO

#Routine to read the data files used throughout this work
pklin, params = util_tools.many_files(list(pklin))
psmooth, params = util_tools.many_files(list(psmooth))
Olin, params = util_tools.many_files(list(Olin))

#Our k-region of interest
kmin, kmax = 0.02, 0.51

fontsize = 28 #Text fontsize

#Using for loop since there is code common to the three data sets
for i, data in enumerate([pklin[0], psmooth[0], Olin[0]]):
    fig, ax = plt.subplots()
    #ax.spines['top'].set_visible(False), ax.spines['right'].set_visible(False)
    ax.set_xscale('log'), ax.set_yscale('log')
    x, y = np.array(data[0]), np.array(data[1])
    idx = np.where(np.logical_and(x<=kmax, x>=kmin))
    x = x[idx]
    y = y[idx]

    ax.plot(x, y, color='teal', linewidth=3)

    logx = np.log10(x)
    logy = np.log10(y)
    xticks = np.logspace(min(logx), max(logx), 4, base=10)
    yticks = np.logspace(min(logy), max(logy), 4, base=10)

    if i==0:
        yticks = [int(x) for x in np.round(yticks, -2)]
        ylabel = r'$ P(k) [$Mpc$^3 h^{-3}]$' 
        name  = 'Pklin' #The filename
        print('First plot done!')
    elif i==1:
        yticks = [int(x) for x in np.round(yticks, -2)]
        ylabel = r'$ P_{smooth}(k) [$Mpc$^3 h^{-3}]$' 
        name = 'Psm'
        print('Second plot done!')
    elif i==2:
        yticks = [x for x in np.round(yticks, 3)]
        ylabel = r'$ O_{lin}(k)$'
        name = 'Olin'
        print('Third plot done!')

    ax.set_xticks([], minor=True)
    ax.set_xticks(xticks)
    ax.set_xticklabels(np.round(xticks, 2), fontsize=fontsize)

    ax.set_yticks([], minor=True)
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks, fontsize=fontsize)

    ax.set_xlabel(r'$k [h $Mpc$^{-1}] $', fontsize=fontsize) 
    ax.set_ylabel(ylabel, fontsize=fontsize)

    fig.set_tight_layout(True)
    plt.savefig(f'../figs/{name}.pdf')
   
plt.show()

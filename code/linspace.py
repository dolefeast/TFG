import util_tools
import pandas as pd
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from pathlib import Path

#TODO:
    #Given a directory full of Psmooth, return Olins and P_nobao.

#files = Path('/home/santi/TFG/outputs_santi/class_output/').glob('*pk.dat')
files = list(Path('/home/santi/TFG/outputs_santi/class_outputs').glob('*fixed_Ob_*pk.dat'))
#hector_files = list(Path('/home/santi/TFG/lrg_eboss/model/Olinkirkby_eboss_comb_z070_matterpower_15.txt').glob('Olin*'))
df, params = util_tools.many_files(files)      #all pk smooth as outputs from class

name_tag = input("Something to add to the filename?")
print('Creating Olin data files...')
for data, omega in zip(df, params):
    Om, OL, Ok = omega
    k_in, pk_in = np.array(data[0]), np.array(data[1])
    x2 = np.linspace(k_in[0], k_in[-1], len(k_in))
    y2 = np.interp(x2, k_in, pk_in)
    tag = f'_{name_tag}_linspace_Om0{int(100*Om)}_OL0{int(100*OL)}.txt'
    pklin_df = pd.DataFrame([x2.T, y2.T]).T
    print(f'\tPrinting to ...pklin_{tag}')
    pklin_df.to_csv(f'/home/santi/TFG/outputs_santi/linspace_class/pklin_{tag}', sep='\t', index=False, header=None)
    pk_nobao = util_tools.remove_bao(k_in, pk_in)
    pk_nobao = np.interp(x2, k_in, pk_nobao)
    pk_nobao_df = pd.DataFrame([x2.T, pk_nobao.T]).T
    print(f'\tPrinting to... psmlin_{tag}')
    pk_nobao_df.to_csv(f'/home/santi/TFG/outputs_santi/linspace_class/psmlin_{tag}', sep='\t', index=False, header=None)
    olin = util_tools.calculate_olin(k_in, pk_in)
    olin = np.interp(x2, k_in, olin)
    olin_df = pd.DataFrame([x2.T, olin.T]).T
    print(f'\tPrinting to ...Olin_{tag}')
    olin_df.to_csv(f'/home/santi/TFG/outputs_santi/linspace_class/Olin_{tag}', sep='\t', index=False, header=None)

plt.show()
print('Done!')

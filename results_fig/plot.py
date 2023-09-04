import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re


df_hector_fixed_Om = pd.read_csv('./results_fig5.txt', header=0, sep='\s+')#, names=col_names)
df_hector_fixed_OL = pd.read_csv('./results_fig6.txt', header=0, sep='\s+')#, names=col_names)
df_santi_fixed_OL = pd.read_csv('./logfilemcmc_phase2_run3_Om016_OL069.txt', header=0, sep='\s+')#, names=col_names)
df_santi_fixed_Om = pd.read_csv('./logfilemcmc_phase2_fixed_ratio_Om031_OL069_run3.txt', header=0, sep='\s+')#, names=col_names)

fig, ax = plt.subplots()

offset=0.001
ax.errorbar(df_hector_fixed_Om["Ok"], df_hector_fixed_Om["apara"], yerr=df_hector_fixed_Om["apara_err"], fmt='^', label='hector apara fixed Om')
ax.errorbar(df_hector_fixed_Om["Ok"]+offset, df_hector_fixed_Om["aperp"], yerr=df_hector_fixed_Om["aperp_err"], fmt='x', label='hector aperp fixed Om')
ax.errorbar(df_hector_fixed_OL["Ok"]+2*offset, df_hector_fixed_OL["apara"], yerr=df_hector_fixed_OL["apara_err"], fmt='o', label='hector apara fixed OL')
ax.errorbar(df_hector_fixed_OL["Ok"]+3*offset, df_hector_fixed_OL["aperp"], yerr=df_hector_fixed_OL["aperp_err"], fmt='*', label='hector aperp fixed OL')
ax.errorbar(df_santi_fixed_Om["Ok"]+4*offset, df_santi_fixed_Om["apara"], yerr=df_santi_fixed_Om["apara_err"], fmt='^', label='santi apara fixed Om')
ax.errorbar(df_santi_fixed_Om["Ok"]+5*offset, df_santi_fixed_Om["aperp"], yerr=df_santi_fixed_Om["aperp_err"], fmt='x', label='santi aperp fixed Om')
ax.errorbar(df_santi_fixed_OL["Ok"]+6*offset, df_santi_fixed_OL["apara"], yerr=df_santi_fixed_OL["apara_err"], fmt='o', label='santi apara fixed OL')
ax.errorbar(df_santi_fixed_OL["Ok"] + 7*offset, df_santi_fixed_OL["aperp"], yerr=df_santi_fixed_OL["aperp_err"], fmt='*', label='santi aperp fixed OL')

ax.legend(loc='best')
plt.show()

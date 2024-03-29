from pathlib import Path
import sys
import re
import util_tools
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as ct
import scipy as sp
import matplotlib
#matplotlib.use('pgf') #Saves the output as pgf
matplotlib.rcParams['axes.unicode_minus'] = False #Latex format stuff
FONTSIZE = 18
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
    'lines.linewidth': 2, 
    'legend.fontsize': FONTSIZE/1.8, 
    'font.size': FONTSIZE, 
    'xtick.labelsize': 'x-large', 
    'ytick.labelsize': 'x-large'
})

zmax = 0.698
H0 = 67.6
rs = 147.784

elinewidth=1.7
capsize=3
capthick=1.5


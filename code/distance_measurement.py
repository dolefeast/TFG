#!/usr/bin/env python3

#A script that runs an interactive plotter 
#in any dir.

import util_tools
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
import sys

#Normalize everything to its max value so the study 
#of the curves' behavious is easier
save = "-s" in sys.argv or "-sv" in sys.argv or "-vs" in sys.argv
view = "-v" in sys.argv or "-sv" in sys.argv or "-vs" in sys.argv

if save:
    fig_name = f'/home/santi/TFG/figs/{str(Path.cwd().parent.stem)}_{str(Path.cwd().stem)}_DH_DM_plot.pdf'
    print(fig_name)
else:
    fig_name = None

filtr = input('Give a filter... (Nothing or "*" for no filter)')
if filtr == '':
    filtr = '*'

files = list(Path('.').glob('*'+filtr+'*'))

util_tools.plot_DH_DM(files, fig_name=fig_name, save=save, view=view)

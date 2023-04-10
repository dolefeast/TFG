import matplotlib as mpl
from pathlib import Path
import util_tools
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as ct
import scipy as sp

files = Path('~/TFG/outputs_santi/linspace_class/').glob('*')
print(files)


df, _ = util_tools.many_files(files)      #all pk smooth as outputs from class


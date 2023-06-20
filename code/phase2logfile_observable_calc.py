from __init__ import * 

files = list(Path('/home/santi/TFG/outputs_santi/phase2/logfiles_phase2').glob('*3rd*'))

util_tools.plot_DH_DM(files, show=True)

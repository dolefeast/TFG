import numpy as np
import itertools
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import scipy as sp
import scipy.constants as ct
import re
from scipy.stats import norm

Om_flat = 0.31 
OL_flat = 0.69
zmax = 0.698
H0 = 67.6
h = H0/100
r_d = 147.784

def select_files(files):
    """
    Input: files
    Returns: file_name

    given a list of file names, allows the user to select the file in which they are interested. returns the name of such file.  the file names list is under the name files. 
    """
    for i, file_name in enumerate(files, start=1):
        print(f' [{i}]: \t {file_name}')

    while True:
        try:
            file_index = int(input('Choose the files to be opened... '))
            break
        except ValueError:
            print('That input was not valid')
        except IndexError:
            print('Too big of an input!')

    file_name = list(files)[file_index-1]

    print('\nAbriendo {}:\n'.format(file_name))
    return file_name

def open_file(file_name):
    """
    Input: file_name
    Returns: (pd.DataFrame, list)
    returns a Dataframe with the given file_name, and the names of the columns
    """
###    with open(file_name, 'r') as f:
#        lines = f.readlines()
#        for i, row in enumerate(lines):
#            if r'#' not in row:
#                break
#
#    names = lines[i-1][1:].split()
#    names_copy = []

#    for i, char in enumerate(names):
#        try:
#            int(char[0]) #The first character should always be 1:foo
#        except ValueError: #Means it doesn't have the structure 1:foo
#            names_copy[-1] += char #In which case it should be a part of the last element of names
#        else:
#            names_copy.append(char) #This means it had the structure 1:foo

    df = pd.read_csv(file_name, names=None, sep='\s+', header=None, comment='#')
    return df

def many_files(files, openfiles=None):
    """Parameters:
        files: A list of file names to chose the files from
    Returns: 
        df_list: a list of pandas dataframes made from the chosen files
        parameters: a list of the curvature parameters of each file
        
    Stops asking for files when inserting a blank space"""
    n_files = len(list(files))
    try:
        assert n_files > 0
    except AssertionError:
        raise AssertionError(f"There are no files in that directory.")

    df_list = []
    param_list = []
    if openfiles != None or n_files==1:
        p = re.compile('Om[0-9]*_O[mL][0-9]*')
        for file_name in files:
                parameters = p.findall(str(file_name))
                if len(parameters) == 0:
                    parameters = str(file_name)
                    print(f'\tOpening {parameters}...')
                param_list.append(get_params(parameters[0]))
                df_list.append(open_file(file_name))

        #return [df_list[0], param_list[0]]
        return [df_list, param_list]

    for i, file_name in enumerate(files, start=1):
        print(f' [{i}]: \t {file_name}')
    print('The exception handling is weak. Handle with care!')
    file_index = input('Choose the files to open("all" to add all files): ')
    while file_index != 'all':
        try:
            file_index = int(file_index)
            file_name = files[file_index-1]
            df_list.append(open_file(file_name))
            p = re.compile('Om[0-9]*_OL[0-9]*')
            parameters = p.findall(str(file_name))
            if len(parameters) == 0:
                parameters = str(file_name)
            param_list.append(get_params(parameters[0]))
        except ValueError:
            print('Opened [Om, OL]:', param_list)
            return df_list, param_list
        except UnboundLocalError:
            print(f'El archivo {parameters} está vacío')
            continue
        except IndexError:
            print(f'The index {file_index} is too big!')
            file_index = input('Choose next file to open: ')
            continue
        file_index = input('Choose next file to open: ')
    else:
        print('Opening all displayed files...')
        return many_files(files, openfiles='all')

def alpha_from_logfile(filename):
    """
    parameters:
        the pathfix of the logfile to be calculated
    returns 
        (Ok, (apara, aparasigma), (aperp, aperpsigma))
    """
    namestr = str(filename)
    if 'logfile' not in namestr:
        raise TypeError(f'{filename} was not a logfile! ')
    omegas = get_params(namestr)
    out = []
    out.append(omegas) #Append all the omegas as (Om OL Ok)


    with filename.open() as open_data:
        p = re.compile('[0-9].[0-9]*e\+?-?[0-9]*')
        indices = [9, 10, 61] #Indices for alpha parallel, alpha perp and chi2
        open_data = list(open_data)
        for i in indices:
            matches = p.findall(open_data[i])
            out.append([float(x) for x in matches])
    return out

def iterable_output(func):
    def wrapper(zmax, Om, OL):
        # Check if Om is iterable
        Ok = 1 - Om - OL
        result = []
        try:
            iter(Om) #Om is the array
            for om in Om:
                result.append(func(zmax, om, OL))
            return np.array(result)
        except TypeError:
            try:
                iter(OL) #OL is the array
                for ol in OL:
                    result.append(func(zmax, Om, ol))
                return np.array(result)
            except TypeError: #Both were scalars
                return func(zmax, Om, OL)
    
    return wrapper
def get_params(filename):
    """
    input:
        filename: the parameters in the name of each dataframe
    output:
        [Om, OL, Ok]: a list containing the parameters in the name
        """
    p = re.compile('O[mL][0-9]*')
    params = re.findall(p, filename)
    params = re.split('\D', "".join(params))
    omegas_list = []
    for i in params:
        try:
            omegas_list.append(int(i)/100) #Because the regex find is not done properly
        except ValueError:
            continue
    try:
        Om, OL = omegas_list
        Ok = 1 - Om - OL
    except ValueError: 
        Om, OL, Ok = [0.31, 0.69, 0.]
    phase = 0
    if 'phase2' in filename:
        Om = (Om, Om_flat)
        OL = (OL, OL_flat)
        phase = 2
    elif 'phase3' in filename:
        Om = (Om_flat, Om)
        OL = (OL_flat, OL)
        phase = 3
    elif 'phase4' in filename:
        Om = (Om, Om)
        OL = (OL, OL)
        phase = 4

    return [Om, OL, Ok, phase]


def calculate_rd(Omega_m):
    omega_nu = 0.06/93.14 * Omega_m/0.31
    omega_b = 0.0481 * h**2* Omega_m/0.31 
    omega_c = 0.2604 * h**2* Omega_m/0.31
    r_d = 55.154 * np.exp(-72.3 * (omega_nu + 0.0006)**2)/(omega_c + omega_b)**0.25351 / omega_b**0.12807
    return r_d

@iterable_output
def H(z, Om, OL):
    return H0*np.sqrt(Om*(1+z)**3 + (1-Om-OL)*(1+z)**2 + OL)
@iterable_output
def DH_fid(z, Om, OL):
    return ct.c/1000/H(z, Om, OL)
@iterable_output
def DC_fid(z, Om, OL):
    return sp.integrate.quad(DH_fid, 0, z, args=(Om, OL))[0]
@iterable_output
def DM_fid(z, Om, OL):
    DC = DC_fid(z, Om, OL)
    DH = DH_fid(z, Om, OL)
    Ok = 1 - Om - OL
    if Ok>0:
        k =  DH/np.sqrt(Ok)
        return k*np.sinh(np.sqrt(Ok)*DC/DH)
    elif Ok<0:
        k =  DH/np.sqrt(np.abs(Ok))
        return k*np.sin(np.sqrt(np.abs(Ok))*DC/DH)
    elif not Ok:
        return DC

def flip(items, ncol):
    return itertools.chain(*[items[i::ncol] for i in range(ncol)])

def plot_DH_DM(files, save=False, view=False, fig_name=None, n_points = 500, markersize = 10, elinewidth=3, capsize=5, capthick=3, fontsize = 20, linewidth = 3, color = 'teal', reduce_ticks = 2, calculate_chi = True, table = False):
    """Plots the 3 rows by 2 columns graphic of the brass outputs
    parameters:
        files: a pathlib list of LOGFILES 
        save=False: if save=True it saves the figure with the name  ../figs/parent/dir_name
        view=False: if view=True it shows the figure
        n_points=500: number of points with which it calculates the curve of continuous Ok
    """
    fid = [0, 0, 'c', 't', 'c, t']
    Om_list = [] #List of tuples (Om_ref, Om_fid)
    OL_list = []
    Ok_list = []
    a_para = []
    a_perp = []
    chi_list = []
    for data in files:
        out = alpha_from_logfile(data)
        try:
            Om_list.append(out[0][0]) 
            OL_list.append(out[0][1]) 
            Ok_list.append(out[0][2]) 
            phase = out[0][3]
            a_para.append(out[1])
            a_perp.append(out[2])
            chi_list.append(out[3])
        except TypeError:
            print(f'{data} was not a logfile!')
            continue

    fid_tag = fid[phase]
    Om_list = np.array(Om_list)
    OL_list = np.array(OL_list)
    Om_ref_cont = np.linspace(max(Om_list[:,0]), min(Om_list[:,0]), n_points) #This could be done smarter by just interpolating. 
    OL_ref_cont = np.linspace(max(OL_list[:,0]), min(OL_list[:,0]), n_points) #No need to check for maxs ans mins and blablabla
    Om_fid_cont = np.linspace(max(Om_list[:,1]), min(Om_list[:,1]), n_points)
    OL_fid_cont = np.linspace(max(OL_list[:,1]), min(OL_list[:,1]), n_points)

    Ok_min, Ok_max = min(Ok_list), max(Ok_list)
    Ok_cont = np.linspace(Ok_min, Ok_max, n_points)
    Ok_rang = Ok_max - Ok_min

    fig, axes = plt.subplots(3 + calculate_chi, 1, sharex=True, figsize=(10, 7))
    
    DH_list = []
    DM_list = []

    if 'changing_Om' in str(Path.cwd().parent):
        #xlabel = f'$\left[ \Omega_k\\right]^{{{fid_tag}}}$'
        xlabel = f'$\left[ \Omega_k\\right]^{{{fid_tag}}}_{{\Omega_\Lambda = 0.69}}$'
        r_d = calculate_rd(Om_fid_cont)
        axes[1].plot(Ok_cont, DH_fid(zmax, Om_ref_cont, OL_flat)/r_d, color=color, linewidth=linewidth, label=r'$D_H^c/r_d^t$') 
        axes[1].plot(Ok_cont, DM_fid(zmax, Om_ref_cont, OL_flat)/r_d, color='coral', linewidth=linewidth, label=r'$D_M^c/r_d^t$') 
    elif 'changing_OL' in str(Path.cwd().parent):
        #xlabel = f'$\left[ \Omega_k\\right]^{{{fid_tag}}}$'
        xlabel = f'$\left[ \Omega_k\\right]^{{{fid_tag}}}_{{\Omega_m = 0.31}}$'
        r_d = calculate_rd(Om_fid_cont)
        axes[1].plot(Ok_cont, DH_fid(zmax, Om_flat, OL_ref_cont)/r_d, color=color, linewidth=linewidth, label=r'$D_H^c/r_d^t$') 
        axes[1].plot(Ok_cont, DM_fid(zmax, Om_flat, OL_ref_cont)/r_d, color='coral', linewidth=linewidth, label=r'$D_M^c/r_d^t$') 
    else: 
        print("Warning: This directory belongs to no fixed Om or OL!\n\tNo fid plot is generated.")

    try: 
        axes[3].set_xlabel(xlabel)
        axes[3].plot(Ok_cont, 0*Ok_cont + 87,  color="#000080")
    except IndexError:
        pass

    for Om, OL, Ok, apara, aperp, chi in zip(Om_list, OL_list, Ok_list, a_para, a_perp, chi_list):
        Om_ref, Om_fid = Om
        OL_ref, OL_fid = OL
        r_d = calculate_rd(Om_fid)
        current_DH = DH_fid(zmax, Om_ref, OL_ref) * apara[0]/r_d
        current_DH_std = DH_fid(zmax, Om_ref, OL_ref) * apara[1]/r_d
        current_DM = DM_fid(zmax, Om_ref, OL_ref) * aperp[0]/r_d
        current_DM_std = DM_fid(zmax, Om_ref, OL_ref) * aperp[1]/r_d
        offset = 0.001
        axes[0].errorbar(Ok, apara[0], yerr=apara[1], fmt='x', #Can be done better by just changing plotrc
                   elinewidth=elinewidth, capsize=capsize, capthick=capthick,
                  color=color, linewidth=linewidth, markersize=markersize)
        axes[0].errorbar(Ok + offset, aperp[0], yerr=aperp[1], fmt='x',
                   elinewidth=elinewidth, capsize=capsize, capthick=capthick,
                  color='coral', linewidth=linewidth, markersize=markersize)
        axes[2].errorbar(Ok, current_DH,  yerr=current_DH_std, fmt='x',
                   elinewidth=elinewidth, capsize=capsize, capthick=capthick, 
                  color=color, linewidth=linewidth, markersize=markersize)
        axes[2].errorbar(Ok + offset, current_DM,  yerr=current_DM_std, fmt='x',
             elinewidth=elinewidth, capsize=capsize, capthick=capthick, 
            color='coral', linewidth=linewidth, markersize=markersize)
        try:
            axes[3].plot(Ok, chi, 'o', color="#000080")
        except UnboundLocalError:
            pass
            
            
        DH_list.append((current_DH, current_DH_std))
        DM_list.append((current_DM, current_DM_std))
    
    axes[0].errorbar(Ok, apara[0], yerr=apara[1], fmt='x', #Can be done better by just changing plotrc
               elinewidth=elinewidth, capsize=capsize, capthick=capthick,
              color=color, linewidth=linewidth, markersize=markersize, label=r'$\alpha_{\parallel}$')
    axes[0].errorbar(Ok, aperp[0], yerr=aperp[1], fmt='x',
               elinewidth=elinewidth, capsize=capsize, capthick=capthick,
              color='coral', linewidth=linewidth, markersize=markersize, label=r'$\alpha_{\perp}$')
    axes[2].errorbar(Ok, current_DH,  yerr=current_DH_std, fmt='x',
               elinewidth=elinewidth, capsize=capsize, capthick=capthick, 
              color=color, linewidth=linewidth, markersize=markersize, label=r'$D_H/r_d$')
    axes[2].errorbar(Ok, current_DM,  yerr=current_DM_std, fmt='x',
         elinewidth=elinewidth, capsize=capsize, capthick=capthick, 
        color='coral', linewidth=linewidth, markersize=markersize, label=r'$D_M/r_d$')


    axes[0].set_ylabel(r'$\alpha_{\parallel}$', fontsize=fontsize)
    axes[0].set_ylabel(r'$\alpha_{\perp}$', fontsize=fontsize)
    axes[1].set_ylabel(r'$D_H^c/r_d^t$', fontsize=fontsize)
    axes[1].set_ylabel(r'$D_M^c/r_d^t$', fontsize=fontsize)
    axes[2].set_ylabel(r'$D_H/r_d$', fontsize=fontsize)
    axes[2].set_ylabel(r'$D_M/r_d$', fontsize=fontsize)
    axes[2].set_xlabel(xlabel, fontsize=fontsize)
    axes[2].set_xlabel(xlabel, fontsize=fontsize)

    if table:
        dic = {"Ok": Ok_list, "Om": Om_list[:, phase-2], "OL": OL_list[:, phase-2], 
                "apara": [a[0] for a in a_para],
                "apara_err": [a[1] for a in a_para],
                "aperp": [a[0] for a in a_perp],
                "aperp_err": [a[1] for a in a_perp],
                "chi2": [x[0] for x in chi_list]
                }
        df = pd.DataFrame(dic)
        tablename = f"~/TFG/results_fig/{files[0]}"
        print(f"Saving dataFrame to {tablename}")
        df.to_csv(tablename, header=True, index=False, sep='\t')
    for ax in axes[:-2]:
        ax.legend(loc='best')
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(flip(handles, 2), flip(labels, 2), loc=9, ncol=2)
    plt.tight_layout()
    if save: 
        print(f'Saving figure to {fig_name}...')
        plt.savefig(fig_name)
    if view: 
        plt.show()


def remove_bao(k_in, pk_in):
    # De-wiggling routine by Mario Ballardini

    # This k range has to contain the BAO features:
    k_ref=[2.8e-2, 4.5e-1]

    # Get interpolating function for input P(k) in log-log space:
    _interp_pk = sp.interpolate.interp1d( np.log(k_in), np.log(pk_in),
                                             kind='quadratic', bounds_error=False )
    interp_pk = lambda x: np.exp(_interp_pk(np.log(x)))

    # Spline all (log-log) points outside k_ref range:
    idxs = np.where(np.logical_or(k_in <= k_ref[0], k_in >= k_ref[1]))[0]
    _pk_smooth = sp.interpolate.UnivariateSpline( np.log(k_in[idxs]),
                                                     np.log(pk_in[idxs]), k=3, s=0 )
    pk_smooth = lambda x: np.exp(_pk_smooth(np.log(x)))

    # Find second derivative of each spline:
    fwiggle = sp.interpolate.UnivariateSpline(k_in, pk_in / pk_smooth(k_in), k=3, s=0)
    derivs = np.array([fwiggle.derivatives(_k) for _k in k_in]).T
    d2 = sp.interpolate.UnivariateSpline(k_in, derivs[2], k=3, s=1.0)

    # Find maxima and minima of the gradient (zeros of 2nd deriv.), then put a
    # low-order spline through zeros to subtract smooth trend from wiggles fn.
    wzeros = d2.roots()
    wzeros = wzeros[np.where(np.logical_and(wzeros >= k_ref[0], wzeros <= k_ref[1]))]
    wzeros = np.concatenate((wzeros, [k_ref[1],]))
    wtrend = sp.interpolate.UnivariateSpline(wzeros, fwiggle(wzeros), k=3, s=0)

    # Construct smooth no-BAO:
    idxs = np.where(np.logical_and(k_in > k_ref[0], k_in < k_ref[1]))[0]
    pk_nobao = pk_smooth(k_in)
    pk_nobao[idxs] *= wtrend(k_in[idxs])

    # Construct interpolating functions:
    ipk = sp.interpolate.interp1d( k_in, pk_nobao, kind='linear',
                                      bounds_error=False, fill_value=0. )

    pk_nobao = ipk(k_in)

    return pk_nobao

def calculate_olin(k_in, pk_in):
    """
    Parameters
        a pandas dataframe with 2 columns, the k_in and pk_in
    Output
        a pandas dataframe with 2 columns, the k_in and the pk_in/pk_nobao
    """
    pk_nobao = remove_bao(k_in, pk_in)
    olin = pk_in/pk_nobao
    return olin

def interpolate(k_in, pk_in):
    pk_interpolate = sp.interpolate.interp1d(k_in, pk_in, kind='linear',
                                            bounds_error=False, fill_value=0.)
    return pk_interpolate

def split_array(arr, n): 
    # Calculate length of subarray 
    length = len(arr)//n 
  
    # Create list for 2D array 
    splitted_arr = [] 
  
    # Split array into n subarrays 
    for i in range(n): 
        # Create subarray 
        subarr = arr[i*length:(i+1)*length] 
  
        # Append subarray to list 
        splitted_arr.append(subarr) 
  
    # Return list 
    return splitted_arr 

rustico_path = list(Path('/home/santi/TFG/DATA/rustico_output').glob('Power_*'))
class_output = list(Path('/home/santi/TFG/class_public/output').glob('*.dat'))
model = list(Path('/home/santi/TFG/lrg_eboss/model/').glob('*'))
#mcmc_output = list(Path('/home/santi/TFG/lrg_eboss/output').glob('mcmc*.txt'))
mcmc_output = []
for file in Path('/home/santi/TFG/lrg_eboss/output').glob('mcmc*.txt'):
    if '__' in str(file):
        continue
    mcmc_output.append(file)

if __name__ == "__main__":
    print("You are not running the file you should be running!")

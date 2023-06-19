import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import scipy as sp
import scipy.constants as ct
import re
from scipy.stats import norm

zmax = 0.698
H0 = 67.6
rs = 147.784

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
        print('But not here')

def alpha_from_logfile(filename):
    """
    parameters:
        the pathfix of the logfile to be calculated
    returns 
        (Ok, (apara, aparasigma), (aperp, aperpsigma))
    """
    p = re.compile('logfile')
    namestr = filename.stem
    find = p.findall(namestr)
    if len(find) == 0:
        raise TypeError(f"{filename.stem} is not a logfile!")
    omegas = get_params(namestr)
    out = []
    Ok = omegas[-1]
    out.append(Ok)
    with filename.open() as open_data:
        for i, line in enumerate(open_data):
            if i == 9:
                p = re.compile('[0-9].[0-9]*e\+?-?[0-9]*')
                matches = p.findall(line)
                out.append([float(x) for x in matches])
            elif i == 10:
                p = re.compile('[0-9].[0-9]*e\+?-?[0-9]*')
                matches = p.findall(line)
                out.append([float(x) for x in matches])
            elif i>10: break
    return out

def iterable_output(func):
    def wrapper(zmax, Ok):
        try:
            iter(Ok)
            result = []
            for ok in Ok:
                result.append(func(zmax, ok))
            return np.array(result)
        except TypeError:
            return func(zmax, Ok)
    return wrapper

def get_params(param_name):
    """
    input:
        param_name: the parameters in the name of each dataframe
    output:
        [Om, OL]: a tuple containing the parameters in the name
        """
    p = re.compile('O[a-zA-Z][0-9]*')
    params = re.findall(p, param_name)
    params = re.split('\D', "".join(params))
    OmOL = []
    for i in params:
        try:
            OmOL.append(int(i)/100)
        except ValueError:
            continue
    OmOL.append(round(1 - sum(OmOL), 2))
    return OmOL 

def plot_DH_DM(files, fig_name=None, save=False,  view=False, n_points = 500, markersize = 10, elinewidth=3, capsize=5, capthick=3, fontsize = 20, linewidth = 3, color = 'teal', reduce_ticks = 2):
    Ok_list = []
    a_para = []
    a_perp = []

    for data in files:
        print('Opening file', str(data), '...')
        out = alpha_from_logfile(data)
        print(out, '')
        Ok_list.append(out[0])
        a_para.append(out[1])
        a_perp.append(out[2])
        
    print(Ok_list, a_para, a_perp, sep='\n')

    Ok_min, Ok_max = min(Ok_list), max(Ok_list)
    Ok_cont = np.linspace(Ok_min, Ok_max, n_points)
    Ok_rang = Ok_max - Ok_min
    H = lambda z, Ok, Om: H0*np.sqrt(Om*(1+z)**3 + Ok*(1+z)**2 + 1-Ok-Om)


    @iterable_output
    def DH_fid(z, Ok):
        return ct.c/1000/H(z, Ok, Om)
    @iterable_output
    def DC_fid(z, Ok, Om):
        return sp.integrate.quad(DH_fid, 0, z, args=(Ok, Om))[0]
    @iterable_output
    def DM_fid(z, Ok, Om):
        DC = DC_fid(z, Ok, Om)
        DH = DH_fid(z, Ok, Om)
        if Ok>0:
            k =  DH/np.sqrt(Ok)
            return k*np.sinh(np.sqrt(Ok)*DC/DH)
        elif Ok<0:
            k =  DH/np.sqrt(np.abs(Ok))
            return k*np.sin(np.sqrt(np.abs(Ok))*DC/DH)
        elif not Ok:
            return DC

    #---Changing z->d (phase 2)
    #DM_fid = np.array([DM(zmax, Ok, DC) for Ok, DC in zip(Ok_cont, DC_fid)])
    #---Changing template (phase 3)

    fig, axes = plt.subplots(3, 2, sharex=True, figsize=(10, 7))
    
    DH_list = []
    DM_list = []

    for Ok, apara, aperp in zip(Ok_list, a_para, a_perp):
        current_DH = DH_fid(zmax, Ok, Om) * apara[0]/rs
        current_DH_std = DH_fid(zmax, Ok, Om) * apara[1]/rs
        current_DM = DM_fid(zmax, Ok, Om) * aperp[0]/rs
        current_DM_std = DM_fid(zmax, Ok) * aperp[1]/rs
        axes[0,0].errorbar(Ok, apara[0], yerr=apara[1], fmt='x',
                     elinewidth=elinewidth, capsize=capsize, capthick=capthick,
                    color=color, linewidth=linewidth, markersize=markersize)
        axes[0,1].errorbar(Ok, aperp[0], yerr=aperp[1], fmt='x',
                     elinewidth=elinewidth, capsize=capsize, capthick=capthick,
                    color=color, linewidth=linewidth, markersize=markersize)
        axes[2,0].errorbar(Ok, current_DH,  yerr=current_DH_std, fmt='x',
                     elinewidth=elinewidth, capsize=capsize, capthick=capthick, 
                    color=color, linewidth=linewidth, markersize=markersize)
        axes[2,1].errorbar(Ok, current_DM,  yerr=current_DM_std, fmt='x',
             elinewidth=elinewidth, capsize=capsize, capthick=capthick, 
            color=color, linewidth=linewidth, markersize=markersize)
        DH_list.append((current_DH, current_DH_std))
        DM_list.append((current_DM, current_DM_std))
       # print("""For Ok = {} 
       #       D_H/r_d = {} \\pm {} 
       #       D_M/r_d = {} \\pm {} """.format(Ok,*[round(x, 2) for x in (current_DH,current_DH_std,current_DM,current_DM_std)]))
        

    axes[1,0].plot(Ok_cont, DH_fid(zmax, Ok_cont)/rs, color=color, linewidth=linewidth) #Multiply by 0 is phase2
    axes[1,1].plot(Ok_cont, DM(zmax, Ok_cont)/rs, color=color, linewidth=linewidth)
    axes[0,0].set_ylabel(r'$\alpha_{\parallel}$', fontsize=fontsize)
    axes[0,1].set_ylabel(r'$\alpha_{\perp}$', fontsize=fontsize)
    axes[1,0].set_ylabel(r'$\left[ D_H/r_d\right]^{fid}$', fontsize=fontsize)
    axes[1,1].set_ylabel(r'$\left[ D_M/r_d\right]^{fid}$', fontsize=fontsize)
    axes[2,0].set_ylabel(r'$D_H/r_d$', fontsize=fontsize)
    axes[2,1].set_ylabel(r'$D_M/r_d$', fontsize=fontsize)
    axes[2,0].set_xlabel(r'$\left[ \Omega_k\right]^{fid}$', fontsize=fontsize)
    axes[2,1].set_xlabel(r'$\left[ \Omega_k\right]^{fid}$', fontsize=fontsize)
    for ax in axes.ravel():
        xticks = [i/100 for i in range(-20, 21, 10)]
        xlabel = xticks.copy()
        ax.set_xticks(xticks)
        ax.set_xticklabels([round(tick, 2) for tick in xticks], fontsize=fontsize)
        yticks = ax.get_yticks()[::reduce_ticks]
        ylabel = ax.get_ylabel()[::reduce_ticks]
        ax.set_yticks(yticks)
        ax.set_yticklabels([round(tick, 2) for tick in yticks], fontsize=fontsize)

    plt.tight_layout()
    if save: 
        plt.savefig(fig_name)
    if view: 
        plt.show()


def calculate_avg_and_std(data):
    """
    input: data should be ((value1, standard_deviation1), (value2, standard_deviation2), ...)
    output: sum(data_i/std_i**2)/sum(1/std_i**2)
    """
    total_avg = 0
    total_std = 0
    for value, std in data:
        total_avg += value/std**2
        total_std += 1/std**2
    return round(total_avg/total_std, 2), round(1/total_std, 2)

def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return (average, np.sqrt(variance))

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
def remove_bao_substraction(k_in, pk_in):
    """
    Parameters:
        k_in array
        pk_in array, of same length
    Returns:
        pk_nobao: an array of equal length without the BAO.
    """
    # This is copied from the code in https://github.com/brinckmann/montepython_public
    # De-wiggling routine by Mario Ballardini
    # This k range has to contain the BAO features:
    # changed it to regular substraction.
    k_ref=[2.8e-2, 4.5e-1]

    # Get interpolating function for input P(k) in log-log space:
    _interp_pk = sp.interpolate.interp1d(np.log(k_in), np.log(pk_in),
                                             kind='quadratic', bounds_error=False)
    interp_pk = lambda x: np.exp(_interp_pk(np.log(x)))

    # Spline all (log-log) points outside k_ref range:
    idxs = np.where(np.logical_or(k_in <= k_ref[0], k_in >= k_ref[1]))
    _pk_smooth = sp.interpolate.UnivariateSpline(np.log(k_in[idxs]),
                                                     np.log(pk_in[idxs]), k=3, s=0)
    pk_smooth = lambda x: _pk_smooth(np.log(x)) #Used to be an np.exp around everything
                                                #without it, looks more promising

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
    idxs = np.where(np.logical_and(k_in > k_ref[0], k_in < k_ref[1]))
    pk_nobao = pk_smooth(k_in)
    pk_nobao[idxs] *= wtrend(k_in[idxs])

    # Construct interpolating functions:
    ipk = sp.interpolate.interp1d(k_in, pk_nobao, kind='linear',
                                      bounds_error=False, fill_value=0.)

    pk_nobao = ipk(k_in)
    return pk_nobao

def fft_bao(k_in, pk_in):
    """Returns the k²P(k)sin(kx)/(kx) integral for given k and P(k) functions
    input: k_in (array), pk_in(array)
    returns: (R, Xi), both np.arrays. Xi is the integral solved for every r
    """
    #Done on a discrete array
    R = np.arange(len(k_in))
    Xi = []
    for r in R:
        Xi.append(sp.integrate.trapz(k_in * pk_in * np.sin(k_in*r)/r, k_in))
    return R, np.array(Xi)

def interpolate_fft(k_in, pk_in, R, kind='linear', n_points=10000):
    """
    Parameters:
        k_in(array): input of k values
        pk_in(array): input of P(k) values
        R: array to map the R values to
    Returns:
        integrate: An integral of k²P(k)sin(kx)/kx from 0 to inf,
                    interpolating P(k) to a continuous function.
    """
    kmin, kmax = k_in[0], k_in[-1]
    pk_interpolate = sp.interpolate.interp1d(k_in, pk_in, kind='linear',
                                            bounds_error=False, fill_value=0.)
    #k_interpolate = sp.interpolate.interp1d(k_in, k_in, kind='linear',
#                                            bounds_error=False, fill_value=0.)
    k_interpolate = np.linspace(kmin, kmax, n_points)
    Xi = []
    for r in R:
       Xi.append(sp.integrate.trapz(k_interpolate * pk_interpolate(k_interpolate) * np.sin(k_interpolate * r)/r, k_interpolate))
    return R, Xi

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

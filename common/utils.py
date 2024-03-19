import os
import numpy as np

def get_nprocs():
    # get the number of physical cores per socket
    nprocs = int(os.popen("lscpu | grep 'Core(s) per socket' | awk '{print $4}'").read().strip())
    # get the number of sockets
    nsockets = int(os.popen("lscpu | grep 'Socket(s)' | awk '{print $2}'").read().strip())
    nprocs *= nsockets

    # we need to limit the number of values if nprocs is too large
    if nprocs > 32:
        # at least the power of 2 values
        llist = [2**i for i in range(0, np.log2(nprocs).astype(int) + 1)] + [nprocs]
        # add the multiples of 8 up to nprocs, that are not in the list
        llist.extend([i for i in range(8, nprocs, 8) if i not in llist])
        # the values must be unique and sorted
        return np.unique(llist)
    else:
        # list of number of threads to use: 1, 2, ... nprocs, with step 2
        return np.insert(np.arange(2, nprocs+1, 2), 0, 1)
    
def get_best_level(results, cloud, nth, key='total'):
    """
        Get the best level for the given parameters
    """
    # get the best level or chunk for the given parameters
    return min(results[cloud], key=lambda x: results[cloud][x][nth][key])

def get_best_level2(ires, nth, pos=2):
    """
        Get the best level for the given parameters
    """
    # get the best level or chunk for the given parameters
    return min(ires, key=lambda x: ires[x][nth][pos])

def get_best_comb(ires, num_threads, pos=2):
    """
        Get the best level, num processes and time for the given parameters
    """
    # create a flattened list of all combinations of x and y along with the time
    flattened = [(x, y, ires[x][y][pos]) for x in ires for y in num_threads]
    # find the minimum based on Total time; this return bestlevel, bestthreads, time
    return min(flattened, key=lambda item: item[2])

def get_best_comb_3d(ires, pos=2):
    """
        Get the best configuration for the given parameters
    """
    # results[name][i][minrad][16][phase]
    xdim = list(ires.keys())
    ydim = list(ires[xdim[0]].keys())
    zdim = list(ires[xdim[0]][ydim[0]].keys())
    # create a flattened list of all combinations of x and y along with the time
    flattened = [(x, y, z, ires[x][y][z][pos]) for x in xdim for y in ydim for z in zdim]
    # find the minimum based on Total time
    return min(flattened, key=lambda item: item[-1])

def get_best_optimization(df, key='Total'):
    """This function returns the best optimization for the given dataframe.

    Parameters
    ----------
    df
        dataframe with the results obatined from All_Optimizations-<hostname>.csv. See process_all_optim.ipynb for more details.

    Returns
    -------
        best: df with the best optimization
        best_label: label of the best optimization
    """
    # get all the different optimizations
    optimizations = df['Optimization'].unique()
    # select best initial case
    best = df.loc[df['Optimization'] == optimizations[0], 'TimeTree':'Total'].copy()
    # get the best time for the initial case
    best_time = best[key].mean()
    # keep the label of the initial case
    best_label = optimizations[0]
    # iterate over all the optimizations
    for opt in optimizations[1:]:
        # select the optimization
        opt_df = df.loc[df['Optimization'] == opt, 'TimeTree':'Total'].copy()
        # get the time for the optimization
        opt_time = opt_df[key].mean()
        # if the time is better than the best time, update
        if opt_time < best_time:
            best = opt_df
            best_time = opt_time
            best_label = opt
    
    return best, best_label

def tokenize_sycl_cpu(filename):
    experiment ={}

    # the executable is the key to continue
    exe = "empty"

    with open(filename) as f:
        for line in f:
            tokens = line.split()
            if "Running:" in tokens:
                # this line format is: Running: ../bin/owm-sycl-cpu ../bin/data/INAER_2011_Alcoy_Core.xyz 4 8
                # tokens position: 1:executable, 2:input, 3:maxNum, 4:threads
                exe=tokens[1].split("/")[2]
                # continue if the executable is not CPU
                if "cpu" not in exe:
                    continue
                input=tokens[2].split("/")[3][:-4]
                cloud = {
                    "INAER_2011_Alcoy_Core": "Alcoy",
                    "BABCOCK_2017_Arzua_3B": "Arzua",
                    "V21_group1_densified_point_cloud": "BrionF",
                    "V19_group1_densified_point_cloud": "BrionU",
                }[input]
                maxNum=int(tokens[3])
                # get the number of threads
                nth = int(tokens[4])
                if exe not in experiment:
                    experiment[exe]={}
                if cloud not in experiment[exe]:
                    experiment[exe][cloud]={}
                if maxNum not in experiment[exe][cloud]:
                    experiment[exe][cloud][maxNum]={}
                experiment[exe][cloud][maxNum][nth]=[]
            # continue if the executable is not CPU
            if "cpu" not in exe:
                continue
            if 'Construction' in tokens:
                experiment[exe][cloud][maxNum][nth].append(float(tokens[5]))
            if "Stage1" in tokens:
                experiment[exe][cloud][maxNum][nth].append(float(tokens[5]))
            if "Stage2" in tokens:
                experiment[exe][cloud][maxNum][nth].append(float(tokens[5]))
            if "Total" in tokens and "KERNEL" in tokens:
                experiment[exe][cloud][maxNum][nth].append(float(tokens[5]))
            if 'TIME' in tokens:
                experiment[exe][cloud][maxNum][nth].append(float(tokens[6]))

    return experiment

def get_nested_values(mydict, pos=3):
    """Return the nested values of a dictionary
    """
    if pos == 0:
        return mydict
    return get_nested_values(list(mydict.values())[0], pos-1)
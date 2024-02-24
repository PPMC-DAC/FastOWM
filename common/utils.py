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
        llist = [2**i for i in range(1, np.log2(nprocs).astype(int) + 1)] + [nprocs]
        # add the multiples of 8 up to nprocs, that are not in the list
        llist.extend([i for i in range(8, nprocs, 8) if i not in llist])
        # the values must be unique and sorted
        return np.unique(llist)
    else:
        # list of number of threads to use: 1, 2, ... nprocs, with step 2
        return np.insert(np.arange(2, nprocs+1, 2), 0, 1)
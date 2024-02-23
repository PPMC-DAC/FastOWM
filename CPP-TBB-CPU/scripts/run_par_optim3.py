import os
import time
import numpy as np
from datetime import datetime

# get the number of physical cores
nprocs = int(os.popen("lscpu | grep 'Core(s) per socket' | awk '{print $4}'").read().strip())

#Sliding window size
Wsize = 10
#Gris size for stage 3 (filling empty areas of Bsize*Bsize with a minimum point)
Bsize = 20
#Overlap for the sliding window. Displacement will be Wsize(1-Overlap)=10m*0.2=2m
Overlap = 0.8
#num_threads for the openmp implementations of stage1 and stage3
num_threads = np.insert(np.linspace(2, nprocs, nprocs//2, dtype=int, endpoint=True), 0, 1)
#number of times the OWM is executed
nreps = 5

executable_list=['../bin/o3memoA', '../bin/o3memo', '../bin/o3memoB']
inputs=[
    "../bin/data/AlcoyH",
    "../bin/data/ArzuaH",
    "../bin/data/BrionFH",
    "../bin/data/BrionUH",
    ]

# get the hostname
hostname = os.popen("hostname").read().strip()
# set the output file
output_list = [f'o3_memoizationA_{hostname}.out', f'o3_memoization_{hostname}.out', f'o3_memoizationB_{hostname}.out']

levels = list(range(2,10))
mR=0.1
maxNumber=32

# zip the output and executable lists
for output, executable_par in zip(output_list, executable_list):

    start = time.time()
    print("Start : %s" % time.ctime())

    with open(output, "a") as f:
        # stamp the start time
        f.write(f'Start: {datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}')
        # iterate over the clouds, levels and number of threads
        for cloud in inputs:
            for lev in levels:
                for nth in num_threads:
                    print("Running: {} -i {} -W {} -B {} -O {} -n {} -l {} -r {} -s {} -L {}".format(executable_par,cloud,Wsize,Bsize,Overlap,nth,nreps,mR,maxNumber,lev))
                    # save the configuration in the file
                    f.write("\n\nRunning: {} -i {} -W {} -B {} -O {} -n {} -l {} -r {} -s {} -L {}\n\n".format(executable_par,cloud,Wsize,Bsize,Overlap,nth,nreps,mR,maxNumber,lev))
                    # flush the buffer
                    f.flush()
                    # excecute the command and save the output to the file
                    os.system("%s -i %s -W %d -B %d -O %f -n %d -l %d -r %f -s %d -L %d| tee -a %s" % (executable_par, cloud, Wsize, Bsize, Overlap, nth, nreps,mR, maxNumber, lev, output))

        end = time.time()
        f.write(f'End: {datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}\n')
        f.write("Total Execution time: {} hours\n".format((end - start)/3600))

    print("End : %s" % time.ctime())
    print("Total Execution time: %f hours" % ((end - start)/3600))
    # copy the output file to the results folder
    os.system(f"cp {output} ../../Results/{output.replace('.out', '.txt')}")